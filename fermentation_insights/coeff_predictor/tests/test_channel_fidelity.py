# -*- coding: utf-8 -*-
"""Round 5 (channel / price-intensity fidelity): the channel-only surface
loss (LossConfig.ch_space), the power gate (ModelConfig.power_gate with the
structural input struct_power_export) and its data check."""
import numpy as np
import pandas as pd
import pytest
import torch

from fermentation_insights.coeff_predictor.features import (
    FeatureSpec, STRUCTURAL, POWER_STRUCTURAL, power_export_factor)
from fermentation_insights.coeff_predictor.model import ModelConfig, make_layout, mpsp_equation
from fermentation_insights.coeff_predictor import train as T
from fermentation_insights.coeff_predictor.tests._fake_tables import make_fake_tables
from fermentation_insights.coeff_predictor.tests.conftest import build


@pytest.fixture(scope='module')
def tabs():
    # LA has no electricity coproduct; on glucose its fake electricity
    # intensity is positive (no lignin)
    return make_fake_tables(n_samples_per_config=4, products=['HP', 'EtOH', 'LA'], feedstocks=['glucose'], seed=7)


def _model(tabs, seed=0, power_export=True, **kw):
    samples, points, _ = tabs
    spec = FeatureSpec(downstream=True, power_export=power_export).fit(samples)
    batch = spec.transform(samples)
    cfg = ModelConfig(**{**dict(width=16, depth=2, dropout=0., water_gate=True, channel_inputs='physics',
                                masked_width=16, r_structure='physics'), **kw})
    layout = make_layout(cfg, samples.columns)
    scales = T.fit_target_scales(samples, points, layout, batch, spec.continuous_mask(), T.fs_prior_z(spec, batch))
    torch.manual_seed(seed)
    m = T.build_model(spec, layout, scales, cfg, torch.float64).eval()
    X = torch.tensor(batch.X)
    rest = [torch.tensor(np.asarray(a, float)) for a in (batch.log_r, batch.prices, batch.w_sugar, batch.g_max)]
    return m, spec, layout, scales, X, rest, samples


def test_power_export_factor():
    df = pd.DataFrame(dict(coproducts=['electricity', 'DDGS+corn_oil', 'none', np.nan, 'DDGS+corn_oil+electricity',
                                       '']))
    f = power_export_factor(df)
    np.testing.assert_array_equal(f[[0, 1, 2, 4]], [1., 0., 0., 1.])
    # unknown (missing / empty) -> NaN, i.e. ungated after the FeatureSpec's unknown -> 1
    assert np.isnan(f[3]) and np.isnan(f[5])
    assert POWER_STRUCTURAL not in STRUCTURAL


def test_feature_spec_power_export_is_opt_in(tabs):
    samples = tabs[0]
    s0 = FeatureSpec(downstream=True).fit(samples)
    assert s0.structural == list(STRUCTURAL) and s0.structural_index(POWER_STRUCTURAL) is None
    s1 = FeatureSpec(downstream=True, power_export=True).fit(samples)
    assert s1.structural == list(STRUCTURAL) + [POWER_STRUCTURAL]
    assert s1.columns[:-1] == s0.columns           # appended last, the rest unchanged
    X0, X1 = s0.transform(samples).X, s1.transform(samples).X
    np.testing.assert_array_equal(X1[:, :-1], X0)
    la = samples['fermentation_product'].eq('LA').to_numpy()
    np.testing.assert_array_equal(X1[:, -1], np.where(la, 0., 1.))
    # never an encoder input, in no feature group
    assert s1.n_features - 1 not in s1.encoder_index()
    assert all(s1.n_features - 1 not in v for v in s1.column_groups().values())
    r = FeatureSpec.from_dict(s1.to_dict())
    assert r.power_export and r.structural == s1.structural
    r0 = FeatureSpec.from_dict(s0.to_dict())
    assert not r0.power_export and r0.structural == list(STRUCTURAL)
    assert not FeatureSpec(power_export=True).power_export      # needs downstream


def test_power_gate(tabs):
    m, spec, layout, scales, X, rest, samples = _model(tabs, power_gate=True)
    assert m.power_index == spec.structural_index(POWER_STRUCTURAL)
    j = layout['price_keys'].index('el')
    gate = X[:, m.power_index].numpy()
    la = samples['fermentation_product'].eq('LA').to_numpy().copy()
    assert (gate[la] == 0).all() and (gate[~la] == 1).all()
    # an ungated twin (same seed): identical parameters, same outputs except
    # the electricity intensity of the gated rows
    m0 = _model(tabs, power_gate=False)[0]
    assert m0.power_index is None
    with torch.no_grad():
        # random weights large enough to give negative raw intensities
        for p, p0 in zip(m.parameters(), m0.parameters()):
            p.normal_(0., 1.)
            p0.copy_(p)
        o, o0 = m(X, *rest), m0(X, *rest)
    q, q0 = o['intensities'][:, j], o0['intensities'][:, j]
    assert (q0[la] < 0).any()                      # the twin would export power
    assert torch.all(q[la] >= 0)
    torch.testing.assert_close(q[~la], q0[~la])
    other = [i for i in range(len(layout['price_keys'])) if i != j]
    torch.testing.assert_close(o['intensities'][:, other], o0['intensities'][:, other])
    torch.testing.assert_close(o['channels'], o0['channels'])
    torch.testing.assert_close(o['g'], o0['g'])
    # the structural input reaches nothing else
    X2 = X.clone()
    X2[:, m.power_index] = 1.
    with torch.no_grad():
        o2 = m(X2, *rest)
    torch.testing.assert_close(o2['intensities'][:, j], q0)
    # exact price linearity and dMPSP/dr <= 0 survive the gate (r-structured electricity)
    assert 'electricity' in m.rs_names
    log_r, prices, w, gm = rest
    with torch.no_grad():
        p2 = prices.clone()
        p2[:, 1] *= 2.
        o3 = m(X, log_r, p2, w, gm)
        torch.testing.assert_close(o3['price_channels'][:, j], 2 * o['price_channels'][:, j])
        qs = [m(X, log_r + d, prices, w, gm)['intensities'][:, j] for d in (-1., 0., 1.)]
    assert torch.all(qs[0] >= qs[1] - 1e-12) and torch.all(qs[1] >= qs[2] - 1e-12)


def test_fs_b_only(tabs):
    """The pure-b' feedstock structure: a' = c' = d' = 0 exactly, b' as in the
    twin without it (same parameters), every other output unchanged."""
    m, spec, layout, scales, X, rest, samples = _model(tabs, fs_b_only=True, fs_prior=True)
    m0 = _model(tabs, fs_b_only=False, fs_prior=True)[0]
    with torch.no_grad():
        for p, p0 in zip(m.parameters(), m0.parameters()):
            p.normal_(0., 1.)
            p0.copy_(p)
        o, o0 = m(X, *rest), m0(X, *rest)
    j = layout['price_keys'].index('fs')
    q, q0 = o['intensities'][:, j], o0['intensities'][:, j]
    assert torch.all(q[:, [0, 2, 3]] == 0) and torch.any(q0[:, [0, 2, 3]] != 0)
    torch.testing.assert_close(q[:, 1], q0[:, 1])
    other = [i for i in range(len(layout['price_keys'])) if i != j]
    torch.testing.assert_close(o['intensities'][:, other], o0['intensities'][:, other])
    torch.testing.assert_close(o['channels'], o0['channels'])
    assert ModelConfig().fs_b_only is False
    d = ModelConfig().to_dict()
    d.pop('fs_b_only')
    assert ModelConfig.from_dict(d).fs_b_only is False


def test_clip_exempt(tabs):
    """clip_exempt keeps the listed columns unclipped (exact structural
    quantities), every other continuous column keeps its clip bounds."""
    col = 'final_product_stoich_mass_ratio'
    m, spec, layout, scales, X, rest, samples = _model(tabs, clip_inputs=0.0, clip_exempt=[col])
    m0 = _model(tabs, clip_inputs=0.0)[0]
    j = spec.columns.index(col)
    assert float(m.x_lo[j]) <= -1e29 and float(m.x_hi[j]) >= 1e29
    assert float(m0.x_hi[j]) < 1e29
    other = [i for i in range(spec.n_features) if i != j]
    torch.testing.assert_close(m.x_lo[other], m0.x_lo[other])
    torch.testing.assert_close(m.x_hi[other], m0.x_hi[other])
    assert ModelConfig(clip_exempt=col).clip_exempt == [col] and ModelConfig(clip_exempt=[]).clip_exempt is None
    d = ModelConfig().to_dict()
    d.pop('clip_exempt')
    assert ModelConfig.from_dict(d).clip_exempt is None


def test_power_gate_data_check(tabs):
    samples, points, _ = tabs
    la = samples['fermentation_product'].eq('LA').to_numpy()
    rep = T.check_constraints(samples, points).set_index('constraint')
    assert rep.loc['power_gate', 'enabled'] and rep.loc['power_gate', 'n_cases'] == la.sum()
    bad = samples.copy()
    bad.loc[la, 'int__el__a'] = -5.
    rep2 = T.check_constraints(bad, points)
    assert not rep2.set_index('constraint').loc['power_gate', 'enabled']
    mc, _ = T.apply_constraint_report(rep2, ModelConfig(power_gate=True), T.LossConfig())
    assert mc.power_gate is False
    # per-point simulated intensities win over the fitted surface
    pts = points.copy()
    ids = set(samples.loc[la, 'sample_id'].astype(str))
    pts['q__el'] = 1.
    pts.loc[pts['sample_id'].astype(str).isin(ids), 'q__el'] = -1.
    rep3 = T.check_constraints(samples, pts).set_index('constraint')
    assert not rep3.loc['power_gate', 'enabled']


def test_power_gate_config_defaults():
    assert ModelConfig().power_gate is False
    assert ModelConfig(kind='plain', power_gate=True).power_gate is False
    d = ModelConfig().to_dict()
    d.pop('power_gate')
    assert ModelConfig.from_dict(d).power_gate is False


def test_power_gate_ensemble_roundtrip(tabs, tmp_path):
    samples, points, _ = tabs
    cfg = ModelConfig(width=16, depth=2, dropout=0., water_gate=True, channel_inputs='physics', masked_width=16,
                      power_gate=True)
    ens = T.train_ensemble(samples, points, 1, cfg, T.LossConfig(),
                           T.TrainConfig(epochs=3, batch_size=16, min_epochs=1, patience=2),
                           calib_frac=0., spec_kwargs=dict(downstream=True))
    assert ens.spec.power_export and ens.mcfg.power_gate and ens.members[0].power_index is not None
    _, o1 = ens.predict_samples(samples)
    la = samples['fermentation_product'].eq('LA').to_numpy()
    j = ens.layout['price_keys'].index('el')
    assert (o1['intensities'][:, la, j] >= 0).all()          # (members, n, J, 4)
    ens2 = T.Ensemble.load(ens.save(str(tmp_path / 'ck')))
    _, o2 = ens2.predict_samples(samples)
    np.testing.assert_allclose(o2['coef'], o1['coef'], rtol=1e-6, atol=1e-8)
    assert ens2.mcfg.power_gate and ens2.spec.power_export
    # without the downstream FeatureSpec the gate has no input and is inert
    ens3 = T.train_ensemble(samples, points, 1, cfg, T.LossConfig(),
                            T.TrainConfig(epochs=1, batch_size=16, min_epochs=1, patience=1), calib_frac=0.)
    assert ens3.members[0].power_index is None and not ens3.spec.power_export


def _batch(fitted, layout, dtype=torch.float64):
    s, p = fitted['samples'], fitted['points']
    return T.make_tensors(fitted['spec'], s, p, layout, fitted['batch'], dtype=dtype)


def test_channel_surface_loss_keeps_total_in_coefficient_space(fitted):
    model, _, _, _, layout, _ = build(fitted)
    b = _batch(fitted, layout)
    with torch.no_grad():
        out = model(*[b[k] for k in ('X', 'log_r', 'prices', 'w_sugar', 'g_max')])
        Lc = T.data_loss(model, out, b, T.LossConfig())
        Ls = T.data_loss(model, out, b, T.LossConfig(coef_space='surface'))
        Lh = T.data_loss(model, out, b, T.LossConfig(ch_space='surface'))
        Lcc = T.data_loss(model, out, b, T.LossConfig(ch_space='coef'))
        Lsc = T.data_loss(model, out, b, T.LossConfig(coef_space='surface', ch_space='coef'))
    torch.testing.assert_close(Lh['tot'], Lc['tot'])
    torch.testing.assert_close(Lh['ch'], Ls['ch'])
    torch.testing.assert_close(Lh['mpsp'], Lc['mpsp'])
    for k in ('tot', 'ch', 'g', 'mpsp', 'data'):
        torch.testing.assert_close(Lcc[k], Lc[k])        # ch_space 'coef' = the default
    torch.testing.assert_close(Lsc['tot'], Ls['tot'])
    torch.testing.assert_close(Lsc['ch'], Lc['ch'])
    assert not torch.allclose(Lh['data'], Lc['data'])


def test_channel_surface_loss_is_the_fidelity_metric(fitted):
    """Perfect totals, one channel 10 % of the total numerator too high at
    every point -> ch = 0.01 (sum over channels of the mean squared error
    relative to the total numerator)."""
    model, _, _, _, layout, _ = build(fitted)
    b = _batch(fitted, layout)
    q = b['q']
    out = dict(coef=b['coef'].clone(), g=b['g'].clone(), channels=b['ch'].clone(), intensities=q.clone(),
               price_channels=b['prices'][:, model.price_index][:, :, None] * q)
    L0 = T.data_loss(model, out, b, T.LossConfig(ch_space='surface'))
    assert float(L0['ch']) < 1e-12 and float(L0['tot']) < 1e-12
    out['channels'] = out['channels'].clone()
    out['channels'][:, 0] = out['channels'][:, 0] + 0.1 * b['coef']
    L = T.data_loss(model, out, b, T.LossConfig(ch_space='surface'))
    Nt = T._surface(b['coef'][:, None], b['Y'], b['T'])[:, 0]
    big = (Nt.abs() >= 0.05 * model.scales_np['num_scale']) | (b['mask'] == 0)
    assert bool(big.all())                                # no denominator floor in the fake data
    assert abs(float(L['ch']) - 0.01) < 1e-6
    assert float(L['tot']) < 1e-12


def test_training_with_channel_surface_loss(fitted):
    from fermentation_insights.coeff_predictor.tests.conftest import FAST
    s, p = fitted['samples'], fitted['points']
    ens = T.train_ensemble(s, p, 1, lcfg=T.LossConfig(ch_space='surface', w_ch=3.),
                           tcfg=T.TrainConfig(**FAST), calib_frac=0.)
    h = ens.histories[0]
    assert h[-1]['train_loss'] < h[0]['train_loss'] and np.isfinite(h[-1]['train_ch'])
    assert T.LossConfig.from_dict(ens.lcfg.to_dict()).ch_space == 'surface'


# %% pipeline plumbing (experiments.py) and the fidelity metrics

def test_mask_channel_targets(tabs):
    samples = tabs[0]
    m = T.mask_channel_targets(samples, {'feedstock': ['LA']})
    la = samples['fermentation_product'].eq('LA').to_numpy()
    fs_cols = [c for c in samples.columns if c.startswith(('ch__feedstock__', 'int__fs__'))]
    assert len(fs_cols) == 8
    assert m.loc[la, fs_cols].isna().all().all()
    pd.testing.assert_frame_equal(m.loc[~la], samples.loc[~la])
    other = [c for c in samples.columns if c not in fs_cols]
    pd.testing.assert_frame_equal(m[other], samples[other])
    # a non-price channel: only its own targets; no matching rows: unchanged
    m2 = T.mask_channel_targets(samples, {'closure': ['HP_glucose']})
    hp = samples['config_id'].eq('HP_glucose').to_numpy()
    assert m2.loc[hp, [f'ch__closure__{k}' for k in 'abcd']].isna().all().all()
    assert not m2.loc[hp, fs_cols].isna().any().any()
    pd.testing.assert_frame_equal(T.mask_channel_targets(samples, {'feedstock': ['none']}), samples)
    # masked rows still train (totals and the other channels are kept)
    lay = make_layout(ModelConfig(), m.columns)
    _, _, ch, q = T._targets(m, lay)
    assert np.isnan(q[la, lay['price_keys'].index('fs')]).all() and np.isfinite(q[~la]).all()


def test_make_job_round5_settings():
    from fermentation_insights.coeff_predictor import experiments as E
    from fermentation_insights.coeff_predictor.evaluate import _configs
    ctx = dict(data_hash='d', split_hash='s', protocol=E.PROTOCOL)
    ctx5 = dict(ctx, lk_extra=dict(E.CHANNEL_LOSS), target_mask={k: list(v) for k, v in E.DEFECTIVE_TARGETS.items()},
                mk_extra=dict(power_gate=True))
    rest = (['1', '2'], {'test': ['3']}, 2, dict(width=8, depth=2), dict(w_ch=0.1), {}, None, 0.0, 0)
    j0 = E.make_job('cv', 't', 'pinn', *rest, ctx)
    j5 = E.make_job('cv', 't', 'pinn', *rest, ctx5)
    assert j0['key'] != j5['key'] and not j0['target_mask'] and j0['lk'] == dict(w_ch=0.1)
    assert j5['lk']['ch_space'] == 'surface' and j5['lk']['w_ch'] == E.CHANNEL_LOSS['w_ch']
    assert j5['target_mask'] == {'feedstock': list(E.DEFECTIVE_TARGETS['feedstock'])}
    assert j5['mk']['power_gate']
    # the plain MLP and the baselines get neither the channel loss nor the mask
    jp = E.make_job('cv', 't', 'plain', *rest, ctx5)
    jr = E.make_job('cv', 't', 'ridge', *rest, ctx5)
    assert 'ch_space' not in jp['lk'] and not jp['target_mask'] and not jr['target_mask']
    # the channel ablation still switches the channel loss off
    ja = E.make_job('cv', 't', 'pinn-no_channels', *rest, ctx5)
    mc, lc = _configs('pinn-no_channels', ja['mk'], ja['lk'])
    assert lc.w_ch == 0. and not mc.channels and ja['target_mask']


def test_hpo_candidates_with_a_fixed_channel_weight():
    from fermentation_insights.coeff_predictor import experiments as E
    c0 = E.hpo_candidates(16, 0)
    c1 = E.hpo_candidates(16, 0, fixed=E.hpo_fixed(dict(lk_extra=dict(E.CHANNEL_LOSS))))
    assert len(c1) == 16 and all(c['w_ch'] == E.CHANNEL_LOSS['w_ch'] for c in c1)
    strip = lambda c: {k: v for k, v in c.items() if k != 'w_ch'}
    # same draws: the other dimensions agree up to the first collision
    assert [strip(c) for c in c1[:8]] == [strip(c) for c in c0[:8]]
    assert len({tuple(sorted(c.items())) for c in c1}) == 16
    assert E.hpo_fixed({}) == {} and E.hpo_candidates(16, 0, fixed={}) == c0


def test_cli_round5_flags():
    from fermentation_insights.coeff_predictor import experiments as E
    a = E.parse_args([])
    assert not (a.no_power_gate or a.v2_channel_loss or a.no_target_mask)
    a = E.parse_args(['--v2-model'])
    assert a.no_power_gate and a.v2_channel_loss and a.no_target_mask and a.no_clip_exempt and not a.fs_b_only
    assert E.parse_args(['--fs-b-only']).fs_b_only and E.FS_B_ONLY is False
    a = E.parse_args(['--no-power-gate'])
    assert a.no_power_gate and not a.v2_channel_loss


def test_run_job_masks_defective_training_targets(tabs, monkeypatch):
    from fermentation_insights.coeff_predictor import experiments as E
    from fermentation_insights.coeff_predictor import evaluate as EV
    samples, points, _ = tabs
    seen = {}

    def stub(method, s, p, *a, **k):
        seen['s'] = s
        raise RuntimeError('stop')

    monkeypatch.setattr(EV, 'fit_method', stub)
    saved = dict(E._W)
    try:
        S = samples.assign(sample_id=samples['sample_id'].astype(str))
        P = points.assign(sample_id=points['sample_id'].astype(str))
        E._W.clear()
        E._W.update(samples=S.set_index('sample_id', drop=False), points=P,
                    pts={k: v for k, v in P.groupby('sample_id', sort=False)})
        ctx = dict(data_hash='d', split_hash='s', protocol=E.PROTOCOL, target_mask={'feedstock': ['LA']})
        ids = S['sample_id'].tolist()
        j = E.make_job('cv', 't', 'pinn', ids, {'test': ids[:2]}, 1, dict(width=8, depth=2), {}, {}, None, 0.,
                       0, ctx)
        r = E._run_job(j)
        assert 'stop' in r['error']
        la = seen['s']['fermentation_product'].eq('LA').to_numpy()
        assert la.any() and seen['s'].loc[la, 'int__fs__b'].isna().all()
        assert seen['s'].loc[~la, 'int__fs__b'].notna().all()
    finally:
        E._W.clear()
        E._W.update(saved)


def _fid_points(samples, points):
    """Points with simulated per-point intensities q__<key> = intensity
    surface / g (consistent with the fake samples' targets)."""
    P = points.copy()
    S = samples.set_index('sample_id')
    for key in ('fs', 'el', 'ng'):
        C = S.loc[P['sample_id'], [f'int__{key}__{k}' for k in 'abcd']].to_numpy(float)
        g = S.loc[P['sample_id'], 'g'].to_numpy(float)
        y, t = P['y'].to_numpy(float), P['t'].to_numpy(float)
        P[f'q__{key}'] = (C[:, 0] + C[:, 1] / y + C[:, 2] / t + C[:, 3] / (y * t)) / g
    return P


def _fid_inputs(samples):
    lay = make_layout(ModelConfig(), samples.columns)
    names = ['feedstock', 'electricity', 'natural_gas'] + lay['channels']
    _, g, ch, q = T._targets(samples, lay)
    prices = samples[['feedstock_price', 'electricity_price', 'natural_gas_price']].to_numpy(float)
    return lay, names, g, ch, q, prices


def test_channel_fidelity_metrics(tabs):
    from fermentation_insights.coeff_predictor.evaluate import channel_fidelity, fidelity_summary
    samples, points, _ = tabs
    samples = samples.assign(channels_ok=True)
    P = _fid_points(samples, points)
    lay, names, g, ch, q, prices = _fid_inputs(samples)
    C = np.concatenate([prices[:, :, None] * q, ch], 1)
    f = channel_fidelity(samples, P, C, names, q, lay['price_keys'], g)
    vals = f[[c for c in f.columns if c.startswith(('ch_', 'sens_', 'resp_'))]].to_numpy(float)
    assert np.nanmax(np.abs(vals)) < 1e-9 and np.nanmax(np.abs(f['bias_fs'])) < 1e-9
    # doubled feedstock intensity (and its channel): sensitivity +100 %, response = the feedstock share
    q2, C2 = q.copy(), C.copy()
    q2[:, 0] *= 2.
    C2[:, 0] *= 2.
    f2 = channel_fidelity(samples, P, C2, names, q2, lay['price_keys'], g)
    np.testing.assert_allclose(f2['sens_fs'], 1., rtol=1e-9)
    np.testing.assert_allclose(f2['bias_fs'], 1., rtol=1e-9)
    np.testing.assert_allclose(f2['resp_fs'], f2['share_fs'], rtol=1e-9)
    assert (f2['ch_feedstock'] > 0).all() and np.nanmax(np.abs(f2['ch_electricity'])) < 1e-9
    # a 10 % higher g lowers every sensitivity by 1/1.1
    f3 = channel_fidelity(samples, P, C, names, q, lay['price_keys'], g * 1.1)
    np.testing.assert_allclose(f3['bias_el'], 1 / 1.1 - 1, rtol=1e-6)
    s = fidelity_summary(f2)
    assert abs(float(s['sens_fs'].iloc[0]) - 1.) < 1e-9 and abs(float(s['bias_fs'].iloc[0]) - 1.) < 1e-9
    # channels_ok False -> NaN channel errors
    bad = samples.assign(channels_ok=False)
    assert channel_fidelity(bad, P, C, names, q, lay['price_keys'], g)['ch_feedstock'].isna().all()


def test_fidelity_tables_mark_defective_configurations(tabs):
    from fermentation_insights.coeff_predictor import experiments as E
    from fermentation_insights.coeff_predictor.evaluate import channel_fidelity
    samples, points, _ = tabs
    masked = E.DEFECTIVE_TARGETS['feedstock'][0]           # a biorefinery whose feedstock targets are masked
    samples = samples.assign(channels_ok=True, biorefinery=np.where(samples['fermentation_product'].eq('LA'),
                                                                     masked, 'HP'))
    P = _fid_points(samples, points)
    lay, names, g, ch, q, prices = _fid_inputs(samples)
    q2 = q.copy()
    q2[:, 0] *= 1.5
    C = np.concatenate([prices[:, :, None] * q2, ch], 1)
    f = channel_fidelity(samples, P, C, names, q2, lay['price_keys'], g).assign(scheme='config')
    summ, per = E.fidelity_tables(f)
    assert abs(float(summ['sens_fs'].iloc[0]) - 0.5) < 1e-9
    assert abs(float(summ['sens_fs_valid'].iloc[0]) - 0.5) < 1e-9
    assert set(per['config_id']) == set(samples['config_id']) and (per['n_samples'] == 4).all()
    la = f['biorefinery'].eq(masked)
    f2 = f.copy()
    f2.loc[la, 'sens_fs'] = 9.
    s2, _ = E.fidelity_tables(f2)
    assert float(s2['sens_fs'].iloc[0]) > 0.5 and abs(float(s2['sens_fs_valid'].iloc[0]) - 0.5) < 1e-9
