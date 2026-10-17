# -*- coding: utf-8 -*-
"""Hard productivity structure (ModelConfig.r_structure): the coefficients of
the r-dependent outputs are theta0(x) + theta1(x) (r_ref / r)^beta with
theta1 >= 0 and theta0, theta1 independent of r, so that dMPSP/dr <= 0 holds
by construction; works with the masked model, the price-explicit heads, the
water gate and the shared encoder; checkpoints round-trip and old ones load."""
import warnings

import numpy as np
import pytest
import torch

from fermentation_insights.coeff_predictor.features import FeatureSpec
from fermentation_insights.coeff_predictor.model import (
    ModelConfig, make_layout, resolve_r_structure, R_STRUCTURE_PRESETS, INTENSITY_CHANNEL, G_KEY, mpsp_equation)
from fermentation_insights.coeff_predictor import train as T
from fermentation_insights.coeff_predictor.tests._fake_tables import make_fake_tables


@pytest.fixture(scope='module')
def ds_tables():
    return make_fake_tables(n_samples_per_config=4, products=['HP', 'EtOH', 'SA', 'BrothAcid'],
                            feedstocks=['glucose', 'cornstover'], seed=5)


def _model(tables, seed=0, downstream=True, scramble=0.5, **kw):
    """Untrained model whose weights are scrambled (so that theta1, the gate and
    every head are far from their initial values) + inputs."""
    samples, points, _ = tables
    spec = FeatureSpec(downstream=downstream).fit(samples)
    batch = spec.transform(samples)
    base = dict(width=16, depth=2, dropout=0., masked_width=16, r_structure='physics')
    if downstream:
        base.update(water_gate=True, channel_inputs='physics')
    cfg = ModelConfig(**{**base, **kw})
    layout = make_layout(cfg, samples.columns)
    scales = T.fit_target_scales(samples, points, layout, batch, spec.continuous_mask())
    torch.manual_seed(seed)
    m = T.build_model(spec, layout, scales, cfg, torch.float64).eval()
    if scramble:
        with torch.no_grad():
            for p in m.parameters():
                p.add_(scramble * torch.randn_like(p))
    X = torch.tensor(batch.X)
    lr = torch.tensor(batch.log_r)
    rest = [torch.tensor(np.asarray(a, float)) for a in (batch.prices, batch.w_sugar, batch.g_max)]
    Y, Tt, M, _, mask = T.pad_points(samples, points)
    return m, layout, X, lr, rest, [torch.tensor(a) for a in (Y, Tt, mask)]


def _slot_values(m, layout, o):
    """{output name: (n, 4) coefficient tensor} of every output (intensities
    by their channel name, before the price)."""
    d = {INTENSITY_CHANNEL[k]: o['intensities'][:, j] for j, k in enumerate(layout['price_keys'])}
    d.update({c: o['channels'][:, i] for i, c in enumerate(layout['channels'])})
    return d


def test_config_resolution_and_compatibility():
    assert resolve_r_structure(None) == () and resolve_r_structure([]) == ()
    assert resolve_r_structure('physics') == R_STRUCTURE_PRESETS['physics']
    assert resolve_r_structure(['capex_fermentation']) == ('capex_fermentation',)
    with pytest.raises(ValueError, match='preset'):
        ModelConfig(r_structure='capex_fermentation')
    with pytest.raises(ValueError, match='r_beta'):
        ModelConfig(r_structure='physics', r_beta=0.)
    with pytest.raises(ValueError, match='r_beta_bounds'):
        ModelConfig(r_structure='physics', r_beta_bounds=(0., 1.))
    with pytest.raises(ValueError, match='outside'):
        ModelConfig(r_structure='physics', r_beta_learn=True, r_beta=2., r_beta_bounds=(0.5, 1.5))
    assert ModelConfig(kind='plain', r_structure='physics').r_structured == ()
    assert ModelConfig(r_structure='physics').r_structured == ('capex_fermentation', 'electricity')
    assert ModelConfig().r_structured == ()
    # JSON round trip (tuples become lists) and checkpoints saved before the option
    c = ModelConfig(r_structure=('capex_fermentation',), r_beta_learn=True, r_beta_bounds=[0.6, 1.4])
    d = c.to_dict()
    assert d['r_structure'] == ['capex_fermentation']
    c2 = ModelConfig.from_dict(dict(d, r_beta_bounds=list(d['r_beta_bounds'])))
    assert c2 == c and c2.r_beta_bounds == (0.6, 1.4)
    old = {k: v for k, v in ModelConfig().to_dict().items() if not k.startswith('r_')}
    assert ModelConfig.from_dict(old).r_structured == ()


@pytest.mark.parametrize('kw', [dict(), dict(r_beta_learn=True, r_beta=0.8), dict(downstream=False),
                                dict(price_explicit=False), dict(residual_width=8)],
                         ids=['masked', 'masked-learned-beta', 'shared-encoder', 'no-price', 'hybrid'])
def test_theta_independent_of_r_and_monotone(ds_tables, kw):
    m, layout, X, lr, rest, (Y, Tt, mask) = _model(ds_tables, **kw)
    assert m.n_rs == 2, m.rs_names
    assert set(m.rs_names) == {'capex_fermentation', 'electricity'}
    lr_leaf = lr.clone().requires_grad_(True)
    o = m(X, lr_leaf, *rest)
    th0, th1 = o['r_theta0'], o['r_theta1']
    # theta1 >= 0 and exactly r-independent theta0, theta1 (no graph path from log r)
    assert (th1 >= 0).all()
    for t in (th0, th1):
        gr = torch.autograd.grad(t.sum(), lr_leaf, retain_graph=True, allow_unused=True)[0]
        assert gr is None or torch.all(gr == 0)
    o2 = m(X, lr + 3.7, *rest)
    assert torch.equal(o2['r_theta0'], th0) and torch.equal(o2['r_theta1'], th1)
    # structured outputs = theta0 + theta1 (r_ref / r)^beta; the others do not depend on r
    vals, vals2 = _slot_values(m, layout, o), _slot_values(m, layout, o2)
    beta = m.r_beta_values()
    rho = torch.exp(beta[None, :] * (np.log(m.r_ref) - lr[:, None]))
    for i, name in enumerate(m.rs_names):
        th = th0[:, i] + th1[:, i] * rho[:, i, None]
        torch.testing.assert_close(vals[name], th, rtol=1e-10, atol=1e-12)
    for name, v in vals.items():
        if name not in m.rs_names:
            torch.testing.assert_close(vals2[name], v, rtol=0, atol=0)
    torch.testing.assert_close(o2['g'], o['g'], rtol=0, atol=0)
    # dMPSP/dr <= 0 everywhere (far outside the training range as well) ...
    r_grid = np.log(np.geomspace(1e-3, 1e3, 25))
    M_prev = None
    for lg in r_grid:
        oo = m(X, torch.full_like(lr, float(lg)), *rest)
        Mr = mpsp_equation(oo['coef'][:, None, :], oo['g'][:, None], Y, Tt)
        if M_prev is not None:
            d = (Mr - M_prev)[mask > 0]
            assert torch.all(d <= 1e-9 * M_prev.abs()[mask > 0].clamp(min=1.)), float(d.max())
        M_prev = Mr
    # ... and the autograd residual of the loss is exactly zero
    lr_leaf = lr.clone().requires_grad_(True)
    o = m(X, lr_leaf, *rest)
    dM = T.mpsp_dlogr(o, lr_leaf, Y, Tt)
    assert torch.all(dM[mask > 0] <= 1e-12)
    if kw.get('r_beta_learn'):
        lo, hi = m.cfg.r_beta_bounds
        assert torch.all((beta >= lo) & (beta <= hi))
        g = torch.autograd.grad(o['coef'].sum(), m.r_beta_raw, allow_unused=True)[0]
        assert g is not None and torch.any(g != 0)


def test_productivity_reaches_only_structured_outputs(ds_tables):
    m, layout, X, lr, rest, _ = _model(ds_tables)
    lr_leaf = lr.clone().requires_grad_(True)
    o = m(X, lr_leaf, *rest)
    for name, v in list(_slot_values(m, layout, o).items()) + [(G_KEY, o['g'])]:
        gr = torch.autograd.grad(v.sum(), lr_leaf, retain_graph=True, allow_unused=True)[0]
        dr = 0. if gr is None else float(gr.abs().sum())
        assert (dr > 0) == (name in m.rs_names), name
    # the masked sub-networks of the structured outputs have no productivity inputs
    for name in m.rs_names:
        assert not m.output_inputs[name]['productivity'] and 'productivity' not in m.output_inputs[name]['groups']
    assert not any(m.sub_r)
    # the shared encoder variant has no productivity input either
    m2, *_ = _model(ds_tables, downstream=False)
    assert m2.encoder[0].in_features == len(m2.enc_index)


def test_initial_model_and_layout_variants(ds_tables):
    # unscrambled: theta1 = the training slope (>= r_slope_floor x std)
    m, layout, X, lr, rest, _ = _model(ds_tables, scramble=0.)
    o = m(X, torch.full_like(lr, np.log(m.r_ref)), *rest)
    torch.testing.assert_close(o['r_theta1'], m.rs_std.expand_as(o['r_theta1']), rtol=0.2, atol=0.)
    sd = np.array([m.int_std[s].numpy() if s < m.n_price else m.ch_std[s - m.n_price].numpy()
                   for s in m.rs_slots])
    assert np.all(m.rs_std.numpy() >= 0.01 * sd * (1 - 1e-6))
    # sum identity holds with the structure
    s = o['channels'].sum(1) + o['price_channels'].sum(1)
    torch.testing.assert_close(s, o['coef'], rtol=0, atol=1e-12)
    # without resolved channels the 'rest' head carries the fermentation capital structure
    m3, *_ = _model(ds_tables, channels=False)
    assert 'rest' in m3.rs_names and 'electricity' in m3.rs_names
    # feedstock intensity cannot be both structured and r-structured
    with pytest.raises(ValueError, match='structured_fs'):
        _model(ds_tables, r_structure=['feedstock'])
    # an r_structure output the model does not have: structure off (with a warning)
    with pytest.warns(UserWarning, match='r_structure'):
        m4, *_ = _model(ds_tables, r_structure=['no_such_channel'])
    assert m4.n_rs == 0 and not m4.r_free


def test_r_slopes_recover_within_configuration_slope():
    rng = np.random.RandomState(0)
    n = 60
    grp = np.repeat(np.arange(6), 10)
    log_r = rng.uniform(-2, 1, n)
    rho = np.exp(log_r.mean() - log_r)
    # per-configuration offsets (cross-configuration variation) + common slope 0.7 on rho (coef a),
    # no r-dependence (coef c), slope -0.2 (coef d)
    v = np.zeros((n, 1, 4))
    v[:, 0, 0] = 5. * grp + 0.7 * rho
    v[:, 0, 1] = rng.normal(size=n)
    v[:, 0, 2] = 3. * grp
    v[:, 0, 3] = -0.2 * rho + 1. * grp
    s = T.r_slopes(v, log_r, grp)
    np.testing.assert_allclose(s[0, [0, 2, 3]], [0.7, 0., -0.2], atol=1e-10)
    assert abs(s[0, 1]) < 0.5


def test_checkpoint_roundtrip(ds_tables, tmp_path):
    samples, points, _ = ds_tables
    cfg = ModelConfig(width=16, depth=2, dropout=0., water_gate=True, masked_width=16, channel_inputs='physics',
                      r_structure='physics', r_beta_learn=True, r_beta=0.9)
    ens = T.train_ensemble(samples, points, 1, cfg, T.LossConfig(colloc_ratio=0.5),
                           T.TrainConfig(epochs=4, batch_size=16, min_epochs=1, patience=3),
                           calib_frac=0., spec_kwargs=dict(downstream=True))
    m = ens.members[0]
    assert 'r_beta_raw' in m.state_dict() and not any(k.startswith('rs_') for k in m.state_dict())
    assert float(m.r_beta_values()[0]) != pytest.approx(0.9, abs=1e-9)   # beta was trained
    _, o1 = ens.predict_samples(samples)
    assert np.isfinite(o1['coef']).all()
    re = T.Ensemble.load(ens.save(str(tmp_path / 'ck')))
    assert re.mcfg == ens.mcfg and re.members[0].rs_names == m.rs_names
    _, o2 = re.predict_samples(samples)
    np.testing.assert_allclose(o2['coef'], o1['coef'], rtol=1e-6, atol=1e-8)
    from fermentation_insights.coeff_predictor.predict import CoefficientPredictor
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        o = CoefficientPredictor(re).predict(samples.iloc[:3])
    assert len(o) == 3 and np.isfinite(o[['a', 'b', 'c', 'd', 'g']].to_numpy(float)).all()


def test_old_checkpoint_without_option_loads(ds_tables, tmp_path):
    """A checkpoint whose meta.json predates r_structure loads unchanged."""
    import json
    import os
    samples, points, _ = ds_tables
    cfg = ModelConfig(width=16, depth=2, dropout=0., water_gate=True, masked_width=16, channel_inputs='physics')
    ens = T.train_ensemble(samples, points, 1, cfg, T.LossConfig(),
                           T.TrainConfig(epochs=2, batch_size=16, min_epochs=1, patience=2),
                           calib_frac=0., spec_kwargs=dict(downstream=True))
    path = ens.save(str(tmp_path / 'old'))
    with open(os.path.join(path, 'meta.json')) as fh:
        meta = json.load(fh)
    meta['model_config'] = {k: v for k, v in meta['model_config'].items() if not k.startswith('r_')}
    with open(os.path.join(path, 'meta.json'), 'w') as fh:
        json.dump(meta, fh)
    re = T.Ensemble.load(path)
    assert re.mcfg.r_structured == () and re.members[0].n_rs == 0
    _, o1 = ens.predict_samples(samples)
    _, o2 = re.predict_samples(samples)
    np.testing.assert_allclose(o2['coef'], o1['coef'], rtol=1e-6, atol=1e-8)
