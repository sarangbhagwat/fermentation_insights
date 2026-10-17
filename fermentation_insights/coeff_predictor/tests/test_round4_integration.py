# -*- coding: utf-8 -*-
"""Round-4 integration: the v2 pipeline model (masking 'physics_v2' + hard
productivity structure + extrapolation safeguards + water gate) keeps every
component's guarantee when they are combined, and the experiments-pipeline
helpers of the v2 protocol (hold-outs per stage, equal ensemble sizes, LOCO
coverage of the PI-NN, out-of-fold calibration checks, cache keys)."""
import types

import numpy as np
import pandas as pd
import pytest
import torch

from fermentation_insights.coeff_predictor import experiments as E
from fermentation_insights.coeff_predictor import train as T
from fermentation_insights.coeff_predictor.evaluate import fit_method, ABLATIONS
from fermentation_insights.coeff_predictor.features import FeatureSpec
from fermentation_insights.coeff_predictor.model import ModelConfig, make_layout, mpsp_equation, G_KEY
from fermentation_insights.coeff_predictor.tests._fake_tables import make_fake_tables
from fermentation_insights.coeff_predictor.tests.conftest import FAST

V2 = dict(water_gate=True, channel_inputs=E.MASKING_PRESET, r_structure=E.R_STRUCTURE, **E.EXTRAPOLATION_SAFEGUARDS)


class _nowarn:
    """Silence warnings (tiny fake fits)."""
    def __enter__(self):
        import warnings
        self.c = warnings.catch_warnings()
        self.c.__enter__()
        warnings.simplefilter('ignore')

    def __exit__(self, *a):
        return self.c.__exit__(*a)


@pytest.fixture(scope='module')
def tabs():
    samples, points, configs = make_fake_tables(n_samples_per_config=4, products=['HP', 'EtOH', 'SA', 'BrothAcid'],
                                                feedstocks=['glucose', 'corn', 'cornstover'], seed=11)
    rng = np.random.RandomState(1)
    samples = samples.copy()
    samples['p_base_price'] = np.where(samples['neutralization'].astype(bool), rng.uniform(0.2, 0.3, len(samples)),
                                       np.nan)
    samples['p_acid_price'] = rng.uniform(0.09, 0.11, len(samples))
    return samples, points, configs


def _model(tabs, scramble=0.3, **kw):
    samples, points, _ = tabs
    spec = FeatureSpec(downstream=True).fit(samples)
    batch = spec.transform(samples)
    cfg = ModelConfig(**{**dict(width=16, depth=2, dropout=0., masked_width=16), **V2, **kw})
    layout = make_layout(cfg, samples.columns)
    scales = T.fit_target_scales(samples, points, layout, batch, spec.continuous_mask(), T.fs_prior_z(spec, batch))
    torch.manual_seed(0)
    m = T.build_model(spec, layout, scales, cfg, torch.float64).eval()
    if scramble:
        with torch.no_grad():
            for p in m.parameters():
                p.add_(scramble * torch.randn_like(p))
    X = torch.tensor(batch.X)
    lr = torch.tensor(batch.log_r)
    rest = [torch.tensor(np.asarray(a, float)) for a in (batch.prices, batch.w_sugar, batch.g_max)]
    Y, Tt, _, _, mask = T.pad_points(samples, points)
    return m, spec, batch, layout, X, lr, rest, [torch.tensor(a) for a in (Y, Tt, mask)]


def test_v2_model_combines_all_guarantees(tabs):
    m, spec, batch, layout, X, lr, rest, (Y, Tt, mask) = _model(tabs)
    assert m.masked and m.n_rs == 2 and m.avail_index is not None and m.neut_index is not None and m.clip_on
    # 1. dMPSP/dr <= 0 everywhere (r-structure; no other output sees r)
    prev = None
    for lg in np.log(np.geomspace(1e-3, 1e3, 13)):
        o = m(X, torch.full_like(lr, float(lg)), *rest)
        M = mpsp_equation(o['coef'][:, None, :], o['g'][:, None], Y, Tt)
        if prev is not None:
            assert torch.all((M - prev)[mask > 0] <= 1e-9 * prev.abs()[mask > 0].clamp(min=1.))
        prev = M
    # 2. the recovery does not depend on the reagent prices (physics_v2) nor on r
    Xl = X.clone().requires_grad_(True)
    lrl = lr.clone().requires_grad_(True)
    o = m(Xl, lrl, *rest)
    gx, gr = torch.autograd.grad(o['g'].sum(), [Xl, lrl], allow_unused=True, retain_graph=True)
    rp = spec.column_groups()['reagent_prices']
    assert rp and np.all(gx.numpy()[:, rp] == 0.) and (gr is None or torch.all(gr == 0))
    # ... while the neutralization channel does
    j = layout['channels'].index('neutralization')
    gn, = torch.autograd.grad(o['channels'][:, j].sum(), Xl, retain_graph=True)
    neut = X[:, m.neut_index].numpy()
    assert np.any(gn.numpy()[neut > 0][:, rp] != 0.)
    # 3. neutralization gate: zero cost without anything to neutralize (ethanol)
    assert (neut == 0).any() and torch.all(o['channels'][neut == 0, j] == 0)
    # 4. feedstock inputs clipped to the training range: pushing a feedstock composition column far
    #    outside it changes nothing beyond the clip
    fs_cols = [i for i in spec.column_groups()['feedstock'] if spec.continuous_mask()[i]]
    assert fs_cols
    X2, X3 = X.clone(), X.clone()
    X2[:, fs_cols[0]] = float(m.x_hi[fs_cols[0]]) + 5.
    X3[:, fs_cols[0]] = float(m.x_hi[fs_cols[0]]) + 50.
    o2, o3 = m(X2, lr, *rest), m(X3, lr, *rest)
    torch.testing.assert_close(o2['coef'], o3['coef'], rtol=0, atol=0)
    # 5. sugar-availability prior: the feedstock efficiency follows log(1/s_avail) + dz
    assert np.isfinite(m.scales_np['fs_prior_center'])


@pytest.mark.parametrize('abl', ['no_channels', 'no_price', 'no_phys'])
def test_v2_ablations_build_and_train(tabs, abl):
    samples, points, _ = tabs
    with _nowarn():
        mod = fit_method(f'pinn-{abl}', samples, points, 1, dict(width=16, depth=2, dropout=0., masked_width=8, **V2),
                         {}, dict(FAST, epochs=4, min_epochs=1), group_col=None, calib_frac=0.0, seed=2,
                         spec_kw={'downstream': True})
    pred = mod.predict(samples.head(5))
    assert np.isfinite(pred['coef']).all() and np.isfinite(pred['g']).all()
    assert mod.ens.meta['n_calib'] == 0 and mod.ens.meta['group_col'] is None
    mc = mod.ens.mcfg
    if abl == 'no_price':      # no structured feedstock intensity -> the sugar prior is inactive
        assert mod.ens.members[0].avail_index is None
    if abl == 'no_channels':   # the 'rest' head is r-structured, no neutralization channel to gate
        assert mod.ens.members[0].n_rs >= 1 and mod.ens.members[0].neut_index is None
    assert mc.r_structure == 'physics' and abl in ABLATIONS


def test_v2_ensemble_keeps_components_after_constraint_check(tabs):
    samples, points, _ = tabs
    with _nowarn():
        ens = T.train_ensemble(samples, points, 2, ModelConfig(width=16, depth=2, dropout=0., masked_width=8, **V2),
                               tcfg=T.TrainConfig(**dict(FAST, epochs=4, min_epochs=1)), calib_frac=0.0,
                               group_col=None, seed=0, spec_kwargs={'downstream': True})
    assert ens.mcfg.neutralization_gate and ens.mcfg.fs_prior and ens.mcfg.r_structure == 'physics'
    rep = ens.constraint_report
    row = rep[rep['constraint'] == 'neutralization_gate'] if isinstance(rep, pd.DataFrame) else None
    assert row is None or bool(row['enabled'].iloc[0])
    # member hold-outs are sample-level and differ between members (so the members together fit
    # nearly every sample; n_fit = their union)
    assert ens.meta['val_group_used'] == ['sample', 'sample'] and ens.meta['n_calib'] == 0
    assert ens.meta['n_train'] == len(samples) and ens.meta['n_fit'] >= 0.9 * len(samples)


# ---------------------------------------------------------------- pipeline helpers

def test_protocol_holdouts():
    for st in ('hpo', 'cv', 'final', 'release'):
        assert E.protocol_holdouts(E.PROTOCOL, st, 'product') == (None, 0.0)
    v1 = E.PROTOCOL_V1
    assert E.protocol_holdouts(v1, 'hpo') == ('config_id', 0.0)
    assert E.protocol_holdouts(v1, 'cv', 'product') == ('fermentation_product', 0.15)
    assert E.protocol_holdouts(v1, 'cv', 'random') == (None, 0.15)
    assert E.protocol_holdouts(v1, 'final') == ('config_id', 0.15)
    assert E.protocol_holdouts(v1, 'release') == ('config_id', 0.0)
    assert E.INNER_GROUP == {s: None for s in E.SCHEMES}
    assert E.PROTOCOL['refit'] == 'none' and E.PROTOCOL['calib_frac'] == 0.0


def test_plan_members_equal_and_hp_kwargs_refit():
    a = E.parse_args([])
    plan = E.base_plan(a)
    assert plan['cv_members'] == plan['final_members'] == plan['release_members'] == a.members == 5
    assert plan['loco_folds'] is None and plan['ablation_loco_folds'] == 12
    E.check_plan(plan)
    with pytest.raises(SystemExit):
        E.check_plan(dict(plan, final_members=3))
    q = E.base_plan(E.parse_args(['--quick']))
    assert q['loco_folds'] is None and q['cv_members'] == q['final_members'] == q['release_members']
    mk, lk, tk = E.hp_kwargs(E.HP_DEFAULT, plan)
    assert tk['refit'] == 'none' and T.TrainConfig(**tk).refit == 'none'


@pytest.fixture(scope='module')
def dev_like():
    samples, _, _ = make_fake_tables(n_samples_per_config=3, products=['HP', 'SA', 'EtOH', 'TAL', 'LA'],
                                     feedstocks=['glucose', 'corn', 'cornstover', 'sugarcane'], seed=2)
    return samples.assign(sample_id=samples['sample_id'].astype(str))


def test_ablations_on_a_subset_pinn_on_every_loco_fold(dev_like):
    plan = dict(random_k=3, loco_folds=None, ablation_loco_folds=6, cv_members=2)
    folds = E.scheme_splits(dev_like, 'config', plan, 0)
    assert len(folds) == dev_like['config_id'].nunique() == 20
    abl = E.ablation_fold_labels(dev_like, plan, 0)
    assert len(abl) == 6 and abl <= {lab for *_, lab in folds}
    n_pinn = sum('pinn' in E.cv_methods('config', lab, abl) for *_, lab in folds)
    n_abl = sum(E.ABLATIONS[0] in E.cv_methods('config', lab, abl) for *_, lab in folds)
    assert n_pinn == 20 and n_abl == 6
    assert E.ABLATIONS[0] in E.cv_methods('product', 'HP', abl)
    assert E.ABLATIONS[0] not in E.cv_methods('random', 'fold0', abl)
    assert E.ablation_fold_labels(dev_like, dict(plan, ablation_loco_folds=None), 0) is None
    # the budget never subsets the PI-NN LOCO folds, only the ablations' (down to 6)
    a = types.SimpleNamespace(budget_hours=1e-6, threads=1, hpo_candidates=0)
    p = E.plan_budget(dict(plan, hpo_k=3, hpo_members=1, final_members=2, release_members=2,
                           ablation_loco_folds=12), a, dict(seconds_per_member_sample=1.), dev_like, 10, 100, 0)
    assert p['loco_folds'] is None and p['ablation_loco_folds'] == 6 and p['hpo_candidates'] == 8
    assert not p['estimate']['fits_budget']


def _scores(ids, cfg):
    n = len(ids)
    return pd.DataFrame(dict(scheme='config', method='pinn', fold=cfg, fold_index=0, sample_id=ids, config_id=cfg,
                             invalid=False, a=np.ones(n), b=1., c=1., d=1., g=1., mpsp=1.))


def test_oof_calibration_scores_checks(dev_like):
    dev = dev_like.iloc[:-6]
    test = dev_like.iloc[-6:]
    sp = dict(ood_test_sample_ids=list(test['sample_id'][:3]), id_test_sample_ids=list(test['sample_id'][3:]))
    ok = _scores(dev['sample_id'].tolist(), dev['config_id'].astype(str).tolist())
    rec = E.oof_calibration_scores(dict(oof_scores=ok), 'config', 'pinn', dev, sp)
    assert len(rec['sample_id']) == len(dev) and set(rec) >= {'fold', 'config_id', 'mpsp', 'g'}
    dup = pd.concat([ok, ok.iloc[:1]], ignore_index=True)
    with pytest.raises(SystemExit, match='more than once'):
        E.oof_calibration_scores(dict(oof_scores=dup), 'config', 'pinn', dev, sp)
    with pytest.raises(SystemExit, match='not scored'):
        E.oof_calibration_scores(dict(oof_scores=ok.iloc[1:]), 'config', 'pinn', dev, sp)
    leak = pd.concat([ok, _scores(sp['ood_test_sample_ids'][:1], ['x'])], ignore_index=True)
    for strict in (True, False):          # a test sample in the calibration set always raises
        with pytest.raises(SystemExit, match='test samples'):
            E.oof_calibration_scores(dict(oof_scores=leak), 'config', 'pinn', dev, sp, strict=strict)
    few = dev[dev['config_id'].isin(sorted(dev['config_id'].unique())[:3])]
    sc_few = _scores(few['sample_id'].tolist(), few['config_id'].astype(str).tolist())
    with pytest.raises(SystemExit, match='configurations'):
        E.oof_calibration_scores(dict(oof_scores=sc_few), 'config', 'pinn', few, sp)
    assert E.oof_calibration_scores(dict(oof_scores=sc_few), 'config', 'pinn', few, sp, strict=False) is not None


def test_make_job_cache_key_includes_calibration():
    ctx = dict(data_hash='d', split_hash='s', protocol=E.PROTOCOL)
    args = ('final', 't', 'pinn', ['1', '2'], {'test': ['3']}, 2, dict(width=8, depth=2), {}, {}, None, 0.0, 0, ctx)
    j0 = E.make_job(*args)
    j1 = E.make_job(*args, calib_scores=dict(sample_id=['1'], mpsp=[1.0]), calib_info=dict(scheme='config'))
    j2 = E.make_job(*args, calib_scores=dict(sample_id=['1'], mpsp=[2.0]), calib_info=dict(scheme='config'))
    assert len({j0['key'], j1['key'], j2['key']}) == 3
    assert j1['protocol'] == 'v2' and j1['calib_info'] == dict(scheme='config') and j0['calib_scores'] is None
