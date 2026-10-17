# -*- coding: utf-8 -*-
"""Training, splitting, constraint checks, checkpoints and prediction."""
import warnings

import numpy as np
import pandas as pd
import pytest
import torch

from fermentation_insights.coeff_predictor import train as T
from fermentation_insights.coeff_predictor.model import ModelConfig, mpsp_equation_np, make_layout
from fermentation_insights.coeff_predictor.predict import CoefficientPredictor
from fermentation_insights.coeff_predictor.features import (MissingSpecError, ApplicabilityWarning,
                                                            PRICE_COLUMNS, REGIME_FIELDS)
from fermentation_insights.coeff_predictor.tests.conftest import FAST


@pytest.fixture(scope='module')
def ensemble(tables):
    samples, points, _ = tables
    return T.train_ensemble(samples, points, n_members=2, mcfg=ModelConfig(width=24, depth=2, dropout=0.),
                            tcfg=T.TrainConfig(**{**FAST, 'epochs': 40, 'patience': 15}),
                            group_col='config_id', seed=0)


def test_training_reduces_loss(ensemble):
    for h in ensemble.histories:
        assert h[-1]['train_loss'] < 0.6 * h[0]['train_loss']
        vals = [r['val_data'] for r in h]
        assert min(vals) < vals[0] and int(np.argmin(vals)) > 0  # best state was recorded
    assert ensemble.constraint_report is not None
    assert set(ensemble.conformal) >= {'a', 'b', 'c', 'd', 'g', 'mpsp'}


def test_constraint_report_on_fake_data(tables):
    samples, points, _ = tables
    rep = T.check_constraints(samples, points).set_index('constraint')
    for c in ('mono_y', 'mono_t', 'positive', 'nonneg_fs', 'nonneg_ng', 'g_le_gmax', 'fs_eta_le_1',
              'channel_sum', 'price_linear', 'mono_r'):
        assert rep.loc[c, 'frac_satisfied'] == 1.0 and rep.loc[c, 'enabled'], c
    assert rep.loc['mono_r', 'n_cases'] > 0
    # violate g <= g_max in 10 % of samples -> hard bound disabled
    bad = samples.copy()
    idx = bad.index[::10]
    bad.loc[idx, 'g'] = bad.loc[idx, 'final_product_stoich_mass_ratio'] * 1.05
    rep2 = T.check_constraints(bad, points)
    m, l = T.apply_constraint_report(rep2, ModelConfig(), T.LossConfig())
    assert not m.g_bound and m.fs_eta_bound
    # ... unless the excess is rounding-level: bound relaxed by a margin
    bad.loc[idx, 'g'] = bad.loc[idx, 'final_product_stoich_mass_ratio'] * 1.002
    m, l = T.apply_constraint_report(T.check_constraints(bad, points), ModelConfig(), T.LossConfig())
    assert m.g_bound and 0.002 < m.g_margin < 0.01
    # violate monotonicity in y -> residual disabled
    bad = samples.copy()
    bad.loc[bad.index[::5], 'b'] = -10.
    bad.loc[bad.index[::5], 'd'] = 0.
    m, l = T.apply_constraint_report(T.check_constraints(bad, points), ModelConfig(), T.LossConfig())
    assert not l.residuals['mono_y'] and 'mono_y' not in l.active() and 'mono_t' in l.active()


def test_checkpoint_round_trip_and_predict(ensemble, tables, tmp_path):
    samples, points, _ = tables
    path = ensemble.save(str(tmp_path / 'ckpt'))
    P = CoefficientPredictor.load(path)
    rec = samples.iloc[3].to_dict()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p1 = CoefficientPredictor(ensemble).predict(rec)
        p2 = P.predict(rec)
    for k in ('a', 'b', 'c', 'd', 'g', 'a_lo', 'g_hi', 'ch__capex_fermentation__a', 'int__fs__b'):
        assert np.isclose(p1[k].iloc[0], p2[k].iloc[0], rtol=1e-6), k
    ch = [c for c in p2 if c.startswith('ch__') and c.endswith('__a')]
    assert np.isclose(sum(p2[c].iloc[0] for c in ch), p2['a'].iloc[0], rtol=1e-5)
    assert p2['g'].iloc[0] <= rec['final_product_stoich_mass_ratio']
    assert p2['a_lo'].iloc[0] <= p2['a'].iloc[0] <= p2['a_hi'].iloc[0]
    # mpsp() equals the members' equation evaluated at the inputs
    m = P.mpsp(rec, y=[0.3, 0.5], t=[40., 100.], warn=False)
    batch, outs = P.ensemble.predict_samples(pd.DataFrame([rec]))
    ref = mpsp_equation_np(outs['coef'][:, 0, :], outs['g'][:, 0], 0.3, 40.).mean()
    assert np.isclose(m['MPSP'].iloc[0], ref, rtol=1e-6)
    assert (m['MPSP_lo'] <= m['MPSP']).all() and (m['MPSP'] <= m['MPSP_hi']).all()


def test_predict_by_name_missing_fields_and_domain_warning(ensemble):
    P = CoefficientPredictor(ensemble)
    spec = dict(feedstock='corn', fermentation_product='LA', final_product='LA', r=0.8,
                feedstock_price=0.2, neutralization=True, neutralizing_agent='Ca(OH)2',
                aerobic=False, fermentation_mode='batch', separation_principle='esterification')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = P.predict(spec)
    assert np.isfinite(out[['a', 'b', 'c', 'd', 'g']].to_numpy()).all()
    assert 'theoretical_yield' in out['filled_from_registry'].iloc[0]
    with pytest.raises(MissingSpecError, match='feedstock_price'):
        P.predict({k: v for k, v in spec.items() if k != 'feedstock_price'})
    with pytest.raises(MissingSpecError, match='r'):
        P.predict(dict(feedstock='unknown_stuff', r=None))
    far = dict(spec, r=500., feedstock_price=50.)
    with pytest.warns(ApplicabilityWarning):
        o = P.predict(far)
    assert not o['in_domain'].iloc[0]


def _ids(samples, idx):
    return set(samples['sample_id'].iloc[idx])


def test_splitters_do_not_leak(tables):
    samples, points, _ = tables
    n = len(samples)
    folds = T.make_splits(samples, 'random', 5, seed=1)
    allte = np.concatenate([te for _, te, _ in folds])
    assert sorted(allte) == list(range(n))
    for tr, te, _ in folds:
        assert not (_ids(samples, tr) & _ids(samples, te))
        ptr = points[points.sample_id.isin(_ids(samples, tr))]
        pte = points[points.sample_id.isin(_ids(samples, te))]
        assert len(ptr) + len(pte) == len(points)
        assert not set(map(tuple, ptr[['sample_id', 'y', 't']].to_numpy())) & \
            set(map(tuple, pte[['sample_id', 'y', 't']].to_numpy()))
    for scheme, col in (('config', 'config_id'), ('feedstock', 'feedstock'),
                        ('product', 'fermentation_product')):
        splits = T.make_splits(samples, scheme)
        assert len(splits) == samples[col].nunique()
        for tr, te, label in splits:
            held = set(samples[col].iloc[te])
            assert held == {label}
            assert not held & set(samples[col].iloc[tr])
            assert len(tr) + len(te) == n
    (tr, te, _), = T.make_splits(samples, 'productivity')
    rel = np.log(samples['r'] / samples['baseline_productivity']).to_numpy()
    assert len(te) > 0 and np.all((rel[te] < rel[tr].min()) | (rel[te] > rel[tr].max()))
    # inner grouped split keeps configurations together
    fit, val = T.inner_split(samples, np.arange(n), 0.2, 0, 'config_id')
    assert not set(samples.config_id.iloc[fit]) & set(samples.config_id.iloc[val])


def test_conformal_quantile():
    s = np.arange(1, 101, dtype=float)
    assert T._conf_q(s, 0.1) == 91.
    # too few scores for the level: the valid quantile is +inf (not the max)
    with pytest.warns(UserWarning, match='too few'):
        assert T._conf_q(np.array([1., 2.]), 0.1) == np.inf
    with pytest.warns(UserWarning):
        assert T._conf_q(np.arange(1., 6.), 0.1) == np.inf      # n=5: k=6 > 5
    assert T._conf_q(np.arange(1., 10.), 0.1) == 9.             # n=9: k=9
    # +inf scores are kept (not covered), only NaN is dropped
    q, n, lev = T._conf_info(np.r_[np.arange(1., 19.), np.inf, np.inf, np.nan], 0.1)
    assert n == 20 and q == np.inf          # rank 19 of 20 hits an inf score
    q, n, lev = T._conf_info(np.r_[np.arange(1., 20.), np.inf], 0.1)   # n=20, k=19
    assert n == 20 and q == 19. and np.isclose(lev, 19 / 21)
    assert np.isnan(T._conf_q(np.array([]), 0.1))


def test_plain_and_ablation_training_runs(tables):
    samples, points, _ = tables
    tc = T.TrainConfig(**{**FAST, 'epochs': 8})
    for mcfg, lcfg in ((ModelConfig(kind='plain', width=16, depth=1), T.LossConfig(physics=False, collocation=False)),
                       (ModelConfig(width=16, depth=1, price_explicit=False, channels=False), T.LossConfig())):
        ens = T.train_ensemble(samples, points, 1, mcfg, lcfg, tc, calib_frac=0.1, seed=1)
        batch, outs = ens.predict_samples(samples.iloc[:4])
        assert outs['coef'].shape == (1, 4, 4) and np.isfinite(outs['coef']).all()


def test_inner_split_grouping_and_fallbacks(tables):
    samples, _, _ = tables
    two = samples[samples['feedstock'].isin(['glucose', 'corn'])].reset_index(drop=True)
    fit, val, used = T.inner_split(two, np.arange(len(two)), 0.15, 0, 'feedstock', return_group=True)
    assert used == 'feedstock' and len(val) and not set(two.feedstock.iloc[fit]) & set(two.feedstock.iloc[val])
    one = samples[samples['feedstock'] == 'glucose'].reset_index(drop=True)
    with pytest.warns(UserWarning, match="grouped by 'config_id'"):
        fit, val, used = T.inner_split(one, np.arange(len(one)), 0.15, 0, 'feedstock', return_group=True)
    assert used == 'config_id' and not set(one.config_id.iloc[fit]) & set(one.config_id.iloc[val])
    cfg1 = one[one.config_id == one.config_id.iloc[0]].reset_index(drop=True)
    with pytest.warns(UserWarning, match='random sample-level'):
        _, _, used = T.inner_split(cfg1, np.arange(len(cfg1)), 0.2, 0, 'feedstock', return_group=True)
    assert used == 'sample'
    assert T.SCHEME_GROUP['feedstock'] == 'feedstock' and T.SCHEME_GROUP['product'] == 'fermentation_product'


def test_collocation_prices_nonnegative_and_in_box(fitted):
    samples, points, spec, batch = fitted['samples'], fitted['points'], fitted['spec'], fitted['batch']
    tb = T.make_tensors(spec, samples, points, make_layout(ModelConfig(), samples.columns), batch)
    tb['prices'][::2, 2] = 0.   # rows without a natural-gas price
    info = T._spec_info(spec, tb, torch.float32)
    gen = torch.Generator().manual_seed(0)
    for f0, f1 in ((1., 1.), (0.5, 1.5)):
        cb = T._collocation(tb, torch.arange(len(samples)), None, info,
                            T.LossConfig(colloc_price_range=(f0, f1)), gen)
        P = cb['prices'].numpy()
        lo, hi = info['price_range'][2]
        pos = P[:, 2][P[:, 2] > 0]
        assert (P >= 0).all() and (P[:, 2] == 0).any() and lo > 0
        assert (pos >= lo * f0 * (1 - 1e-6)).all() and (pos <= hi * f1 * (1 + 1e-6)).all()


def test_conformal_meta_and_price_explicit_ad(ensemble, tables):
    samples, _, _ = tables
    conf, meta = ensemble.conformal, ensemble.meta
    assert meta['calib_group_used'] == 'config_id' and meta['n_calib_groups'] >= 2
    assert 0 < meta['n_fit'] <= meta['n_train'] and set(meta['val_group_used']) == {'config_id'}
    for k in ('a', 'g', 'mpsp'):
        assert conf[f'{k}_n'] > 0 and 0 < conf[f'{k}_level'] <= 1
    assert conf['mpsp_n'] == conf['n']          # one MPSP score per calibration sample
    assert conf['calibration_grouping'] == 'config_id' and not conf['coverage_guaranteed']
    st = ensemble.spec.price_stats
    assert all(v['std'] >= 0.05 * abs(v['mean']) - 1e-15 for v in st.values())
    P = CoefficientPredictor(ensemble)
    rec = dict(samples.iloc[3].to_dict(), **ensemble.spec.price_fill)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p0 = P.predict(rec)
        p1 = P.predict(dict(rec, natural_gas_price=10 * rec['natural_gas_price']))
    assert np.isclose(p1['ad_ratio'].iloc[0], p0['ad_ratio'].iloc[0])   # prices not in the AD distance
    assert p1['price_outside_range'].iloc[0] == 'natural_gas_price' and p0['price_outside_range'].iloc[0] == ''


def test_invalid_inputs_and_per_row_reporting(ensemble, tables):
    samples, _, _ = tables
    P = CoefficientPredictor(ensemble)
    rec = samples.iloc[3].to_dict()
    for bad in (dict(r=0.), dict(r=-1.), dict(fs_w_fermentable_sugar_eq=0.),
                dict(final_product_stoich_mass_ratio=0.)):
        with pytest.raises(MissingSpecError, match='invalid'):
            P.predict(dict(rec, **bad))
    b = ensemble.spec.transform(pd.DataFrame([rec, dict(rec, r=0.)]), strict=False)
    assert list(b.invalid) == [False, True]
    part = {k: v for k, v in rec.items() if k not in ('electricity_price', 'p_electricity_price')}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = P.predict([rec, part])
    assert 'electricity_price' not in out['imputed'].iloc[0] and 'electricity_price' in out['imputed'].iloc[1]
    with pytest.warns(UserWarning, match='electricity_price'):
        P.mpsp(part, y=0.5 * rec['inflection_yield'], t=50.)
    with pytest.raises(ValueError):
        P.mpsp(rec, y=[0., 0.3], t=50., warn=False)
    ty = rec['theoretical_yield']
    with pytest.warns(ApplicabilityWarning, match='theoretical'):
        m = P.mpsp(rec, y=[0.5 * ty, 1.5 * ty], t=50.)
    assert np.isfinite(m['MPSP'].iloc[0]) and np.isnan(m['MPSP'].iloc[1])
    # every feedstock-price alias is checked (a NaN column must not shadow another)
    base = {k: v for k, v in rec.items() if k not in ('feedstock_price', 'p_feedstock_price', 'fs_price')}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        assert len(P.predict([dict(base, p_feedstock_price=0.2), dict(base, fs_price=0.2)])) == 2
    spec = dict(feedstock='corn stover', fermentation_product='LA', final_product='LA', r=0.8,
                feedstock_price=0.08)
    with pytest.raises(MissingSpecError, match="did you mean 'cornstover'"):
        P.predict(spec)
    assert 'r' in P.field_units and 'wet' in P.field_units['feedstock_price']


def test_predict_by_name_fills_regime_and_contextual(ensemble, tables):
    samples, _, _ = tables
    P = CoefficientPredictor(ensemble)
    by_name = dict(feedstock='corn', fermentation_product='HP', final_product='AcrylicAcid', r=0.5,
                   feedstock_price=0.2, electricity_price=0.07, natural_gas_price=0.25)
    hp = samples[samples.fermentation_product == 'HP'].iloc[0]
    explicit = dict(by_name, **{c: hp[c] for c in REGIME_FIELDS if c in hp})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        o1, o2 = P.predict(by_name), P.predict(explicit)
    assert {'aerobic', 'neutralizing_agent', 'fermentation_mode'} <= set(o1['filled_from_registry'].iloc[0].split(', '))
    np.testing.assert_allclose(o1[list('abcdg')].to_numpy(), o2[list('abcdg')].to_numpy(), rtol=1e-6)
    assert 'p_IRR' in o1['contextual_defaults'].iloc[0] and 'feedstock_capacity' in o1['contextual_defaults'].iloc[0]
    # applicable parameter omitted -> product baseline (indicator 0); not applicable -> indicator 1
    df, det = P._frame([by_name, dict(by_name, fermentation_product='EtOH', final_product='EtOH')], details=True)
    b = P.spec.transform(df, strict=False)
    j = P.spec.columns.index('missing__p_boiler_efficiency')
    assert b.X[0, j] == 0. and b.X[1, j] == 1.
    # unknown product: omitted categoricals -> unknown level, reported
    new = dict(feedstock='corn', fermentation_product='Mystery', final_product='Mystery', r=0.5,
               feedstock_price=0.2, theoretical_yield=0.6, inflection_yield=0.5,
               final_product_stoich_mass_ratio=1.0)
    with pytest.warns(UserWarning, match='not in the registry'):
        o = P.predict(new)
    assert {'fermentation_mode', 'separation_principle'} <= set(o['imputed'].iloc[0].split(', '))
    assert 'p_IRR' in o['contextual_defaults'].iloc[0]


def test_no_channels_breakdown_sums_to_totals(tables):
    samples, points, _ = tables
    ens = T.train_ensemble(samples, points, 1, ModelConfig(width=16, depth=1, channels=False),
                           T.LossConfig(w_ch=0.), T.TrainConfig(**{**FAST, 'epochs': 5}), seed=2)
    assert ens.channel_names()[-1] == 'rest'
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p = CoefficientPredictor(ens).predict(samples.iloc[2].to_dict())
    for k in 'abcd':
        ch = [c for c in p if c.startswith('ch__') and c.endswith(f'__{k}')]
        assert len(ch) == 4 and np.isclose(sum(p[c].iloc[0] for c in ch), p[k].iloc[0], rtol=1e-6)


def test_user_recovery_overrides_g_keeps_primed_coefficients(ensemble, tables):
    samples, _, _ = tables
    P = CoefficientPredictor(ensemble)
    rec = samples.iloc[3].to_dict()
    assert 'recovery' not in rec
    gmax = rec['final_product_stoich_mass_ratio']
    gu, gu2 = 0.5 * gmax, 0.7 * gmax
    y, t = [0.3, 0.5], [40., 100.]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p0, p1 = P.predict(rec), P.predict(rec, recovery=gu)
        m0, m1 = P.mpsp(rec, y=y, t=t), P.mpsp(rec, y=y, t=t, recovery=gu)
        u1 = P.unprimed(rec, recovery=gu)
    keep = [c for c in p0 if c[0] in 'abcd' and (len(c) == 1 or c[1] == '_')] + \
        [c for c in p0 if c.startswith(('ch__', 'int__'))] + ['ad_ratio']
    np.testing.assert_allclose(p1[keep].to_numpy(float), p0[keep].to_numpy(float), rtol=1e-12)
    assert p0['g_source'].iloc[0] == 'predicted' and np.isclose(p0['g_predicted'].iloc[0], p0['g'].iloc[0])
    assert p1['g_source'].iloc[0] == 'user' and p1['g'].iloc[0] == gu and p1['g_std'].iloc[0] == 0.
    assert np.isclose(p1['g_predicted'].iloc[0], p0['g'].iloc[0])
    assert p1['g_lo'].iloc[0] == p1['g_hi'].iloc[0] == gu
    np.testing.assert_allclose(u1.to_numpy(float)[0], p0[list('abcd')].to_numpy(float)[0] / gu, rtol=1e-12)
    # MPSP = numerator of the predicted primed coefficients / user recovery
    a, b, c, d = (p0[k].iloc[0] for k in 'abcd')
    Y, T = m1['y'].to_numpy(), m1['t'].to_numpy()
    np.testing.assert_allclose(m1['MPSP'], (a + b / Y + c / T + d / (Y * T)) / gu, rtol=1e-9)
    assert (m1['g_source'] == 'user').all() and (m0['g_source'] == 'predicted').all()
    assert (m1['MPSP_lo'] <= m1['MPSP']).all() and (m1['MPSP'] <= m1['MPSP_hi']).all()
    assert not np.allclose(m0['MPSP'], m1['MPSP'])
    # spec field 'recovery' (the argument takes precedence), per-row values, None/NaN = predicted
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p2 = P.predict(dict(rec, recovery=gu))
        p3 = P.predict(dict(rec, recovery=gu), recovery=gu2)
        p4 = P.predict([rec, rec, rec], recovery=[None, gu, float('nan')])
        p5 = P.predict([dict(rec, recovery=None), dict(rec, recovery=gu)])
        p6 = P.predict(rec, recovery=gmax)                         # bound itself is allowed
    assert p2['g'].iloc[0] == gu and p3['g'].iloc[0] == gu2
    assert list(p4['g_source']) == ['predicted', 'user', 'predicted'] and p4['g'].iloc[1] == gu
    assert np.isclose(p4['g'].iloc[0], p0['g'].iloc[0]) and np.isclose(p4['g'].iloc[2], p0['g'].iloc[0])
    assert list(p5['g_source']) == ['predicted', 'user'] and p6['g'].iloc[0] == gmax


@pytest.mark.parametrize('bad', [0., -0.2, 'gmax*1.01', float('inf'), 'abc', [0.5, 0.5]])
def test_user_recovery_validation(ensemble, tables, bad):
    samples, _, _ = tables
    P = CoefficientPredictor(ensemble)
    rec = samples.iloc[3].to_dict()
    if bad == 'gmax*1.01':
        bad = 1.01 * rec['final_product_stoich_mass_ratio']
    match = 'one value per spec row' if isinstance(bad, list) else 'g_max'
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with pytest.raises(ValueError, match=match):
            P.predict(rec, recovery=bad)
        with pytest.raises(ValueError, match=match):
            P.mpsp(rec, y=0.3, t=50., recovery=bad)
        if not isinstance(bad, list):
            with pytest.raises(ValueError, match="spec field 'recovery'"):
                P.predict(dict(rec, recovery=bad))
    assert 'recovery' in P.field_units


def test_unused_regime_fields_not_reported_as_ambiguous(ensemble):
    import copy
    ens = copy.copy(ensemble)
    ens.spec = copy.deepcopy(ensemble.spec)
    P = CoefficientPredictor(ens)
    assert not P.spec.downstream
    reg = P.spec.spec_registry['regime']['HP|AcrylicAcid']
    reg['varying'] = ['separation_principle', 'neutralization', 'separation_ops', 'water_removal_route']
    for c in reg['varying']:
        reg['fields'].pop(c, None)
    q = dict(feedstock='corn', fermentation_product='HP', final_product='AcrylicAcid', r=0.5, feedstock_price=0.2)
    _, det = P._frame(q, details=True)
    assert det[0]['ambiguous'] == ['separation_principle', 'neutralization']
    with pytest.warns(UserWarning, match='not unique'):
        P.predict(dict(q, recovery=None))
