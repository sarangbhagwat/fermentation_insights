# -*- coding: utf-8 -*-
"""Training, splitting, constraint checks, checkpoints and prediction."""
import warnings

import numpy as np
import pandas as pd
import pytest
import torch

from fermentation_insights.coeff_predictor import train as T
from fermentation_insights.coeff_predictor.model import ModelConfig, mpsp_equation_np
from fermentation_insights.coeff_predictor.predict import CoefficientPredictor
from fermentation_insights.coeff_predictor.features import MissingSpecError, ApplicabilityWarning
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
    assert T._conf_q(np.array([1., 2.]), 0.1) == 2.
    assert np.isnan(T._conf_q(np.array([]), 0.1))


def test_plain_and_ablation_training_runs(tables):
    samples, points, _ = tables
    tc = T.TrainConfig(**{**FAST, 'epochs': 8})
    for mcfg, lcfg in ((ModelConfig(kind='plain', width=16, depth=1), T.LossConfig(physics=False, collocation=False)),
                       (ModelConfig(width=16, depth=1, price_explicit=False, channels=False), T.LossConfig())):
        ens = T.train_ensemble(samples, points, 1, mcfg, lcfg, tc, calib_frac=0.1, seed=1)
        batch, outs = ens.predict_samples(samples.iloc[:4])
        assert outs['coef'].shape == (1, 4, 4) and np.isfinite(outs['coef']).all()
