# -*- coding: utf-8 -*-
"""Metrics, baselines and the CV runner (sklearn baselines only: fast)."""
import numpy as np
import pandas as pd
import pytest

from fermentation_insights.coeff_predictor import evaluate as E
from fermentation_insights.coeff_predictor import train as T


def test_metrics_of_oracle_predictions(tables):
    samples, points, _ = tables
    pred = dict(coef=samples[['a', 'b', 'c', 'd']].to_numpy(float), g=samples['g'].to_numpy(float),
                g_max=samples['final_product_stoich_mass_ratio'].to_numpy(float))
    m = {(a, b): v for a, b, v in E.compute_metrics(pred, samples, points)}
    for k in ('a', 'b', 'c', 'd', 'g'):
        assert np.isclose(m[('R2', k)], 1.0) and np.isclose(m[('MAE', k)], 0.0)
    assert np.isclose(m[('mpsp_mape', 'MPSP')], m[('oracle_mape', 'MPSP')])
    assert np.isclose(m[('oracle_mape', 'MPSP')], samples['fit_mape'].mean(), rtol=0.2)
    for k in ('mono_y', 'mono_t', 'positive', 'g_le_gmax'):
        assert m[('violation_rate', k)] == 0.0


def test_baselines_and_cv_runner(tables):
    samples, points, _ = tables
    res = E.run_cv(samples, points, 'config', methods=['ridge', 'knn', 'hgb'], max_folds=2, log=None)
    assert not (res['metric'] == 'error').any(), res[res['metric'] == 'error']
    assert set(res['method']) == {'ridge', 'knn', 'hgb'} and set(res['fold']) == {0, 1, -1}
    pooled = res[(res.fold == -1) & (res.metric == 'R2')]
    assert len(pooled) == 15 and (res[res.fold == -1]['label'] == 'pooled').all()
    assert {'scheme', 'fold', 'label', 'method', 'metric', 'target', 'value'} <= set(res.columns)
    v = res[(res.metric == 'violation_rate') & (res.target == 'mono_r')]
    assert len(v) == 6 and v['value'].between(0, 1).all()
    s = E.summarize(res)
    assert ('config', 'mpsp_mape', 'MPSP') in s.index
    # baselines fit the same subset the networks would (calibration and
    # early-stopping hold-outs excluded) and report it
    f0 = res[(res.fold == 0) & (res.metric == 'R2')]
    assert (f0['n_fit'] < f0['n_train']).all() and (f0['n_calib'] > 0).all()
    assert (s.loc[('config', 'n_folds_ok', 'folds')] == 2).all()


def test_ablation_config_parsing():
    m, l = E._configs('pinn-no_phys+no_channels')
    assert not m.channels and not l.physics and not l.collocation
    m, l = E._configs('plain')
    assert m.kind == 'plain' and not l.physics


def test_summarize_keeps_final_test_and_uses_common_folds():
    rows = [dict(scheme='random', fold=0, method='ridge', metric='R2', target='a', value=0.5),
            dict(scheme='random', fold=-1, method='ridge', metric='R2', target='a', value=0.6),
            dict(scheme='final_test', fold=0, method='ridge', metric='R2', target='a', value=0.7),
            dict(scheme='final_test', fold=0, method='ridge', metric='mpsp_mape', target='MPSP', value=0.1),
            dict(scheme='random', fold=0, method='ridge', metric='mpsp_mape', target='MPSP', value=0.1),
            dict(scheme='random', fold=1, method='ridge', metric='mpsp_mape', target='MPSP', value=0.3),
            dict(scheme='random', fold=0, method='knn', metric='mpsp_mape', target='MPSP', value=0.2),
            dict(scheme='random', fold=1, method='knn', metric='error', target='boom', value=np.nan)]
    with pytest.warns(UserWarning, match='excluded'):
        s = E.summarize(pd.DataFrame(rows))
    assert s.loc[('final_test', 'R2', 'a'), 'ridge'] == 0.7
    assert s.loc[('random', 'R2', 'a'), 'ridge'] == 0.6
    assert s.loc[('random', 'mpsp_mape', 'MPSP'), 'ridge'] == 0.1   # fold 1 dropped: knn failed there
    assert s.loc[('random', 'n_folds_ok', 'folds'), 'knn'] == 1
    assert s.loc[('random', 'n_folds_ok', 'folds'), 'ridge'] == 2


def test_grouped_holdout_fit_subsets_and_empty_splits(tables):
    samples, points, _ = tables
    tr, te = E.grouped_holdout(samples, 0.2)
    assert len(te) and not set(samples.config_id.iloc[tr]) & set(samples.config_id.iloc[te])
    one = samples[samples.config_id == samples.config_id.iloc[0]].reset_index(drop=True)
    with pytest.raises(ValueError):
        E.grouped_holdout(one, 0.2)
    with pytest.raises(SystemExit):
        E.main(['--processed', '.', '--report-final'])
    m = E.fit_method('ridge', samples, points, group_col='config_id', seed=3)
    sp = T.ensemble_splits(samples, 1, 0.15, T.TrainConfig().val_frac, 'config_id', 3)
    assert m.fit_info['n_fit'] == len(sp['members'][0][0]) and m.fit_info['n_calib'] == len(sp['calib'])
    res = E.run_cv(one, points, 'config', methods=['ridge'], log=None)
    assert len(res) == 0       # the only split has an empty training set: skipped, not an error
