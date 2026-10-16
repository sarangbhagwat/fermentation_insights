# -*- coding: utf-8 -*-
"""Metrics, baselines and the CV runner (sklearn baselines only: fast)."""
import numpy as np

from fermentation_insights.coeff_predictor import evaluate as E


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


def test_ablation_config_parsing():
    m, l = E._configs('pinn-no_phys+no_channels')
    assert not m.channels and not l.physics and not l.collocation
    m, l = E._configs('plain')
    assert m.kind == 'plain' and not l.physics
