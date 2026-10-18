# -*- coding: utf-8 -*-
"""Report-only statistics of ``report_checks`` on small hand-made inputs:
family-level paired statistics, LOCO subsets, new-product folds, input
twins, the test set without some configurations and the data's own
constraint violations."""
import numpy as np
import pandas as pd
import pytest

from fermentation_insights.coeff_predictor import report_checks as RC


def _res(values):
    """Tidy LOCO results {(label, method): fold MAPE}."""
    return pd.DataFrame([dict(scheme='config', fold=i, label=lab, method=m, metric='mpsp_mape', target='MPSP',
                              value=v) for i, ((lab, m), v) in enumerate(values.items())])


def test_family_stats_average_folds_within_families():
    fam = {'A1': 'A', 'A2': 'A', 'A3': 'A', 'B1': 'B', 'C1': 'C', 'D1': 'D', 'E1': 'E'}
    pin = dict(A1=.10, A2=.10, A3=.10, B1=.20, C1=.30, D1=.10, E1=.10)
    plain = dict(A1=.11, A2=.12, A3=.13, B1=.25, C1=.29, D1=.15, E1=.20)
    res = _res({**{(k, 'pinn'): v for k, v in pin.items()}, **{(k, 'plain'): v for k, v in plain.items()}})
    t = RC.baseline_family_stats(res, fam, methods=('plain',)).iloc[0]
    assert t['n_folds'] == 7 and t['n_families'] == 5
    d = pd.Series({'A': .02, 'B': .05, 'C': -.01, 'D': .05, 'E': .10})
    assert t['mean'] == pytest.approx(d.mean())
    assert (t['n_pinn_better'], t['n_baseline_better']) == (4, 1)
    assert t['dominant_family'] == 'E'


def test_loco_subset_and_single_config_products():
    s = pd.DataFrame(dict(config_id=['a', 'a', 'b', 'c', 'd'], fermentation_product=['P', 'P', 'P', 'Q', 'R']))
    assert RC.single_config_products(s) == ['c', 'd']
    res = _res({('a', 'pinn'): .1, ('b', 'pinn'): .2, ('c', 'pinn'): .9, ('d', 'pinn'): .7,
                ('a', 'knn'): .2, ('b', 'knn'): .4, ('c', 'knn'): .5, ('d', 'knn'): .5})
    t = RC.loco_subset(res, ['c', 'd'], methods=('pinn', 'knn')).set_index('method')
    assert t.loc['pinn', 'n_folds'] == 2 and t.loc['pinn', 'mean'] == pytest.approx(.15)
    assert t.loc['knn', 'diff_mean'] == pytest.approx(.15) and t.loc['knn', 'n_pinn_better'] == 2
    f = RC.loco_per_family(res, {'a': 'F', 'b': 'F', 'c': 'G', 'd': 'H'}).set_index('family')
    assert f.loc['F', 'n_folds'] == 2 and f.loc['F', 'max_config'] == 'b'


class _Spec:
    """Stand-in for a FeatureSpec: encodes x1, x2 (no productivity / prices)."""

    def transform(self, df, strict=False):
        return df[['x1', 'x2']].to_numpy(float)

    def full_matrix(self, X):
        return X


class _P:
    spec = _Spec()


def test_input_twins_and_test_without_them():
    s = pd.DataFrame(dict(config_id=['u', 'u', 'v', 'w'], sample_kind=['random', 'baseline', 'baseline', 'baseline'],
                          sample_index=[1, 0, 0, 0], sample_id=['u1', 'u0', 'v0', 'w0'],
                          x1=[9., 1., 1., 3.], x2=[9., 2., 2., 2.]))
    tw, near = RC.input_twins(_P(), s, {'u': 'test_ood', 'v': 'development', 'w': 'development'})
    assert len(tw) == 1 and set(tw.iloc[0][['config_id', 'twin']]) == {'u', 'v'}
    assert bool(tw.iloc[0]['crosses_split']) and tw.iloc[0]['n_features_different'] == 0
    assert near['distance'] == pytest.approx(2.)
    S = s.set_index('sample_id', drop=False)
    ev = dict(ids=['u0', 'u1', 'v0'], coef=np.arange(12.).reshape(3, 4), g=np.array([.5, .6, .7]),
              scores=dict(config_id=['u', 'u', 'v'], mpsp=[1., 2., 3.]), rows=[('x', 'y', 1.)])
    out = RC.drop_configs(ev, S, ['u'])
    assert out['ids'] == ['v0'] and out['g'].tolist() == [.7] and out['coef'].shape == (1, 4)
    assert out['scores']['mpsp'] == [3.] and 'rows' not in out


def test_data_violation_rates_of_the_fitted_coefficients():
    s = pd.DataFrame(dict(sample_id=['ok', 'bad'], config_id=['A', 'B'], a=[3., 3.], b=[.5, -.5], c=[10., 10.],
                          d=[0., 0.], g=[.9, 1.2], final_product_stoich_mass_ratio=[1., 1.]))
    pts = {sid: pd.DataFrame(dict(sample_id=sid, y=[.2, .4, .2, .4], t=[20., 20., 80., 80.], MPSP=1.))
           for sid in s['sample_id']}
    S = s.set_index('sample_id', drop=False)
    r = RC.data_violation_rates(S, pts, ['ok', 'bad'])
    assert r['viol_mono_y'] == pytest.approx(.5) and r['viol_mono_t'] == 0 and r['viol_g_le_gmax'] == .5
    byc = RC.data_violations_by_config(S, pts, ['ok', 'bad'])
    assert set(zip(byc['constraint'], byc['config_id'])) == {('mono_y', 'B'), ('g_le_gmax', 'B')}


def test_sim_checks_collect_pairs_as_is_and_corrected_runs(tmp_path):
    import json
    from fermentation_insights.coeff_predictor.data_generation import sim_checks as SC

    def pts(scale):
        out = [dict(point='baseline', y=0.5, t=50., status='ok', MPSP=2. * scale, recovery=0.9, biogas_kg_h=10.,
                    ic_mass_gain_kg_h=5.)]
        for i, y in enumerate((0.2, 0.4, 0.6)):
            for j, t in enumerate((20., 40., 60.)):
                m = 0.5 + 0.2 / y + 5. / t + 1. / (y * t)
                out.append(dict(point=f'grid_{i}{j}', y=y, t=t, status='ok', recovery=0.9, biogas_kg_h=10.,
                                ic_mass_gain_kg_h=5., MPSP=m * (scale if y < 0.3 else 1.)))
        return out
    rec = dict(config_id='X', capacity=dict(baseline=100., points=[
        dict(mult=m, p_feedstock_capacity=100. * m, y=0.5, t=50., r=1., status='ok', MPSP=2. / m ** 0.1,
             recovery=0.9, FCI=1.) for m in (1., 2.)]),
        ic_biogas=dict(units=[['U1', 'lumped']], as_is=pts(1.), fixed=pts(1.1)))
    json.dump(rec, open(tmp_path / 'X.json', 'w'))
    json.dump(dict(config_id='Y', ic_biogas=dict(units=[], skipped='none')), open(tmp_path / 'Y.json', 'w'))
    written = SC.collect(str(tmp_path), str(tmp_path / 'csv'))
    names = sorted(p.split('/')[-1] for p in written)
    assert names == ['capacity_response.csv', 'ic_biogas_bias.csv', 'ic_biogas_points.csv']   # no dup check run
    b = pd.read_csv(tmp_path / 'csv' / 'ic_biogas_bias.csv').set_index('config_id')
    assert b.loc['X', 'dMPSP_baseline_pct'] == pytest.approx(10.)
    assert b.loc['X', 'dMPSP_low_yield_pct'] == pytest.approx(10.)
    assert b.loc['X', 'dMPSP_grid_min_pct'] == pytest.approx(0., abs=1e-9)
    assert not b.loc['Y', 'affected']
    c = pd.read_csv(tmp_path / 'csv' / 'capacity_response.csv')
    assert list(c['mult']) == [1., 2.]
