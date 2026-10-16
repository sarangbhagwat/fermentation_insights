# -*- coding: utf-8 -*-
"""End-to-end tests: synthetic raw JSONL (real harness) -> dataset.build ->
recovery of the known ground truth; loader robustness."""
import json
import os
import shutil

import numpy as np
import pandas as pd
import pytest

from fermentation_insights.coeff_predictor import dataset, physics, synthetic


@pytest.fixture(scope='module')
def raw(tmp_path_factory):
    d = tmp_path_factory.mktemp('raw')
    paths = synthetic.generate(str(d), n_random=3, n_productivity=1, fail_every=17)
    return str(d), paths


@pytest.fixture(scope='module')
def built(raw, tmp_path_factory):
    d, paths = raw
    out = str(tmp_path_factory.mktemp('proc'))
    data = dataset.build([d], out_dir=out)
    return data, out, synthetic.ground_truth(paths)


def test_synthetic_round_trip(built):
    data, out, gt = built
    s = data.samples
    assert len(s) == len(gt) > 0
    m = s.merge(gt, on='sample_id', suffixes=('', '_gt'))
    assert len(m) == len(gt)
    np.testing.assert_allclose(m['phi'], m['phi_true'], rtol=1e-7)
    np.testing.assert_allclose(m['F_fixed'], m['F_fixed_true'], rtol=1e-6)
    # synthetic FOC = F + 0.03*FCI -> FOC/FCI intercept equals F
    np.testing.assert_allclose(m['diag_F_foc_intercept'], m['F_fixed_true'], rtol=1e-6)
    np.testing.assert_allclose(m['g'], m['g_true'], rtol=1e-12)
    for col in dataset.CHANNEL_TARGETS + dataset.INTENSITY_TARGETS:
        scale = np.abs(m[['a', 'b', 'c', 'd']]).max(axis=1) if col.startswith('ch__') \
            else np.abs(m[col + '_gt']).max() + 1e-12
        err = np.abs(m[col] - m[col + '_gt']) / scale
        assert err.max() < 1e-7, (col, err.max())
    # contract columns present
    for c in ('sample_id', 'config_id', 'biorefinery', 'sample_kind', 'seed', 'feedstock',
              'feedstock_type', 'fermentation_product', 'final_product', 'neutralization',
              'final_product_stoich_mass_ratio', 'inflection_yield', 'y_min', 't_max',
              'baseline_productivity', 'r', 'r_mult', 'feedstock_price', 'electricity_price',
              'natural_gas_price', 'a', 'b', 'c', 'd', 'g', 'R2', 'n_points', 'fit_mape',
              'phi', 'F_fixed', 'closure_rel', 'p_IRR', 'bp_capex_scale', 'fp_MW', 'fs_price'):
        assert c in s.columns, c
    assert set(s['sample_kind']) == {'baseline', 'productivity', 'random'}
    assert (s['sum_check'] < 1e-10).all()
    assert s['channels_ok'].all()


def test_saved_tables_reload(built):
    data, out, _ = built
    for n in ('samples', 'points', 'configs'):
        assert os.path.exists(os.path.join(out, f'{n}.csv.gz'))
    re = dataset.load(out)
    assert len(re.samples) == len(data.samples)
    assert len(re.points) == len(data.points) == data.samples['n_points'].sum()
    assert set(re.points['sample_id']) == set(re.samples['sample_id'])
    assert {'y', 't', 'MPSP', 'recovery', 'MPSP_fit', 'm__feedstock', 'm__closure'} <= set(re.points.columns)
    # channel sums reproduce MPSP at the point level after the CSV round trip
    M = re.points[[f'm__{c}' for c in physics.CHANNELS]].sum(axis=1)
    np.testing.assert_allclose(M, re.points['MPSP'], rtol=1e-9)
    assert len(re.configs) == 5 and (re.configs['n_samples_retained'] > 0).all()


def test_kappa_nonlinearity_approximate(tmp_path):
    paths = synthetic.generate(str(tmp_path / 'raw'), names=['SYNB'], n_random=2,
                               n_productivity=0, kappa=0.02)
    data = dataset.build([str(tmp_path / 'raw')])
    gt = synthetic.ground_truth(paths, kappa=0.02)
    m = data.samples.merge(gt, on='sample_id', suffixes=('', '_gt'))
    assert len(m) > 0
    assert (np.abs(m['phi'] / m['phi_true'] - 1) < 0.15).all(), m['phi'] / m['phi_true']
    assert (m['closure_rel'] < 0.03).all() and (m['closure_rel'] > 0).all()
    assert (m['sum_check'] < 1e-10).all()


def _rewrite(path, fn):
    with open(path) as fh:
        lines = fh.readlines()
    with open(path, 'w') as fh:
        fh.writelines(fn(lines))


def test_loader_robustness(raw, tmp_path):
    src, _ = raw
    d = tmp_path / 'raw'
    shutil.copytree(src, d)
    base = dataset.build([src])
    # 1) truncated last line + garbage line
    pa = d / 'SYNA' / 'SYNA.jsonl'
    _rewrite(pa, lambda L: L[:5] + ['{"record": "point", "sampl\n'] + L[5:] + [L[-1][:40]])
    # 2) columns missing in one adapter's points
    pb = d / 'SYNB' / 'SYNB.jsonl'

    def drop(L):
        out = []
        for line in L:
            rec = json.loads(line)
            if rec.get('record') == 'point':
                for k in ('natural_gas_flow', 'base_cost', 'chp__material_cost', 'GWP'):
                    rec.pop(k, None)
            out.append(json.dumps(rec) + '\n')
        return out
    _rewrite(pb, drop)
    # 3) an incomplete (still running) sample at the end of SYNC
    pc = d / 'SYNC' / 'SYNC.jsonl'
    lines = open(pc).readlines()
    end_idx = max(i for i, l in enumerate(lines) if '"sample_end"' in l)
    sid = json.loads(lines[end_idx])['sample_id']
    _rewrite(pc, lambda L: L[:end_idx] + L[end_idx + 1:])
    # 4) a smoke file duplicating a SYNA sample with a restarted attempt
    sa = [json.loads(l) for l in open(src + '/SYNA/SYNA.jsonl')]
    sid_dup = next(r['sample_id'] for r in sa if r['record'] == 'sample_start')
    dup = [r for r in sa if r.get('sample_id') == sid_dup]
    with open(d / 'SYNA' / 'SYNA_smoke.jsonl', 'w') as fh:
        fh.write(json.dumps(dup[0]) + '\n' + json.dumps(dup[1]) + '\n')  # aborted attempt
        for r in dup:
            fh.write(json.dumps(r) + '\n')
    data = dataset.build([str(d)])
    st = data.report['stats']
    assert st['bad_lines'] >= 2
    assert st['duplicate_samples'] == 1 and st['restarted_attempts'] >= 1
    assert data.report['rejections'].get('incomplete', 0) == 1
    assert sid not in set(data.samples['sample_id'])
    assert len(data.samples) == len(base.samples) - 1
    m = data.samples.merge(base.samples, on='sample_id', suffixes=('', '_0'))
    np.testing.assert_allclose(m[['a', 'b', 'c', 'd', 'g']].to_numpy(),
                               m[['a_0', 'b_0', 'c_0', 'd_0', 'g_0']].to_numpy(), rtol=1e-12)
    synb = data.samples[data.samples['biorefinery'] == 'SYNB']
    assert synb['diag_imputed'].str.contains('natural_gas_flow').all()
    assert np.allclose(synb[[f'ch__natural_gas__{k}' for k in 'abcd']], 0)
    assert (synb['sum_check'] < 1e-10).all()
    # include_incomplete keeps the partial sample if it passes the filters
    data2 = dataset.build([str(d)], filters=dataset.QualityFilters(require_complete=False))
    assert len(data2.samples) >= len(data.samples)


def test_quality_filters(raw):
    src, _ = raw
    strict = dataset.build([src], filters=dataset.QualityFilters(min_points=30))
    assert len(strict.samples) == 0
    assert strict.report['rejections']['too_few_points'] == strict.report['stats']['samples_seen']
    data = dataset.build([src], filters=dataset.QualityFilters(max_closure_rel=-1))
    assert len(data.samples) == 0 and data.report['rejections']['high_closure'] > 0


def test_empty_dir(tmp_path):
    data = dataset.build([str(tmp_path)])
    assert len(data.samples) == 0 and data.report['n_retained'] == 0


def test_sample_level_inflection_config_ranking_and_fs(raw, tmp_path):
    src, _ = raw
    d = tmp_path / 'raw'
    shutil.copytree(src, d)
    pa = d / 'SYNA' / 'SYNA.jsonl'
    recs = [json.loads(l) for l in open(pa)]
    cfg0 = next(r for r in recs if r['record'] == 'config')
    sids = [r['sample_id'] for r in recs if r['record'] == 'sample_start'
            and r['sample_id'].split('|')[1] == cfg0['config_id']]
    sid0, sid1 = sids[0], sids[1]
    for r in recs:
        if r.get('sample_id') == sid1 and r['record'] == 'point':
            r['x_inflection_yield'] = 0.7          # the sample's own (lower) inflection
        if r.get('sample_id') == sid1 and r['record'] == 'sample_start':
            r['fs_w_fermentable_sugar_eq'] = 0.5   # feedstock described after sampling
    with open(pa, 'w') as fh:
        fh.writelines(json.dumps(r) + '\n' for r in recs)
    # stale smoke file (sorts after the production file) with a different
    # config record and a smoke-only sample of the same configuration
    stale = dict(cfg0, pKa1=9.99)
    smoke_only = sid0.rsplit('|', 1)[0] + '|999'
    extra = [dict(r, sample_id=smoke_only) for r in recs if r.get('sample_id') == sid0]
    with open(d / 'SYNA' / 'SYNA_smoke.jsonl', 'w') as fh:
        fh.writelines(json.dumps(r) + '\n' for r in [stale] + extra)
    data = dataset.build([str(d)])
    s = data.samples.set_index('sample_id')
    assert data.report['stats']['config_conflicts'] == 1
    assert data.report['rejections'].get('smoke_superseded') == 1 and smoke_only not in s.index
    assert (s.loc[s['config_id'] == cfg0['config_id'], 'pKa1'] == cfg0['pKa1']).all()
    assert s.loc[sid1, 'inflection_yield'] == 0.7
    assert s.loc[sid1, 'inflection_yield_config'] == cfg0['inflection_yield']
    assert (data.points.loc[data.points['sample_id'] == sid1, 'y'] <= 0.7 + 1e-12).all()
    assert s.loc[sid0, 'inflection_yield'] == cfg0['inflection_yield']
    assert s.loc[sid1, 'fs_w_fermentable_sugar_eq'] == 0.5
    assert s.loc[sid0, 'fs_w_fermentable_sugar_eq'] == cfg0['fs_w_fermentable_sugar_eq']
    # g is the median over all ok points (before clipping), as published
    rec1 = [r['recovery'] for r in recs if r.get('sample_id') == sid1 and r['record'] == 'point'
            and r.get('status', 'ok') == 'ok' and r.get('recovery') is not None]
    assert s.loc[sid1, 'g'] == pytest.approx(float(np.median(rec1)), rel=1e-12)


def test_ng_pool_consistent_per_config(built):
    data, _, _ = built
    s = data.samples
    has = s[s['diag_ng_pool'] != 'none']
    assert (has.groupby(['biorefinery', 'config_id'])['diag_ng_pool'].nunique() == 1).all()
    assert s['diag_ng_pool_admissible'].all() and 'diag_ng_pool_vote' in s
    assert (s['diag_fixed_basis'] == 'FOC').all() and not s['diag_F_negative'].any()
