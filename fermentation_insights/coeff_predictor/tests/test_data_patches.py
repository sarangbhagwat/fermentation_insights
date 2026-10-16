# -*- coding: utf-8 -*-
"""Post-hoc data corrections (data_patches) on small synthetic records, and the
raw-reader behaviour they rely on (per-attempt config record, M11)."""
import json
import shutil

import numpy as np
import pandas as pd
import pytest

from fermentation_insights.coeff_predictor import data_patches as dp
from fermentation_insights.coeff_predictor import dataset, synthetic


def _tables(samples):
    """samples: list of (bio, cid, harness_version, config_extra, start_extra,
    [point dicts]) -> (cfg, starts, points) as dataset.read_raw returns them."""
    cfgs, starts, points = {}, [], []
    for k, (bio, cid, hv, cextra, sextra, pts) in enumerate(samples):
        key = (bio, cid, json.dumps(cextra, sort_keys=True))  # one record per version
        crec = cfgs.setdefault(key, dict(biorefinery=bio, config_id=cid, **cextra))
        sid = f'{bio}|{cid}|s1|{k}'
        starts.append(dict(sample_id=sid, harness_version=hv, _config=crec, _complete=True, **sextra))
        for i, p in enumerate(pts):
            points.append(dict(dict(biorefinery=bio, config_id=cid, sample_id=sid, point_id=i,
                                    status='ok'), **p))
    return pd.DataFrame(list(cfgs.values())), pd.DataFrame(starts), pd.DataFrame(points)


def _ng(fu, cu, fm, cm, ngf, ngp=0.3, **kw):
    return dict(ng_flow_utility=fu, ng_cost_utility=cu, ng_flow_material=fm, ng_cost_material=cm,
                natural_gas_flow=ngf, natural_gas_price=ngp, y=0.3, t=50., **kw)


def _apply_twice(cfg, starts, points):
    c1, s1, p1, n1 = dp.apply(cfg, starts, points)
    c2, s2, p2, n2 = dp.apply(c1, s1, p1)
    num = p1.select_dtypes('number').columns
    np.testing.assert_allclose(p2[num].to_numpy(float), p1[num].to_numpy(float), equal_nan=True)
    assert all(v.get('points', 0) == 0 and v.get('config_records', 0) == 0 for v in n2.values()), n2
    assert c2.drop(columns=[dp.MARK], errors='ignore').equals(c1.drop(columns=[dp.MARK], errors='ignore'))
    return c1, s1, p1, n1


def test_D3_natural_gas_pools_per_adapter():
    S = [
        ('HP', 'HP_corn', 2, {}, {}, [_ng(0., 0., 100., 0., 100.)]),
        ('oxalic', 'oxalic_neutral_sugarcane', 2, {}, {},
         [_ng(0., 0., 50., 0., 50., x_BT_natural_gas_cost=15.)]),
        ('isobutanol', 'IBO_EtOH_corn', 2, {}, {}, [_ng(30., 9., 70., 0., 100.)]),
        ('TAL', 'TAL_corn', 2, {}, {}, [_ng(20., 6., 50., 0., 70., x_natural_gas_BT=50.)]),
        ('ethanol', 'ethanol_cornstover', 2, {}, {}, [_ng(0., 0., 40., 0., 40., x_natural_gas_BT=40.)]),
        ('succinic', 'succinic_corn', 2, {}, {},
         [_ng(0., 0., 70., 2.5, 70., x_natural_gas_BT=60., x_natural_gas_drying=10.)]),
        # not affected: lactic (priced boiler gas), ethanol corn, harness v3 rows
        ('lactic', 'lactic_neutral_SHF_cornstover', 2, {}, {}, [_ng(0., 0., 80., 20., 80.)]),
        ('ethanol', 'ethanol_corn', 2, {'adapter_version': 2}, {}, [_ng(10., 3., 0., 0., 10.)]),
        ('HP', 'HP_glucose', 3, {}, {}, [_ng(100., 30., 0., 0., 100.)]),
        # v2 row without the defect's signature (material NG priced): skipped
        ('HP', 'HP_sugarcane', 2, {}, {}, [_ng(0., 0., 100., 25., 100.)]),
    ]
    cfg, starts, points = _tables(S)
    _, _, p, n = _apply_twice(cfg, starts, points)
    exp = {  # config -> (fu, cu, fm, cm)
        'HP_corn': (100., 30., 0., 0.), 'oxalic_neutral_sugarcane': (50., 15., 0., 0.),
        'IBO_EtOH_corn': (100., 30., 0., 0.), 'TAL_corn': (70., 21., 0., 0.),
        'ethanol_cornstover': (40., 12., 0., 0.), 'succinic_corn': (60., 18., 10., 2.5),
        'lactic_neutral_SHF_cornstover': (0., 0., 80., 20.), 'ethanol_corn': (10., 3., 0., 0.),
        'HP_glucose': (100., 30., 0., 0.), 'HP_sugarcane': (0., 0., 100., 25.),
    }
    for cid, e in exp.items():
        r = p[p['config_id'] == cid].iloc[0]
        got = tuple(r[c] for c in dp_cols())
        np.testing.assert_allclose(got, e, err_msg=cid)
        # the total NG cost never changes except by the boiler gas now charged
    assert n['D3']['points'] == 6 and n['D3']['points_skipped'] == 1
    assert 'D3' in p.loc[p['config_id'] == 'HP_corn', dp.MARK].iloc[0]
    assert set(n['D3']['by_config']) == {'HP|HP_corn', 'oxalic|oxalic_neutral_sugarcane',
                                         'isobutanol|IBO_EtOH_corn', 'TAL|TAL_corn',
                                         'ethanol|ethanol_cornstover', 'succinic|succinic_corn'}


def dp_cols():
    return ('ng_flow_utility', 'ng_cost_utility', 'ng_flow_material', 'ng_cost_material')


def _cane_points(credit=0.2593, n=4, vary=False):
    pts = []
    for i in range(n):
        h, Fs, Ff = 8000., 1e4 * (1 + 0.1 * i), 0.97e4 * (1 + 0.1 * i)
        c = credit * (1 + (0.01 * i if vary else 0.))
        pts.append(dict(y=0.4, t=60. + i, MPSP=0.7 + 0.01 * i, operating_hours=h, F_final_product=Ff,
                        x_F_product_stream=Fs, coproduct_sales=c * h * Fs,
                        x_MPSP_denatured=0.66, x_MESP_USD_per_gal_denatured=0.66 * dp.KG_PER_GAL_ETHANOL))
    return pts


def test_D4_RIN_credit():
    cfg, starts, points = _tables([('ethanol', 'ethanol_sugarcane', 2, {}, {}, _cane_points()),
                                   ('ethanol', 'ethanol_sugarcane', 2, {}, {}, _cane_points())])
    old = points.copy()
    _, _, p, n = _apply_twice(cfg, starts, points)
    np.testing.assert_allclose(p['MPSP'], old['MPSP'] + old['coproduct_sales']
                               / (old['operating_hours'] * old['F_final_product']))
    np.testing.assert_allclose(p['x_MPSP_denatured'], 0.66 + 0.2593)
    np.testing.assert_allclose(p['x_MESP_USD_per_gal_denatured'], (0.66 + 0.2593) * dp.KG_PER_GAL_ETHANOL)
    np.testing.assert_allclose(p['x_RIN_credit'], old['coproduct_sales'] / old['operating_hours'])
    assert (p['coproduct_sales'] == 0).all() and n['D4']['points'] == 8
    assert n['D4']['credit_USD_per_kg_product_stream'] == pytest.approx(0.2593)
    # credit per kg not constant -> coproduct sales are not the credit alone: no patch
    cfg, starts, points = _tables([('ethanol', 'ethanol_sugarcane', 2, {}, {}, _cane_points(vary=True))])
    _, _, p, n = dp.apply(cfg, starts, points)
    np.testing.assert_allclose(p['MPSP'], points['MPSP'])
    assert n['D4']['points'] == 0 and n['D4']['points_skipped'] == 4
    # adapter v2 rows carry the fix
    cfg, starts, points = _tables([('ethanol', 'ethanol_sugarcane', 3, {'adapter_version': 2}, {},
                                    _cane_points())])
    _, _, p, n = dp.apply(cfg, starts, points)
    np.testing.assert_allclose(p['MPSP'], points['MPSP'])


def _hp_points(p_cell, ys, xcc=None):
    infl = (1 - p_cell) * 1.0
    return [dict(y=y, t=50., x_lime_fermentation=1000., base_cost=286.335,
                 x_cell_mass_conversion=(xcc if xcc is not None else p_cell) if y < infl else 0.01)
            for y in ys]


def test_D6_D7_HP():
    ys = [0.1, 0.4, 0.70, 0.74, 0.78]
    cfg, starts, points = _tables([
        ('HP', 'HP_corn', 2, dict(theoretical_yield=1.0, inflection_yield=0.755),
         dict(p_base_price=0.2, p_cell_mass_yield=0.28), _hp_points(0.28, ys)),
        ('HP', 'HP_corn', 2, dict(theoretical_yield=1.0, inflection_yield=0.755),
         dict(p_base_price=0.3, p_cell_mass_yield=0.22), _hp_points(0.22, ys, xcc=0.3)),  # basis mismatch
    ])
    _, _, p, n = _apply_twice(cfg, starts, points)
    np.testing.assert_allclose(p['base_cost'], [200.] * 5 + [300.] * 5)
    a = p[p['sample_id'].str.endswith('|0')]
    np.testing.assert_allclose(a['x_inflection_yield'], 0.72)
    assert a['x_above_inflection'].tolist() == [False, False, False, True, True]
    b = p[p['sample_id'].str.endswith('|1')]
    assert b['x_inflection_yield'].isna().all()       # skipped: cell-mass basis does not match
    assert n['D6']['points'] == 10 and n['D7']['points'] == 5 and n['D7']['points_skipped'] == 5
    assert n['D7']['points_above_inflection'] == 2


def test_config_patches_and_unpatchable():
    tal_params = json.dumps([['p_IRR', 0.1, ''],
                             [dp.TAL_OLD_CELL_MASS, 0.339, 'Fermentation: y. lipolytica yield']])
    cfg, starts, points = _tables([
        ('TAL', 'TAL_corn', 2, dict(neutralization=False, neutralizing_agent='',
                                    product_water_solubility=5.41, parameters=tal_params),
         {dp.TAL_OLD_CELL_MASS: 0.35}, [dict(y=0.2, t=20.)]),
        ('succinic', 'succinic_corn', 2, dict(neutralization=False, neutralizing_agent=''), {},
         [dict(y=0.2, t=20.)]),
        ('HP', 'HP_corn', 2, dict(neutralization=False,
                                  separation_principle='acidulation+adsorption|dehydration'), {},
         [dict(y=0.2, t=20.)]),
        ('HP', 'HP_neutral_corn', 2, dict(neutralization=True,
                                          separation_principle='acidulation+adsorption|dehydration'), {},
         [dict(y=0.2, t=20.)]),
        ('TAL', 'TAL_glucose', 3, dict(adapter_version=2, neutralizing_agent='NaOH',
                                       product_water_solubility=9.8), {}, [dict(y=0.2, t=20.)]),
        ('ethanol', 'ethanol_corn', 2, {}, {}, [dict(y=0.2, t=20.)]),
        ('ethanol', 'ethanol_corn', 3, {'adapter_version': 2}, {}, [dict(y=0.2, t=21.)]),
    ])
    c, s, p, n = _apply_twice(cfg, starts, points)
    C = c.set_index('config_id')
    assert C.loc['TAL_corn', 'neutralizing_agent'] == 'NaOH' and not C.loc['TAL_corn', 'neutralization']
    assert C.loc['TAL_corn', 'product_water_solubility'] == 9.8
    assert [q[0] for q in json.loads(C.loc['TAL_corn', 'parameters'])] == ['p_IRR', 'p_cell_mass_yield']
    assert C.loc['succinic_corn', 'neutralizing_agent'] == 'Ca(OH)2'
    assert C.loc['HP_corn', 'separation_principle'] == 'adsorption|dehydration'
    assert C.loc['HP_neutral_corn', 'separation_principle'] == 'acidulation+adsorption|dehydration'
    S = s.set_index('sample_id')
    t0 = 'TAL|TAL_corn|s1|0'
    assert S.loc[t0, 'p_cell_mass_yield'] == 0.35 and np.isnan(S.loc[t0, dp.TAL_OLD_CELL_MASS])
    assert S.loc[t0, '_config']['neutralizing_agent'] == 'NaOH'    # the sample's own record too
    assert set(S.loc[t0, dp.MARK].split(';')) == {'M1', 'M3', 'M12'}
    assert S['_unpatchable'].tolist().count('D2') == 1
    assert S.loc['ethanol|ethanol_corn|s1|5', '_unpatchable'] == 'D2'
    assert n['M1']['configs'] == ['TAL|TAL_corn', 'succinic|succinic_corn']
    assert n['M23']['configs'] == ['HP|HP_corn'] and n['M3']['samples'] == 1
    assert n['D2']['samples_unpatchable'] == 1


def _write(path, recs):
    with open(path, 'w') as fh:
        fh.writelines(json.dumps(r) + '\n' for r in recs)


def test_read_raw_attempt_config_and_resumed_duplicates(tmp_path):
    """M11: a resumed sample appends a second attempt; only the complete
    attempt's points are kept. Each attempt uses the config record preceding
    it in its file (a resumed campaign of a newer adapter appends its own)."""
    cid = 'HP_corn'
    c1 = dict(record='config', biorefinery='HP', config_id=cid, theoretical_yield=1.0)
    c2 = dict(c1, adapter_version=2)

    def sample(k, n, end=True, hv=2):
        sid = f'HP|{cid}|s1|{k}'
        out = [dict(record='sample_start', sample_id=sid, harness_version=hv, p_base_price=0.2)]
        out += [dict(record='point', biorefinery='HP', config_id=cid, sample_id=sid, point_id=i,
                     y=0.1 + 0.1 * i, t=50., status='ok', base_cost=1.0, x_lime_fermentation=10.)
                for i in range(n)]
        if end:
            out.append(dict(record='sample_end', sample_id=sid, n_ok=n))
        return out
    _write(tmp_path / 'HP_corn.jsonl', [c1] + sample(0, 2, end=False) + [c1] + sample(0, 3)
           + [c2] + sample(1, 4, hv=3))
    cfg, starts, points, stats = dataset.read_raw([str(tmp_path)])
    assert stats['restarted_attempts'] == 1
    assert not points.duplicated(['sample_id', 'point_id']).any()
    assert (points['sample_id'] == f'HP|{cid}|s1|0').sum() == 3
    S = starts.set_index('sample_id')
    assert 'adapter_version' not in S.loc[f'HP|{cid}|s1|0', '_config']
    assert S.loc[f'HP|{cid}|s1|1', '_config']['adapter_version'] == 2
    _, _, p, n = dp.apply(cfg, starts, points)
    assert n['D6']['points'] == 3                     # adapter v1 sample only
    np.testing.assert_allclose(p.loc[p['sample_id'].str.endswith('|0'), 'base_cost'], 2.0)
    np.testing.assert_allclose(p.loc[p['sample_id'].str.endswith('|1'), 'base_cost'], 1.0)


def test_build_applies_patches_and_drops_flagged_points(tmp_path):
    """build() reports patch counts, rejects unpatchable samples and drops
    points flagged x_above_inflection (and uses x_inflection_yield)."""
    d = tmp_path / 'raw'
    synthetic.generate(str(d), names=['SYNB'], n_random=2, n_productivity=1)
    pa = next(d.rglob('*.jsonl'))
    recs = [json.loads(line) for line in open(pa)]
    sids = [r['sample_id'] for r in recs if r['record'] == 'sample_start']
    flagged = []
    for r in recs:
        if r.get('sample_id') == sids[0] and r['record'] == 'point' and r.get('status') == 'ok' \
                and len(flagged) < 3:
            r['x_above_inflection'] = True
            flagged.append((r['y'], r['t']))
    _write(pa, recs)
    base = dataset.build([str(d)], patches=False)
    data = dataset.build([str(d)])
    assert data.report['patches'] and 'D3' in data.report['data_patches']
    assert data.report['data_patches']['D3']['points'] == 0     # synthetic: no affected adapter
    assert data.report['points_beyond_inflection_dropped'] >= 3
    P0 = data.points[data.points['sample_id'] == sids[0]]
    assert not any(((P0['y'] == y) & (P0['t'] == t)).any() for y, t in flagged)
    assert len(data.samples) == len(base.samples)
    assert {'harness_version', 'adapter_version', 'data_patches'} <= set(data.samples.columns)
    # an adapter-v1 ethanol_corn configuration is not patchable (D2): rejected
    d2 = tmp_path / 'raw2'
    shutil.copytree(d, d2)
    pb = next(d2.rglob('*.jsonl'))
    cid0 = next(r['config_id'] for r in recs if r['record'] == 'config')
    out = []
    for r in recs:
        if r.get('config_id') == cid0 or str(r.get('sample_id', '')).startswith(f'SYNB|{cid0}|'):
            r = dict(r)
            if 'biorefinery' in r:
                r['biorefinery'] = 'ethanol'
            if 'config_id' in r:
                r['config_id'] = 'ethanol_corn'
            if 'sample_id' in r:
                r['sample_id'] = r['sample_id'].replace(f'SYNB|{cid0}|', 'ethanol|ethanol_corn|')
        out.append(r)
    _write(pb, out)
    data2 = dataset.build([str(d2)])
    n_corn = sum(1 for r in out if r['record'] == 'sample_start'
                 and r['sample_id'].startswith('ethanol|ethanol_corn|'))
    assert data2.report['rejections'].get('unpatchable_D2') == n_corn > 0
    assert not (data2.samples['config_id'] == 'ethanol_corn').any()
