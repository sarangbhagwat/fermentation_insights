"""Frozen final test split of experiments.py."""
import json
import os

import pandas as pd
import pytest

from fermentation_insights.coeff_predictor import experiments as E


def _samples(cfgs, n=40):
    rows = []
    for c in cfgs:
        for i in range(n):
            rows.append(dict(sample_id=f"{c['biorefinery']}|{c['config_id']}|s1|{i}", config_id=c['config_id'],
                             biorefinery=c['biorefinery'], feedstock=c['feedstock']))
    return pd.DataFrame(rows)


def test_planned_list():
    ids = [c['config_id'] for c in E.PLANNED_CONFIGS]
    assert len(ids) == len(set(ids)) == 41
    for c in ('HP_neutral_hexanol_cornstover', 'TAL_SA_corn', 'succinic_neutral_glucose',
              'lactic_lowpH_SHF_cornstover', 'ethanol_corn', 'oxalic_neutral_sugarcane', 'IBO_EtOH_corn'):
        assert c in ids


def test_split_depends_only_on_plan_and_seed():
    a = E.make_frozen_split(2026)
    b = E.make_frozen_split(2026)
    assert a == b
    tc = [c['config_id'] for c in a['test_configs']]
    assert len(tc) == 4                                  # HP, TAL, succinic, lactic (>= 4 planned configs)
    assert {c['biorefinery'] for c in a['test_configs']} == {'HP', 'TAL', 'succinic', 'lactic'}
    assert len({c['feedstock'] for c in a['test_configs']}) == 4    # distinct feedstocks preferred
    # data availability never enters: applying to any subset keeps the same test configurations
    sub = _samples([c for c in E.PLANNED_CONFIGS if c['biorefinery'] in ('HP', 'ethanol')])
    sp, dev = E.apply_split(sub, a)
    assert [c['config_id'] for c in sp['test_configs']] == tc
    assert set(sp['test_configs_missing']) == set(tc) - set(sub['config_id'])
    assert not set(dev['config_id']) & set(tc)


def test_id_membership_stable_under_growth():
    fr = E.make_frozen_split(7)
    cfgs = [c for c in E.PLANNED_CONFIGS if c['biorefinery'] == 'ethanol']
    small, big = _samples(cfgs, 30), _samples(cfgs, 300)
    sp_s, _ = E.apply_split(small, fr)
    sp_b, _ = E.apply_split(big, fr)
    assert set(sp_s['id_test_sample_ids']) <= set(sp_b['id_test_sample_ids'])
    assert set(sp_s['id_test_sample_ids']) == set(sp_b['id_test_sample_ids']) & set(small['sample_id'])
    frac = sp_b['n_test_id'] / len(big)
    assert 0.06 < frac < 0.14


def test_missing_test_config_fails_at_final(tmp_path):
    fr = E.make_frozen_split(2026)
    sub = _samples([c for c in E.PLANNED_CONFIGS if c['biorefinery'] == 'ethanol'])
    with pytest.raises(SystemExit):
        E.apply_split(sub, fr, require_test_configs=True)
    p = tmp_path / 'split.json'
    p.write_text(json.dumps(fr))
    assert E.load_frozen_split(str(p))['test_configs'] == fr['test_configs']
    bad = dict(fr, test_configs=fr['test_configs'][:1])
    p.write_text(json.dumps(bad))
    with pytest.raises(SystemExit):
        E.load_frozen_split(str(p))


def test_committed_split_follows_rule():
    if not os.path.exists(E.SPLIT_FILE):
        pytest.skip('no committed split.json')
    fr = E.load_frozen_split()
    assert fr['planned_configs'] == [dict(c) for c in E.PLANNED_CONFIGS]


def test_hpo_defaults_and_candidates():
    cands = E.hpo_candidates(6, 0)
    assert cands[0] == E.HP_DEFAULT and E.HP_DEFAULT['clip_inputs'] == -1.0
    assert all(c['clip_inputs'] in (-1.0, 0.0) for c in cands)
    plan = dict(epochs=10, patience=5, min_epochs=2, batch_size=32)
    mk, lk, tk = E.hp_kwargs(cands[1], plan)
    assert mk['clip_inputs'] == cands[1]['clip_inputs'] and tk['early_stop'] in ('data', 'mpsp')
    assert E.HPO_METRIC == 'mpsp_mape_cfg'


def test_config_balanced_metrics():
    import numpy as np
    from fermentation_insights.coeff_predictor.evaluate import compute_metrics
    s = pd.DataFrame(dict(sample_id=['x1', 'x2', 'y1'], config_id=['X', 'X', 'Y'], a=1., b=0., c=0., d=0., g=1.))
    p = pd.DataFrame(dict(sample_id=['x1', 'x2', 'y1'], y=0.5, t=50., MPSP=1., recovery=1.))
    pred = dict(coef=np.array([[1.1, 0, 0, 0], [1.1, 0, 0, 0], [1.4, 0, 0, 0]]), g=np.ones(3))
    m = {(k, t): v for k, t, v in compute_metrics(pred, s, p)}
    assert abs(m[('mpsp_mape', 'MPSP')] - 0.2) < 1e-9            # point-weighted (0.1, 0.1, 0.4)
    assert abs(m[('mpsp_mape_cfg', 'MPSP')] - 0.25) < 1e-9       # configuration-balanced (0.1, 0.4)
    assert abs(m[('mpsp_mape_cfg_gmean', 'MPSP')] - 0.2) < 1e-9  # sqrt(0.1 * 0.4)
