# -*- coding: utf-8 -*-
"""User-API robustness on the shipped checkpoints: boolean spellings (also
pandas-3 'str' columns), the training encoding of a missing neutralizing
agent, a custom feedstock passed with a config_id, by-name queries (domain,
baseline prices), configurations the inputs cannot distinguish, the
uncalibrated-interval warning and the warning texts."""
import json
import os
import warnings

import numpy as np
import pandas as pd
import pytest

from fermentation_insights.coeff_predictor import features as F
from fermentation_insights.coeff_predictor import train as T
from fermentation_insights.coeff_predictor.predict import CoefficientPredictor

K = ['a', 'b', 'c', 'd', 'g']
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_parse_bool_and_to_float_string_dtypes():
    for v, want in ((True, True), ('True', True), (' yes ', True), ('1', True), (1, True), (0., False),
                    ('false', False), ('No', False), (None, None), (float('nan'), None), ('', None)):
        assert F.parse_bool(v, 'x') is want, v
    for bad in ('maybe', 2, 'on'):
        with pytest.raises(ValueError):
            F.parse_bool(bad, 'x')
    for dtype in (object, 'str'):
        s = pd.Series(['True', 'no', None, '1'], dtype=dtype)
        np.testing.assert_array_equal(F._to_float(s).to_numpy(), [1., 0., np.nan, 1.])


def _quiet(f, *a, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return f(*a, **kw)


def _caught(f, *a, **kw):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        out = f(*a, **kw)
    return out, [str(x.message) for x in w]


@pytest.fixture(scope='module')
def shipped():
    proc = os.path.join(HERE, 'data', 'processed', 'samples.csv.gz')
    ck = os.path.join(HERE, 'models', 'release', 'feature_spec.json')
    if not (os.path.exists(proc) and os.path.exists(ck)):
        pytest.skip('processed tables (data/processed/*.csv.gz, tracked in git) or shipped checkpoint '
                    'not available -- restore them from git or rebuild them (README section 8, '
                    '"Data availability")')
    s, _, _ = T.load_tables(os.path.join(HERE, 'data', 'processed'))
    return s, _quiet(CoefficientPredictor.load, os.path.join(HERE, 'models', 'release'))


KNOWN = dict(config_id='HP_neutral_cornstover', feedstock='cornstover', fermentation_product='HP',
             final_product='AcrylicAcid', neutralization=True, neutralizing_agent='Ca(OH)2', r=0.28,
             feedstock_price=0.0517)
SUCC = dict(feedstock='glucose', fermentation_product='SuccinicAcid', final_product='SuccinicAcid',
            neutralizing_agent='Ca(OH)2', r=0.5, feedstock_price=0.3)


def test_boolean_spellings(shipped):
    _, P = shipped
    ref = {v: _quiet(P.predict, dict(SUCC, neutralization=v))[K].to_numpy() for v in (True, False)}
    assert not np.allclose(ref[True], ref[False])
    for v, want in (('True', True), ('yes', True), ('1', True), ('false', False), ('no', False)):
        np.testing.assert_allclose(_quiet(P.predict, dict(SUCC, neutralization=v))[K].to_numpy(), ref[want])
    # a pandas-3 string column
    df = pd.DataFrame([dict(SUCC, neutralization='True')]).astype({'neutralization': 'str'})
    np.testing.assert_allclose(_quiet(P.predict, df)[K].to_numpy(), ref[True])
    with pytest.raises(ValueError, match='neutralization'):
        _quiet(P.predict, dict(SUCC, neutralization='maybe'))


def test_missing_agent_without_neutralization_encodes_as_training(shipped):
    s, _ = shipped
    path = os.path.join(HERE, 'models', 'final_dev')
    split = os.path.join(HERE, 'results', 'split_applied.json')
    if not (os.path.isdir(path) and os.path.exists(split)):
        pytest.skip('final_dev / split not available')
    P = _quiet(CoefficientPredictor.load, path)
    ood = set(json.load(open(split))['ood_test_sample_ids'])
    rows = s[s.sample_id.astype(str).isin(ood)].reset_index(drop=True)
    assert rows['neutralizing_agent'].isna().any()      # e.g. HP_corn: neutralization False, agent NaN
    b, o = P.ensemble.predict_samples(rows, strict=False)
    out = _quiet(P.predict, rows, warn=False)
    np.testing.assert_allclose(out[K].to_numpy(), np.column_stack([o['coef'].mean(0), o['g'].mean(0)]),
                               rtol=1e-6, atol=1e-9)
    ratio, oor = P._ad(b)
    dom = np.isfinite(ratio) & (ratio <= 1) & np.array([not x for x in oor])
    np.testing.assert_array_equal(out['in_domain'].to_numpy(), dom)
    assert not out['outside_range'].str.contains('neutralizing_agent').any()


def _straw(P):
    comp = {d['name']: 0. for d in P.spec.numeric if d['name'].startswith('fs_w_')}
    comp.update(fs_w_glucan=0.30, fs_w_xylan=0.17, fs_w_lignin=0.15, fs_w_ash=0.06, fs_w_water=0.15,
                fs_w_extractives=0.10)
    comp['fs_w_other'] = 1 - sum(v for k, v in comp.items() if k not in ('fs_w_other', F.W_SUGAR))
    comp[F.W_SUGAR] = 0.30 * 180.156 / 162.141 + 0.17 * 150.13 / 132.115
    spec = {**KNOWN, 'feedstock_type': 'lignocellulose', 'feedstock_price': 0.06, **comp}
    spec.pop('feedstock')
    return spec


def test_custom_feedstock_with_config_id_is_a_new_feedstock(shipped):
    _, P = shipped
    p, msgs = _caught(P.predict, _straw(P))
    assert not p['intervals_calibrated'].iloc[0]
    filled = p['filled_from_registry'].iloc[0].split(', ')
    assert 'feedstock' not in filled and not any(f.startswith('fs_') for f in filled)
    assert 'separation_ops' in filled                      # the configuration's downstream fields are kept
    assert any('treated as a NEW feedstock' in m for m in msgs)
    assert any('not in training' in m and 'new feedstock: MPSP' in m for m in msgs)
    # a partial composition is imputed, not taken from corn stover
    part = {k: v for k, v in _straw(P).items() if not k.startswith('fs_w_') or k in
            ('fs_w_glucan', 'fs_w_xylan', 'fs_w_lignin', 'fs_w_water', F.W_SUGAR)}
    q = _quiet(P.predict, part)
    assert 'fs_w_ash' in q['imputed'].iloc[0]
    # the same composition as the configuration's is not a custom feedstock
    reg = P.spec.spec_registry['config']['HP_neutral_cornstover']['values']
    same = {**KNOWN, **{k: v for k, v in reg.items() if k.startswith('fs_w_')}}
    same.pop('feedstock')
    assert _quiet(P.predict, same)['intervals_calibrated'].iloc[0]


def test_by_name_queries_domain_prices_and_twins(shipped):
    s, P = shipped
    b = s[s.sample_kind.astype(str) == 'baseline'].reset_index(drop=True)
    q = _quiet(P.predict, b[['config_id', 'r', 'feedstock_price']])
    full = _quiet(P.predict, b, warn=False)
    np.testing.assert_allclose(q['ad_ratio'].to_numpy(), full['ad_ratio'].to_numpy(), rtol=1e-9, atol=1e-12)
    np.testing.assert_array_equal(q['in_domain'].to_numpy(), full['in_domain'].to_numpy())
    # no prices given: the configuration's baseline electricity / natural-gas prices
    np.testing.assert_allclose(q[K].to_numpy(), full[K].to_numpy(), rtol=1e-6, atol=1e-9)
    tw = P._config_twins()
    assert tw.get('lactic_lowpH_SSCF_cornstover') == ['lactic_lowpH_SHF_cornstover']
    _, msgs = _caught(P.predict, dict(config_id='lactic_neutral_SSCF_cornstover', r=0.3, feedstock_price=0.06))
    assert any('encodes identically' in m and 'SHF' in m for m in msgs)
    _, msgs = _caught(P.predict, dict(config_id='HP_corn', r=0.3, feedstock_price=0.2))
    assert not any('encodes identically' in m for m in msgs)


def test_warning_texts_have_plain_row_numbers(shipped):
    _, P = shipped
    _, msgs = _caught(P.mpsp, KNOWN, y=[0.5, 0.9, 1.5], t=50.)
    assert any('rows [0]' in m for m in msgs) and not any('np.int64' in m for m in msgs)
