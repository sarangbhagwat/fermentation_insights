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


def test_unknown_fields_and_alias_conflicts_are_reported(shipped):
    s, P = shipped
    base = dict(config_id='HP_corn', r=0.3, feedstock_price=0.2)
    ref = _quiet(P.predict, base)[K].to_numpy()
    # a misspelled price is ignored (baseline electricity price) -- with a warning and a suggestion
    out, msgs = _caught(P.predict, dict(base, electricty_price=0.2))
    np.testing.assert_allclose(out[K].to_numpy(), ref)
    assert any("'electricty_price'" in m and "did you mean 'electricity_price'" in m for m in msgs)
    _, msgs = _caught(P.predict, dict(base, foo=1))
    assert any('IGNORED' in m and "'foo'" in m for m in msgs)
    # a price and its alias with different values: the canonical name wins, with a warning
    out, msgs = _caught(P.predict, dict(base, fs_price=0.5))
    np.testing.assert_allclose(out[K].to_numpy(), ref)
    assert any('feedstock_price=0.2 vs. fs_price=0.5' in m for m in msgs)
    _, msgs = _caught(P.predict, dict(base, fs_price=0.2))
    assert not any('alias' in m for m in msgs)
    # training rows (all columns of the processed table) and predict() outputs pass silently
    rows = s.groupby('config_id').head(1).reset_index(drop=True)
    out, msgs = _caught(P.predict, rows)
    assert not any('IGNORED' in m or 'alias' in m for m in msgs), msgs
    _, msgs = _caught(P.predict, pd.concat([rows[['config_id', 'r', 'feedstock_price']], out[K]], axis=1))
    assert not any('IGNORED' in m for m in msgs)


def test_regime_vocabularies_are_validated_and_normalized(shipped):
    _, P = shipped
    ref = _quiet(P.predict, KNOWN)[K].to_numpy()
    for k, v in (('feedstock_type', 'Lignocellulose'), ('feedstock_type', ' lignocellulose '),
                 ('fermentation_mode', 'BATCH'), ('neutralizing_agent', 'ca(oh)2')):
        out, msgs = _caught(P.predict, dict(KNOWN, **{k: v}))
        np.testing.assert_allclose(out[K].to_numpy(), ref, err_msg=f'{k}={v!r}')
        assert out['in_domain'].iloc[0] and not any('not seen in training' in m for m in msgs)
    fed = _quiet(P.predict, dict(KNOWN, fermentation_mode='Fed-batch'))[K].to_numpy()
    np.testing.assert_allclose(fed, _quiet(P.predict, dict(KNOWN, fermentation_mode='fed-batch'))[K].to_numpy())
    for k, v, hint in (('feedstock_type', 'lignocellulosic', 'lignocellulose'), ('fermentation_mode', 'foo', None),
                       ('feedstock_type', 'Lignocellulose-rich', None)):
        with pytest.raises(ValueError, match=k) as e:
            _quiet(P.predict, dict(KNOWN, **{k: v}))
        if hint:
            assert hint in str(e.value)
    # documented / new values not seen in training: unknown level, with a warning naming them
    for k, v in (('fermentation_mode', 'continuous'), ('neutralizing_agent', 'KOH')):
        _, msgs = _caught(P.predict, dict(KNOWN, **{k: v}))
        assert any(f'{k}={v!r} was not seen in training' in m for m in msgs)


@pytest.fixture(scope='module')
def with_ranges(shipped, tmp_path_factory):
    """The release checkpoint with the per-configuration contextual ranges
    (``spec_registry['context_range']``; added to an older checkpoint by
    ``refresh_spec_registry`` from its training samples)."""
    s, P = shipped
    if (P.spec.spec_registry or {}).get('context_range'):
        return P
    import shutil
    from fermentation_insights.coeff_predictor.predict import refresh_spec_registry
    dst = str(tmp_path_factory.mktemp('ck') / 'release')
    shutil.copytree(os.path.join(HERE, 'models', 'release'), dst)
    assert refresh_spec_registry(dst, s)
    return _quiet(CoefficientPredictor.load, dst)


def test_feedstock_capacity_validation_and_reporting(shipped):
    _, P = shipped
    base = dict(config_id='HP_corn', r=0.3, feedstock_price=0.2)
    # a given capacity is used and not reported as defaulted
    out, msgs = _caught(P.predict, dict(base, p_feedstock_capacity=60000.))
    assert 'feedstock_capacity' not in out['contextual_defaults'][0].split(', ')
    assert not any('assumed baseline values: feedstock_capacity' in m for m in msgs)
    ref = _quiet(P.predict, base)
    assert 'feedstock_capacity' in ref['contextual_defaults'][0].split(', ')
    assert not np.allclose(out[K].to_numpy(), ref[K].to_numpy())
    # fs_F_mass (the training tables' name) alone is the capacity, not overridden by the registry
    np.testing.assert_allclose(_quiet(P.predict, dict(base, fs_F_mass=60000.))[K].to_numpy(), out[K].to_numpy())
    # invalid values are rejected (as r <= 0), never silently replaced by the baseline
    for v in (-10, 0, 'n/a', 'abc', float('inf')):
        with pytest.raises(ValueError, match='positive number'):
            _quiet(P.predict, dict(base, p_feedstock_capacity=v))
    with pytest.raises(ValueError, match='positive number'):
        _quiet(P.predict, dict(base, fs_F_mass=-1.))
    # not given: None / NaN / ''
    for v in (None, float('nan'), ''):
        np.testing.assert_allclose(_quiet(P.predict, dict(base, p_feedstock_capacity=v))[K].to_numpy(),
                                   ref[K].to_numpy())


def test_unknown_fields_reach_the_missing_spec_error_and_fs_w_typos(shipped):
    from fermentation_insights.coeff_predictor.predict import MissingSpecError
    _, P = shipped
    with pytest.raises(MissingSpecError, match=r"Unknown field\(s\), ignored: 'confg_id' \(did you mean 'config_id'"):
        _quiet(P.predict, dict(confg_id='HP_corn', r=0.3, feedstock_price=0.2))
    _, msgs = _caught(P.predict, dict(KNOWN, fs_w_watre=0.2))
    assert any('IGNORED' in m and "'fs_w_watre' (did you mean 'fs_w_water'" in m for m in msgs), msgs
    _, msgs = _caught(P.predict, dict(KNOWN, fs_w_water=0.2))
    assert not any('IGNORED' in m for m in msgs), msgs


def test_contextual_parameters_outside_the_configuration_range(with_ranges, shipped):
    s, _ = shipped
    P = with_ranges
    base = dict(config_id='HP_corn', r=0.3, feedstock_price=0.2)
    out, msgs = _caught(P.predict, base)
    assert out['contextual_outside_range'][0] == ''
    assert not any('outside the training range of the configuration' in m for m in msgs)
    cap = float(s.loc[s['config_id'] == 'HP_corn', 'p_feedstock_capacity'].median())
    # 3x the capacity: inside the global range (in_domain unchanged), outside the configuration's
    out, msgs = _caught(P.predict, dict(base, p_feedstock_capacity=3 * cap))
    assert out['contextual_outside_range'][0] == 'feedstock_capacity'
    assert bool(out['in_domain'][0])
    assert any('feedstock_capacity=' in m and 'configuration HP_corn' in m and 'economy of scale' in m
               for m in msgs), msgs
    # the same check for a (pair, feedstock) query and in mpsp()
    spec = dict(feedstock='corn', fermentation_product='HP', final_product='AcrylicAcid', r=0.3,
                feedstock_price=0.2, p_feedstock_capacity=3 * cap)
    assert _quiet(P.predict, spec)['contextual_outside_range'][0] == 'feedstock_capacity'
    m = _quiet(P.mpsp, dict(base, p_feedstock_capacity=3 * cap), y=0.5, t=50.)
    assert m['contextual_outside_range'][0] == 'feedstock_capacity'
    # another contextual parameter, and one never defined for the configuration in training
    assert _quiet(P.predict, dict(base, p_IRR=0.3))['contextual_outside_range'][0] == 'p_IRR'
    out = _quiet(P.predict, dict(config_id='ethanol_cornstover', r=0.3, feedstock_price=0.05,
                                 p_operating_days=200))
    assert out['contextual_outside_range'][0] == 'p_operating_days'
    # training rows are inside their configuration's ranges
    rows = s.groupby('config_id').head(2).reset_index(drop=True)
    assert (_quiet(P.predict, rows)['contextual_outside_range'] == '').all()
