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


def test_capacity_response_check(with_ranges, shipped):
    from fermentation_insights.coeff_predictor import report_checks as RC
    s, _ = shipped
    b = s[(s['config_id'] == 'HP_corn') & (s['sample_kind'] == 'baseline')].iloc[0]
    cap = float(b['p_feedstock_capacity'])
    sim = pd.DataFrame(dict(config_id='HP_corn', mult=[1., 2.], p_feedstock_capacity=[cap, 2 * cap],
                            y=float(b['baseline_yield']), t=float(b['baseline_titer']), r=float(b['r']),
                            status='ok', MPSP=[1.40, 1.33]))
    t = RC.capacity_response(with_ranges, sim, s).set_index('mult')
    assert t.loc[1., 'dMPSP_pred_pct'] == 0. and t.loc[2., 'dMPSP_sim_pct'] == pytest.approx(-5., abs=1e-9)
    assert t.loc[2., 'contextual_outside_range'] == 'feedstock_capacity'
    assert t.loc[1., 'contextual_outside_range'] == ''
    assert np.isfinite(t.loc[2., 'share_of_simulated_change'])


def test_numeric_fields_are_validated(shipped):
    """A given numeric field that is not a number raises (it would otherwise
    be read as missing and encoded as 'not applicable' / imputed)."""
    _, P = shipped
    base = dict(config_id='HP_corn', r=0.3, feedstock_price=0.2)
    for field, v in (('p_IRR', '10%'), ('p_operating_days', '300 days'), ('p_operating_days', 'abc'),
                     ('electricity_price', 'abc'), ('T_fermentation', 'warm'), ('r', 'fast')):
        with pytest.raises(ValueError, match=f'not a number: {field}='):
            _quiet(P.predict, dict(base, **{field: v}))
    # numbers given as strings, omitted values and 'n/a' (not applicable) are accepted
    np.testing.assert_allclose(_quiet(P.predict, dict(base, p_IRR='0.1'))[K].to_numpy(),
                               _quiet(P.predict, dict(base, p_IRR=0.1))[K].to_numpy())
    for v in (None, float('nan'), ''):
        np.testing.assert_allclose(_quiet(P.predict, dict(base, p_IRR=v))[K].to_numpy(),
                                   _quiet(P.predict, base)[K].to_numpy())
    assert 'p_IRR' in _quiet(P.predict, dict(base, p_IRR='n/a'))['not_applicable'][0].split(', ')


def test_unknown_composition_key_is_ignored(shipped):
    """An fs_w_* key the checkpoint does not know is ignored (and reported),
    so it does not turn a by-name configuration into a custom feedstock."""
    _, P = shipped
    base = dict(config_id='HP_corn', r=0.3, feedstock_price=0.2)
    out, msgs = _caught(P.predict, dict(base, fs_w_watre=0.2))
    np.testing.assert_allclose(out[K].to_numpy(), _quiet(P.predict, base)[K].to_numpy())
    assert any('IGNORED' in m and 'fs_w_watre' in m for m in msgs), msgs
    assert not any('custom' in m.lower() and 'feedstock' in m.lower() for m in msgs), msgs


def test_dataframe_index_is_kept(shipped):
    _, P = shipped
    df = pd.DataFrame([dict(config_id='HP_corn', r=0.3, feedstock_price=0.2),
                       dict(config_id='TAL_glucose', r=1., feedstock_price=0.3),
                       dict(config_id='HP_corn', r=0.5, feedstock_price=0.2)], index=[10, 3, 7])
    out = _quiet(P.predict, df)
    assert out.index.tolist() == [10, 3, 7]
    ref = _quiet(P.predict, df.reset_index(drop=True))
    np.testing.assert_allclose(out[K].to_numpy(), ref[K].to_numpy())
    # mpsp() numbers the spec rows positionally
    assert sorted(_quiet(P.mpsp, df, y=0.5, t=50.)['row'].unique()) == [0, 1, 2]


def test_recovery_series_is_aligned_with_the_spec_index(shipped):
    _, P = shipped
    df = pd.DataFrame([dict(config_id='HP_corn', r=0.3, feedstock_price=0.2),
                       dict(config_id='TAL_glucose', r=1., feedstock_price=0.3),
                       dict(config_id='HP_corn', r=0.5, feedstock_price=0.2)], index=['p', 'q', 'r'])
    rec = pd.Series([0.5, 0.6, 0.7], index=['r', 'q', 'p'])
    out = _quiet(P.predict, df, recovery=rec)
    assert out['g'].tolist() == [0.7, 0.6, 0.5]
    np.testing.assert_allclose(_quiet(P.mpsp, df, y=0.5, t=50., recovery=rec)['MPSP'].to_numpy(),
                               _quiet(P.mpsp, df, y=0.5, t=50., recovery=[0.7, 0.6, 0.5])['MPSP'].to_numpy())
    # positional sequences, a Series with the spec's index
    assert _quiet(P.predict, df, recovery=rec.to_numpy())['g'].tolist() == [0.5, 0.6, 0.7]
    assert _quiet(P.predict, df, recovery=rec.sort_index())['g'].tolist() == [0.7, 0.6, 0.5]
    with pytest.raises(ValueError, match='without values for the spec row label'):
        _quiet(P.predict, df, recovery=pd.Series([0.5, 0.6], index=['p', 'q']))


def test_tea_settings_fixed_by_the_configuration(shipped):
    """Contextual TEA parameters a configuration's adapter does not sample
    (NaN in its training rows, registry/tea_settings.json): reported in
    ``tea_fixed``; the configuration's own value is encoded as missing, as in
    training; any other value warns that it is applied relative to the
    imputed value."""
    from fermentation_insights.coeff_predictor import predict as PR
    if not os.path.exists(PR.TEA_SETTINGS_PATH):
        pytest.skip('registry/tea_settings.json not available')
    _, P = shipped
    tab = json.load(open(PR.TEA_SETTINGS_PATH))['configs']
    assert set(tab) == set(P.spec.spec_registry['config'])
    assert tab['ethanol_corn']['fixed']['p_IRR'] == 0.15 and tab['IBO_EtOH_corn']['fixed']['p_IRR'] == 0.15
    assert abs(tab['ethanol_cornstover']['fixed']['p_operating_days'] - 350.4) < 1e-9
    assert not tab['HP_corn']['fixed']
    base = dict(config_id='ethanol_corn', r=1., feedstock_price=0.15)
    ref = _quiet(P.predict, base)
    assert ref['tea_fixed'][0] == 'p_IRR=0.15, p_operating_days=330'
    # by names that resolve to the configuration: the same
    byname = dict(feedstock='corn', fermentation_product='Ethanol', final_product='Ethanol', r=1.,
                  feedstock_price=0.15)
    assert _quiet(P.predict, byname)['tea_fixed'][0] == 'p_IRR=0.15, p_operating_days=330'
    assert _quiet(P.predict, dict(config_id='ethanol_cornstover', r=1., feedstock_price=0.08))['tea_fixed'][0] \
        == 'p_IRR=0.1, p_operating_days=350.4'
    out, msgs = _caught(P.predict, dict(base, p_IRR=0.15))
    np.testing.assert_allclose(out[K].to_numpy(), ref[K].to_numpy())
    assert out['contextual_outside_range'][0] == '' and 'p_IRR' in out['not_applicable'][0]
    assert any('fixed TEA value of ethanol_corn' in m for m in msgs), msgs
    out, msgs = _caught(P.predict, dict(base, p_IRR=0.12))
    assert any('fixes p_IRR at 0.15' in m and 'imputes 0.1' in m for m in msgs), msgs
    assert 'p_IRR' in out['contextual_outside_range'][0]
    # a configuration that samples the parameter: nothing fixed
    out, msgs = _caught(P.predict, dict(config_id='HP_corn', r=0.3, feedstock_price=0.2, p_IRR=0.11))
    assert out['tea_fixed'][0] == '' and not any('TEA' in m for m in msgs), msgs


@pytest.fixture(scope='module')
def with_family_coverage(shipped, tmp_path_factory):
    """The release checkpoint with known per-family coverage metadata
    (``meta['coverage_by_family']``, as the experiment pipeline writes it)."""
    import shutil
    dst = str(tmp_path_factory.mktemp('ck_cov') / 'release')
    shutil.copytree(os.path.join(HERE, 'models', 'release'), dst)
    mp = os.path.join(dst, 'meta.json')
    m = json.load(open(mp))
    m.setdefault('meta', {})['coverage_by_family'] = dict(
        nominal=0.9,
        known_configuration={'TAL': dict(a=0.95, b=1.0, c=0.76, d=1.0, g=0.44, mpsp=0.91), 'HP': dict(g=1.0)},
        new_configuration={'Ethanol': dict(b=0.67, g=1.0, mpsp=0.77), 'TAL': dict(g=0.43, mpsp=0.91)})
    json.dump(m, open(mp, 'w'))
    return _quiet(CoefficientPredictor.load, dst)


def test_family_undercoverage_is_flagged(with_family_coverage):
    P = with_family_coverage
    specs = [dict(config_id='TAL_glucose', r=1., feedstock_price=0.3), dict(config_id='HP_corn', r=0.3,
                                                                              feedstock_price=0.2)]
    out, msgs = _caught(P.predict, specs)
    assert out['intervals_undercovered'].tolist() == ['c 0.76, g 0.44', '']
    assert out['intervals_calibrated'].all()
    assert any('UNDER-cover' in m and 'TAL' in m and 'g 0.44' in m for m in msgs), msgs
    # a user recovery: its g interval is the given value, not flagged
    assert _quiet(P.predict, specs[0], recovery=0.7)['intervals_undercovered'][0] == 'c 0.76'
    # mpsp(): the MPSP key only; a new ethanol configuration (ethanol from
    # glucose: not in training) uses the LOCO coverage
    eth = dict(feedstock='glucose', fermentation_product='Ethanol', final_product='Ethanol', r=1.,
               feedstock_price=0.3)
    m, msgs = _caught(P.mpsp, [specs[0], eth], y=0.4, t=60.)
    assert m.groupby('row')['intervals_undercovered'].first().tolist() == ['', 'mpsp 0.77']
    assert any('new configuration' in x and 'mpsp 0.77' in x for x in msgs), msgs
    assert _quiet(P.predict, eth)['intervals_undercovered'][0] == 'b 0.67'
    # names that resolve to one training configuration encode as it: the
    # known-configuration flags, as by config_id
    eth_corn = dict(eth, feedstock='corn', feedstock_price=0.2)
    for spec in (eth_corn, dict(eth_corn, config_id='ethanol_corn')):
        assert _quiet(P.predict, spec)['intervals_undercovered'][0] == ''
        assert _quiet(P.mpsp, spec, y=0.4, t=60.)['intervals_undercovered'][0] == ''
    tal = dict(feedstock='glucose', fermentation_product='TAL', final_product='TAL', r=1., feedstock_price=0.3)
    assert _quiet(P.predict, tal)['intervals_undercovered'][0] == 'c 0.76, g 0.44'
    # a config_id whose regime differs from the configuration's is a new
    # configuration; a contextual parameter keeps it known
    assert _quiet(P.predict, dict(specs[0], T_fermentation=310.))['intervals_undercovered'][0] == 'g 0.43'
    assert _quiet(P.predict, dict(specs[0], p_IRR=0.11))['intervals_undercovered'][0] == 'c 0.76, g 0.44'


def test_shipped_checkpoints_carry_family_coverage(shipped):
    """The pipeline records the per-family coverage of the shipped quantiles
    (development out-of-fold scores) in both checkpoints; the known TAL
    recovery under-coverage (REPORT.md section 6) is flagged by name."""
    _, P = shipped
    cbf = (P.ensemble.meta or {}).get('coverage_by_family')
    if not cbf:
        pytest.skip('checkpoint without coverage_by_family (written by experiments.py)')
    assert cbf['known_configuration']['TAL']['g'] < 0.6
    assert cbf['new_configuration']['LacticAcid']['mpsp'] < 0.8
    out = _quiet(P.predict, dict(config_id='TAL_glucose', r=1., feedstock_price=0.3))
    assert 'g 0.' in out['intervals_undercovered'][0]
    fd = json.load(open(os.path.join(HERE, 'models', 'final_dev', 'meta.json')))
    assert fd['meta']['coverage_by_family']['known_configuration']['TAL']['g'] < 0.6


@pytest.fixture(scope='module')
def with_y_domain(shipped, tmp_path_factory):
    s, P = shipped
    if (P.spec.spec_registry or {}).get('y_domain'):
        return P
    import shutil
    from fermentation_insights.coeff_predictor.predict import refresh_spec_registry
    dst = str(tmp_path_factory.mktemp('ck_y') / 'release')
    shutil.copytree(os.path.join(HERE, 'models', 'release'), dst)
    assert refresh_spec_registry(dst, s)
    return _quiet(CoefficientPredictor.load, dst)


def test_yields_below_the_training_range_warn(with_y_domain, shipped):
    s, _ = shipped
    P = with_y_domain
    ymin = float(s.loc[s['config_id'] == 'succinic_glucose', 'y_min'].min())
    q = dict(config_id='succinic_glucose', r=0.5, feedstock_price=0.3)
    _, msgs = _caught(P.mpsp, q, y=0.2 * ymin, t=60.)
    assert any('y below the training yield range' in m for m in msgs), msgs
    _, msgs = _caught(P.mpsp, q, y=2 * ymin, t=60.)
    assert not any('y below the training yield range' in m for m in msgs), msgs


def test_training_rows_encode_as_their_configuration(with_family_coverage, shipped):
    """Training rows vary in the sample's own inflection yield and (lactic,
    isobutanol) the feedstock composition; they still encode as their
    configuration (the known-configuration coverage is measured on such
    rows), while an inflection yield outside the configuration's training
    spread makes a new configuration."""
    s, _ = shipped
    P = with_family_coverage
    rows = s[s['sample_kind'] == 'random'].groupby('config_id').head(3).reset_index(drop=True)
    _, _, det = _quiet(P._batch, rows, None, False)
    assert all(c in d['training_configs'] for c, d in zip(rows['config_id'], det))
    tal = s[(s['config_id'] == 'TAL_glucose') & (s['sample_kind'] == 'random')].head(2).reset_index(drop=True)
    out = _quiet(P.predict, tal)
    assert out['intervals_undercovered'].tolist() == ['c 0.76, g 0.44'] * 2
    far = tal.assign(inflection_yield=0.5 * tal['inflection_yield'])
    assert _quiet(P.predict, far)['intervals_undercovered'].tolist() == ['g 0.43'] * 2


def test_cell_mass_yield_sets_the_inflection_yield(shipped):
    """A given p_cell_mass_yield of a known configuration sets its inflection
    yield by the configuration's training relation (HP: exactly
    (1 - p_cell_mass_yield) x theoretical yield), unless inflection_yield is
    given; mpsp()'s fitted-domain check uses the derived value."""
    _, P = shipped
    q = dict(config_id='HP_cornstover', r=0.28, feedstock_price=0.06)
    df, _, det = _quiet(P._batch, dict(q, p_cell_mass_yield=0.2), None, False)
    assert df['inflection_yield'][0] == pytest.approx(0.8 * df['theoretical_yield'][0], rel=1e-3)
    assert det[0]['inflection_derived'][2] == 'HP_cornstover'
    out, msgs = _caught(P.predict, dict(q, p_cell_mass_yield=0.2))
    assert 'inflection_yield (from p_cell_mass_yield)' in out['filled_from_registry'][0].split(', ')
    assert any('derived from the given p_cell_mass_yield' in m for m in msgs), msgs
    # given inflection yield: kept; no cell-mass yield: registry baseline
    df2, _, _ = _quiet(P._batch, dict(q, p_cell_mass_yield=0.2, inflection_yield=0.7), None, False)
    assert df2['inflection_yield'][0] == 0.7
    df3, _, det3 = _quiet(P._batch, q, None, False)
    assert 'inflection_derived' not in det3[0]
    # TAL_glucose at the high end of its cell-mass range: the derived
    # inflection yield is below the registry one, and mpsp() warns above it
    tal = dict(config_id='TAL_glucose', r=0.2, feedstock_price=0.4)
    infl0 = float(_quiet(P._batch, tal, None, False)[0]['inflection_yield'][0])
    dft = _quiet(P._batch, dict(tal, p_cell_mass_yield=0.405), None, False)[0]
    infl = float(dft['inflection_yield'][0])
    assert infl < infl0
    y = 0.5 * (infl + infl0)
    _, msgs = _caught(P.mpsp, dict(tal, p_cell_mass_yield=0.405), y=y, t=30.)
    assert any('y above the inflection yield' in m for m in msgs), msgs
    _, msgs = _caught(P.mpsp, tal, y=y, t=30.)
    assert not any('y above the inflection yield' in m for m in msgs), msgs


def test_derived_inflection_yield_stays_in_the_training_spread(shipped):
    """At the ends of a configuration's cell-mass-yield range the derived
    inflection yield stays inside the training spread, so the row still counts
    as the known configuration; a names-only row that encodes no training
    configuration warns that the registry inflection yield was kept."""
    _, P = shipped
    for cid in ('TAL_glucose', 'TAL_cornstover', 'TAL_SA_cornstover', 'HP_corn'):
        lo, hi = P._inflection_ranges(cid)[:2]
        rows = pd.DataFrame([dict(config_id=cid, r=0.3, feedstock_price=0.2, p_cell_mass_yield=v)
                             for v in np.linspace(lo, hi, 5)])
        _, _, det = _quiet(P._batch, rows, None, False)
        assert all(cid in d['training_configs'] for d in det), cid
        assert all(d.get('inflection_derived') for d in det), cid
    q = dict(feedstock='corn', fermentation_product='HP', final_product='AcrylicAcid', neutralization=False,
             r=0.3, feedstock_price=0.2, p_cell_mass_yield=0.22)
    _, msgs = _caught(P.predict, q)
    assert any('no training relation' in m for m in msgs), msgs
    _, msgs = _caught(P.predict, dict(q, config_id='HP_corn'))
    assert not any('no training relation' in m for m in msgs), msgs
    assert any('derived from the given p_cell_mass_yield' in m for m in msgs), msgs
    # a cell-mass yield whose derived inflection yield is the registry one: no warning
    infl0 = float(_quiet(P._batch, dict(config_id='HP_corn', r=0.3, feedstock_price=0.2), None, False)[0]
                  ['inflection_yield'][0])
    slope, icpt, _ = P._inflection_relation('HP_corn')
    _, msgs = _caught(P.predict, dict(config_id='HP_corn', r=0.3, feedstock_price=0.2,
                                      p_cell_mass_yield=(infl0 - icpt) / slope))
    assert not any('derived from the given p_cell_mass_yield' in m for m in msgs), msgs


def test_given_contextual_parameters_with_unreliable_response_warn(shipped):
    _, P = shipped
    q = dict(config_id='lactic_lowpH_SHF_cornstover', r=0.9, feedstock_price=0.06)
    for extra, want in ((dict(p_acid_price=0.1), True), (dict(p_turbogenerator_efficiency=0.85), True),
                        (dict(p_IRR=0.11), False), (dict(p_operating_days=330.), False),
                        (dict(p_acid_price='n/a'), False), ({}, False)):
        for f, kw in ((P.predict, {}), (P.mpsp, dict(y=0.5, t=80.))):
            _, msgs = _caught(f, dict(q, **extra), **kw)
            got = any('not reliable' in m and 'REPORT.md section 8 (ii-d)' in m for m in msgs)
            assert got == want, (extra, msgs)
    from fermentation_insights.coeff_predictor.predict import CONTEXT_RESPONSE_UNRELIABLE
    assert set(CONTEXT_RESPONSE_UNRELIABLE).isdisjoint({'p_IRR', 'p_operating_days'})
