# -*- coding: utf-8 -*-
"""The user API (``CoefficientPredictor``) must encode a training
configuration exactly as the validated path (``Ensemble.predict_samples``,
behind every CV / test number): full training rows, by-name queries of
training configurations, heterogeneous batches, and the 'not applicable'
state of contextual parameters."""
import json
import os
import warnings

import numpy as np
import pandas as pd
import pytest

from fermentation_insights.coeff_predictor import train as T
from fermentation_insights.coeff_predictor.features import NOT_APPLICABLE
from fermentation_insights.coeff_predictor.model import ModelConfig
from fermentation_insights.coeff_predictor.predict import CoefficientPredictor, refresh_spec_registry
from fermentation_insights.coeff_predictor.tests.conftest import FAST

K = ['a', 'b', 'c', 'd', 'g']
NAMES = ['config_id', 'feedstock', 'fermentation_product', 'final_product', 'r', 'feedstock_price',
         'electricity_price', 'natural_gas_price']
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture(scope='module')
def serving_tables(tables):
    """Fake tables with the structures that broke the API: a regime
    categorical that is NaN for one feedstock of a pair (HP on corn), a
    feedstock-dependent contextual parameter (cornstover runs 180 days) and
    a parameter not applicable for one configuration only (EtOH on corn)."""
    samples, points, configs = tables
    s = samples.copy()
    hp_corn = (s.fermentation_product == 'HP') & (s.feedstock == 'corn')
    s['neutralizing_agent'] = s['neutralizing_agent'].astype(object)
    s.loc[hp_corn, 'neutralizing_agent'] = np.nan
    s.loc[hp_corn, 'neutralization'] = False
    s.loc[s.feedstock == 'cornstover', 'p_operating_days'] = 180.
    s.loc[(s.fermentation_product == 'EtOH') & (s.feedstock == 'corn'), 'p_IRR'] = np.nan
    # structured downstream fields (FeatureSpec(downstream=True)) that vary by product
    hp = s.fermentation_product == 'HP'
    s['separation_ops'] = np.where(hp, 'centrifugation+evaporation+crystallization', 'centrifugation+distillation')
    s['coproducts'] = np.where(hp, 'electricity', 'none')
    s['upgrading_ops'] = 'none'
    return s, points, configs


@pytest.fixture(scope='module')
def predictor(serving_tables, tmp_path_factory):
    s, points, _ = serving_tables
    ens = T.train_ensemble(s, points, n_members=2, mcfg=ModelConfig(width=16, depth=1, dropout=0.),
                           tcfg=T.TrainConfig(**{**FAST, 'epochs': 6, 'patience': 6}), calib_frac=0., seed=1,
                           spec_kwargs=dict(downstream=True))
    path = str(tmp_path_factory.mktemp('serving') / 'ck')
    ens.save(path)
    return CoefficientPredictor.load(path), path


def _validated(P, rows):
    _, o = P.ensemble.predict_samples(rows, strict=False)
    return np.column_stack([o['coef'].mean(0), o['g'].mean(0)])


def _quiet(f, *a, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return f(*a, **kw)


def test_training_rows_reproduce_validated_path(serving_tables, predictor):
    s, _, _ = serving_tables
    P, _ = predictor
    out = _quiet(P.predict, s, warn=False)
    np.testing.assert_allclose(out[K].to_numpy(), _validated(P, s), rtol=1e-9, atol=1e-12)
    # no training row is flagged for an encoding artefact ('__other__' level)
    assert not out['outside_range'].str.contains('__other__').any()


def test_by_name_query_of_each_configuration_reproduces_its_baseline(serving_tables, predictor):
    s, _, _ = serving_tables
    P, _ = predictor
    b = s[s.sample_kind == 'baseline'].reset_index(drop=True)
    ref = _validated(P, b)
    np.testing.assert_allclose(_quiet(P.predict, b[NAMES])[K].to_numpy(), ref, rtol=1e-9, atol=1e-12)
    # names only (each fake configuration is unique for its pair x feedstock)
    q = _quiet(P.predict, b[NAMES[1:]])
    np.testing.assert_allclose(q[K].to_numpy(), ref, rtol=1e-9, atol=1e-12)
    # config_id alone fills the names too
    q1 = _quiet(P.predict, b[['config_id', 'r', 'feedstock_price', 'electricity_price', 'natural_gas_price']])
    np.testing.assert_allclose(q1[K].to_numpy(), ref, rtol=1e-9, atol=1e-12)


def test_not_applicable_and_feedstock_dependent_contextual_parameters(serving_tables, predictor):
    s, _, _ = serving_tables
    P, _ = predictor
    j_irr = P.spec.columns.index('missing__p_IRR')
    etoh_corn = s[(s.fermentation_product == 'EtOH') & (s.feedstock == 'corn')].iloc[0]
    spec = {c: etoh_corn[c] for c in NAMES[1:]}
    df, det = P._frame(spec, details=True)
    assert pd.isna(df.iloc[0].get('p_IRR')) and 'p_IRR' in det[0]['not_applicable']
    assert P.spec.transform(df, missing_categorical='other').X[0, j_irr] == 1.
    # the same pair on another feedstock: applicable (indicator 0)
    df2, _ = P._frame(dict(spec, feedstock='glucose'))
    assert P.spec.transform(df2, missing_categorical='other').X[0, j_irr] == 0.
    # feedstock-dependent days: the combination's own baseline
    hp_cs = dict(feedstock='cornstover', fermentation_product='HP', final_product='AcrylicAcid', r=0.5,
                 feedstock_price=0.08)
    assert _quiet(P._frame, hp_cs)[0]['p_operating_days'].iloc[0] == 180.
    # a new product on cornstover: the feedstock's baseline of its other products
    new = dict(hp_cs, fermentation_product='Mystery', final_product='Mystery', theoretical_yield=0.6,
               inflection_yield=0.5, final_product_stoich_mass_ratio=1.0)
    df3, det3 = _quiet(P._frame, new, details=True)
    assert df3['p_operating_days'].iloc[0] == 180. and not det3[0]['cross_feedstock']
    # explicit 'not applicable' for a new product: NaN, indicator 1
    df4, det4 = _quiet(P._frame, dict(new, p_IRR=NOT_APPLICABLE), details=True)
    assert pd.isna(df4.iloc[0].get('p_IRR')) and 'p_IRR' in det4[0]['not_applicable']
    assert P.spec.transform(df4, missing_categorical='other', strict=False).X[0, j_irr] == 1.


def test_cross_feedstock_default_warns(serving_tables, predictor):
    P, _ = predictor
    reg = P.spec.spec_registry
    saved = reg['feedstock_context'].pop('glucose')
    try:
        spec = dict(feedstock='glucose', fermentation_product='HP', final_product='AcrylicAcid', r=0.5,
                    feedstock_price=0.3)
        key = 'HP|AcrylicAcid|glucose'
        pfe = reg['pair_feedstock'].pop(key)
        try:
            with pytest.warns(UserWarning, match='OTHER feedstocks'):
                P.predict(spec)
        finally:
            reg['pair_feedstock'][key] = pfe
    finally:
        reg['feedstock_context']['glucose'] = saved


def test_heterogeneous_batch_equals_single_rows(serving_tables, predictor):
    s, _, _ = serving_tables
    P, _ = predictor
    assert P.spec.downstream
    b = s[s.sample_kind == 'baseline']
    known = {c: b[b.config_id == b.config_id.iloc[0]].iloc[0][c] for c in NAMES}
    other = {c: b.iloc[-1][c] for c in NAMES[1:]}
    new = dict(feedstock='glucose', fermentation_product='Mystery', final_product='Mystery', r=0.4,
               feedstock_price=0.4, theoretical_yield=0.9, inflection_yield=0.8,
               final_product_stoich_mass_ratio=1.0, neutralizing_agent='Ca(OH)2', aerobic=True,
               separation_ops='centrifugation', coproducts='electricity')
    rows = [known, new, other]
    batch = _quiet(P.predict, rows)[K].to_numpy()
    single = np.vstack([_quiet(P.predict, r)[K].to_numpy() for r in rows])
    np.testing.assert_allclose(batch, single, rtol=1e-12, atol=0)
    # a NaN multi-valued downstream field is 'not given', like an omitted one
    a = _quiet(P.predict, dict(known, separation_ops=np.nan, coproducts=float('nan')))[K].to_numpy()
    np.testing.assert_allclose(a, single[:1], rtol=1e-12, atol=0)


def test_input_errors(predictor):
    from fermentation_insights.coeff_predictor.features import MissingSpecError
    P, _ = predictor
    with pytest.raises(ValueError, match='no specification rows'):
        P.predict([])
    with pytest.raises(MissingSpecError):
        P.predict({})
    spec = dict(config_id='HP_cornX', feedstock='corn', fermentation_product='HP', final_product='AcrylicAcid',
                r=0.5, feedstock_price=0.2)
    with pytest.warns(UserWarning, match="config_id 'HP_cornX' is not in the registry"):
        P.predict(spec)
    # names that contradict config_id: the configuration's baselines are not used, with a warning
    cid = next(iter(P.spec.spec_registry['config']))
    with pytest.warns(UserWarning, match='baseline values were not used'):
        P.predict(dict(spec, config_id=cid, feedstock='cornstover' if 'corn' in cid and 'stover' not in cid
                       else 'corn'))


def test_refresh_spec_registry(serving_tables, predictor, tmp_path):
    import shutil
    s, _, _ = serving_tables
    _, path = predictor
    dst = str(tmp_path / 'old')
    shutil.copytree(path, dst)
    f = os.path.join(dst, 'feature_spec.json')
    d = json.load(open(f))
    reg = d['spec_registry']
    full = dict(reg)
    for k in ('config', 'pair_feedstock', 'feedstock_context', 'pair_descriptors', 'version'):
        reg.pop(k, None)
    json.dump(d, open(f, 'w'))
    with pytest.warns(UserWarning, match='predates'):
        CoefficientPredictor.load(dst)
    with pytest.raises(ValueError, match='do not reproduce'):
        refresh_spec_registry(dst, s.iloc[:20])
    assert refresh_spec_registry(dst, s)
    assert json.load(open(f))['spec_registry'] == json.loads(json.dumps(full))
    assert not refresh_spec_registry(dst, s)


# ---- the shipped checkpoints (skipped without the processed tables) ----

def _real(name):
    proc = os.path.join(HERE, 'data', 'processed', 'samples.csv.gz')
    ck = os.path.join(HERE, 'models', name, 'feature_spec.json')
    split = os.path.join(HERE, 'results', 'split_applied.json')
    if not (os.path.exists(proc) and os.path.exists(ck) and os.path.exists(split)):
        pytest.skip('processed tables (data/processed/*.csv.gz, tracked in git) or shipped checkpoint '
                    'not available -- restore them from git or rebuild them (README section 8, '
                    '"Data availability")')
    s, _, _ = T.load_tables(os.path.join(HERE, 'data', 'processed'))
    P = _quiet(CoefficientPredictor.load, os.path.join(HERE, 'models', name))
    if name == 'final_dev':
        sp = json.load(open(split))
        test = set(sp['ood_test_sample_ids']) | set(sp['id_test_sample_ids'])
        s = s[~s.sample_id.astype(str).isin(test)]
    return s.reset_index(drop=True), P


@pytest.mark.parametrize('name', ['release', 'final_dev'])
def test_shipped_checkpoint_api_matches_validated_path(name):
    s, P = _real(name)
    if (P.spec.spec_registry or {}).get('version', 1) < 2:
        pytest.fail(f'models/{name}: spec registry predates per-configuration defaults')
    out = _quiet(P.predict, s, warn=False)
    np.testing.assert_allclose(out[K].to_numpy(), _validated(P, s), rtol=1e-6, atol=1e-9)
    b = s[s.sample_kind.astype(str) == 'baseline'].reset_index(drop=True)
    assert b.config_id.nunique() == len(b)
    q = _quiet(P.predict, b[NAMES])
    y, t = b['baseline_yield'].to_numpy(float), b['baseline_titer'].to_numpy(float)
    mp = lambda C: (C[:, 0] + C[:, 1] / y + C[:, 2] / t + C[:, 3] / (y * t)) / C[:, 4]
    rel = mp(q[K].to_numpy()) / mp(_validated(P, b)) - 1
    assert np.abs(rel).max() < 1e-6, pd.Series(rel, index=b.config_id)[np.abs(rel) > 1e-6]
