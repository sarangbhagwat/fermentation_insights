# -*- coding: utf-8 -*-
"""Structured downstream specification: table (downstream.py), the S1 config
patch (data_patches), derived features, the water-removal gate of the model
and the evaluation additions (spec/random split, Gram-norm error, noise floor)."""
import json

import numpy as np
import pandas as pd
import pytest
import torch

from fermentation_insights.coeff_predictor import data_patches as dp
from fermentation_insights.coeff_predictor import dataset, downstream as ds, synthetic
from fermentation_insights.coeff_predictor.evaluate import (compute_metrics, gram_errors, noise_floor,
                                                             spec_determined)
from fermentation_insights.coeff_predictor.experiments import PLANNED_CONFIGS
from fermentation_insights.coeff_predictor.features import (FeatureSpec, DS_PHYS, STRUCTURAL,
                                                            apply_spec_registry, build_spec_registry)
from fermentation_insights.coeff_predictor.model import ModelConfig, make_layout, GATED_CHANNELS
from fermentation_insights.coeff_predictor import train as T
from fermentation_insights.coeff_predictor.tests._fake_tables import make_fake_tables


@pytest.fixture(scope='module')
def broth_tables():
    return make_fake_tables(n_samples_per_config=4, products=['HP', 'EtOH', 'SA', 'BrothAcid'],
                            feedstocks=['glucose', 'cornstover'], seed=3)


def test_table_covers_planned_configs_and_validates():
    planned = {c['config_id'] for c in PLANNED_CONFIGS}
    assert planned == set(ds.TABLE) and len(planned) == 41
    for cid, e in ds.TABLE.items():
        assert set(ds.FIELDS) <= set(e) and e['downstream_spec_version'] == ds.SPEC_VERSION
        assert ds.validate(e)['water_removal_route'] in ds.ROUTES
        assert e['final_product_is_broth'] == (e['water_removal_route'] == 'none')
    ox = ds.spec_for('oxalic_neutral_sugarcane')
    assert ox['final_product_mass_fraction'] is None and ox['neutralization_degree'] == 0.91
    assert 'acidulation' in ds.parse_ops(ds.TABLE['HP_neutral_corn']['separation_ops'])
    assert 'acidulation' not in ds.parse_ops(ds.TABLE['HP_corn']['separation_ops'])
    assert 'LLE' in ds.parse_ops(ds.TABLE['HP_hexanol_glucose']['separation_ops'])
    assert ds.TABLE['ethanol_corn']['coproducts'] == 'DDGS+corn_oil'
    with pytest.raises(ValueError):
        ds.validate(dict(separation_ops='evaporation+teleportation'))
    with pytest.raises(ValueError):
        ds.validate(dict(water_removal_route='sublimate'))
    assert ds.parse_ops('a,b+a') == ['a', 'b'] and ds.parse_ops(None) == [] and ds.parse_ops('none') == []


def _tables(recs):
    cfg = pd.DataFrame([dict(r) for r in recs])
    starts = pd.DataFrame([dict(sample_id=f"{r['biorefinery']}|{r['config_id']}|s1|{i}", harness_version=3,
                                _config=dict(r), _complete=True) for i, r in enumerate(recs)])
    points = pd.DataFrame([dict(biorefinery=r['biorefinery'], config_id=r['config_id'], status='ok', y=0.3, t=50.,
                                sample_id=f"{r['biorefinery']}|{r['config_id']}|s1|{i}", point_id=0)
                           for i, r in enumerate(recs)])
    return cfg, starts, points


def test_S1_patch_versioned_and_idempotent():
    recs = [dict(biorefinery='isobutanol', config_id='IBO_EtOH_corn', adapter_version=2,
                 coproducts='Ethanol,DDGS,crude_oil'),
            dict(biorefinery='oxalic', config_id='oxalic_neutral_sugarcane', adapter_version=2,
                 neutralization_degree=0.91, product_form='clarified_broth'),
            dict(biorefinery='SYN', config_id='SYN_x', adapter_version=2)]
    c1, s1, p1, n1 = dp.apply(*_tables(recs))
    C = c1.set_index('config_id')
    assert C.loc['IBO_EtOH_corn', 'coproducts'] == 'DDGS+corn_oil+electricity+ethanol'
    assert C.loc['IBO_EtOH_corn', 'coproducts_recorded'] == 'Ethanol,DDGS,crude_oil'
    assert C.loc['oxalic_neutral_sugarcane', 'water_removal_route'] == 'none'
    assert 'S1' in C.loc['oxalic_neutral_sugarcane', dp.MARK]
    assert n1['S1']['configs'] == ['isobutanol|IBO_EtOH_corn', 'oxalic|oxalic_neutral_sugarcane']
    assert n1['S1']['configs_without_spec'] == ['SYN|SYN_x']
    S = s1.set_index('sample_id')
    assert S.loc['oxalic|oxalic_neutral_sugarcane|s1|1', '_config']['final_product_is_broth'] is True
    assert 'S1' in S.loc['oxalic|oxalic_neutral_sugarcane|s1|1', dp.MARK]
    # second pass: nothing to do
    c2, s2, p2, n2 = dp.apply(c1, s1, p1)
    assert n2['S1']['config_records'] == 0
    assert c2.drop(columns=[dp.MARK]).equals(c1.drop(columns=[dp.MARK]))
    # an older table version is re-patched
    old = [dict(r, downstream_spec_version=0, water_removal_route='evaporate_water') for r in recs[:2]]
    c3, _, _, n3 = dp.apply(*_tables(old))
    assert n3['S1']['config_records'] == 4      # global records + the samples' own copies
    assert c3.set_index('config_id').loc['oxalic_neutral_sugarcane', 'water_removal_route'] == 'none'


def test_apply_to_frame_and_build(tmp_path):
    df = pd.DataFrame(dict(config_id=['HP_corn', 'HP_corn', 'unknown_cfg'], data_patches=['D3', '', '']))
    out, n = ds.apply_to_frame(df)
    assert n == 2 and out.loc[0, 'water_removal_route'] == 'evaporate_water'
    assert out.loc[0, 'data_patches'] == 'D3;S1' and out.loc[2, 'data_patches'] == ''
    assert pd.isna(out.loc[2, 'water_removal_route'])
    out2, n2 = ds.apply_to_frame(out)
    assert n2 == 0 and out2.equals(out)
    assert set(ds.SAMPLE_FIELDS) <= set(dataset.SPEC_COLUMNS)
    # build(): samples of a known configuration carry the fields
    d = tmp_path / 'raw'
    synthetic.generate(str(d), names=['SYNB'], n_random=1, n_productivity=1)
    pa = next(d.rglob('*.jsonl'))
    recs = [json.loads(line) for line in open(pa)]
    cid0 = next(r['config_id'] for r in recs if r['record'] == 'config')
    out = []
    for r in recs:
        r = dict(r)
        if r.get('config_id') == cid0:
            r['config_id'] = 'TAL_glucose'
        if 'sample_id' in r:
            r['sample_id'] = r['sample_id'].replace(f'|{cid0}|', '|TAL_glucose|')
        out.append(r)
    with open(pa, 'w') as fh:
        fh.writelines(json.dumps(r) + '\n' for r in out)
    data = dataset.build([str(d)])
    S = data.samples[data.samples['config_id'] == 'TAL_glucose']
    assert len(S) and (S['water_removal_route'] == 'crystallize').all()
    assert (S['final_product_mass_fraction'] == 0.896).all()
    assert all('S1' in str(v) for v in S['data_patches'])
    assert data.report['data_patches']['S1']['config_records'] >= 1


def test_features_downstream_switch(broth_tables):
    samples, points, _ = broth_tables
    base = FeatureSpec().fit(samples)
    dsp = FeatureSpec(downstream=True).fit(samples)
    assert not any(c.startswith(('ds_', 'struct_')) or c in DS_PHYS for c in base.columns)
    assert any(c.startswith('separation_principle==') for c in base.columns)
    assert not any(c.startswith('separation_principle==') for c in dsp.columns)
    assert any(c.startswith('water_removal_route==') for c in dsp.columns)
    assert 'ds_sep__crystallization' in dsp.columns and 'phys_base_per_product' in dsp.columns
    ns = len(STRUCTURAL)
    assert dsp.columns[-ns:] == list(STRUCTURAL)
    assert dsp.structural_index(STRUCTURAL[0]) == dsp.n_features - ns
    assert dsp.encoder_index() == list(range(dsp.n_features - ns))
    assert base.encoder_index() == list(range(base.n_features))
    b = dsp.transform(samples)
    gate = b.X[:, dsp.structural_index(STRUCTURAL[0])]
    broth = samples['fermentation_product'].eq('BrothAcid').to_numpy()
    assert (gate[broth] == 0).all() and (gate[~broth] == 1).all()
    assert not dsp.continuous_mask()[-ns:].any()
    # derived values
    df = dsp._prepare(samples)
    r = df[broth].iloc[0]
    assert r['phys_water_removal_a'] == 0 and r['phys_water_removal_c'] == 0
    assert np.isclose(r['phys_base_per_product'], 2 * 0.9 * ds.BASE_EQ_MASS['Ca(OH)2'] / 90.03, rtol=1e-3)
    assert np.isclose(r['phys_gypsum_per_product'], 2 * 0.9 * ds.GYPSUM_EQ_MASS / 90.03, rtol=1e-3)
    h = df[samples['fermentation_product'].eq('HP').to_numpy()].iloc[0]
    assert np.isclose(h['phys_water_removal_a'], -(1 - 0.997) / 0.997) and h['phys_evap_water_c'] == 1000
    e = df[samples['fermentation_product'].eq('EtOH').to_numpy()].iloc[0]
    assert e['phys_distill_c'] == 1000 and e['phys_base_per_product'] == 0 and e['ds_cop__electricity'] == 1
    # serialization round trip; a query without downstream fields: imputed and reported
    re = FeatureSpec.from_json(dsp.to_json())
    np.testing.assert_allclose(re.transform(samples).X, b.X)
    q = samples.drop(columns=list(ds.FIELDS)).iloc[:1]
    bq = dsp.transform(q, strict=False, missing_categorical='other')
    assert STRUCTURAL[0] in bq.imputed_rows[0] and bq.X[0, dsp.structural_index(STRUCTURAL[0])] == 1.
    # old specs (no 'downstream' key) load with it off
    d = base.to_dict()
    d.pop('downstream'), d.pop('structural')
    assert FeatureSpec.from_dict(d).downstream is False


def test_registry_fills_downstream_fields(broth_tables):
    samples = broth_tables[0]
    reg = build_spec_registry(samples)
    filled, f = apply_spec_registry(dict(fermentation_product='EtOH', final_product='EtOH'), reg)
    assert filled['water_removal_route'] == 'distill_product' and 'separation_ops' in f
    filled, f = apply_spec_registry(dict(config_id='oxalic_neutral_sugarcane'), {})
    assert filled['water_removal_route'] == 'none' and filled['final_product_is_broth'] is True


def test_water_gate(broth_tables, tmp_path):
    samples, points, _ = broth_tables
    torch.manual_seed(0)
    spec = FeatureSpec(downstream=True).fit(samples)
    batch = spec.transform(samples)
    cfg = ModelConfig(width=16, depth=2, dropout=0., water_gate=True)
    layout = make_layout(cfg, samples.columns)
    scales = T.fit_target_scales(samples, points, layout, batch, spec.continuous_mask())
    m = T.build_model(spec, layout, scales, cfg, torch.float64).eval()
    X = torch.tensor(batch.X)
    inp = [torch.tensor(np.asarray(a, float)) for a in (batch.log_r, batch.prices, batch.w_sugar, batch.g_max)]
    o = m(X, *inp)
    broth = samples['fermentation_product'].eq('BrothAcid').to_numpy()
    gi = [layout['channels'].index(c) for c in GATED_CHANNELS if c in layout['channels']]
    assert len(gi) == 4
    assert torch.all(o['channels'][broth][:, gi, 2:] == 0) and torch.any(o['channels'][~broth][:, gi, 2:] != 0)
    jng = layout['price_keys'].index('ng')
    assert torch.all(o['intensities'][broth][:, jng, 2:] == 0)
    # the structural input never reaches the encoder: ungated outputs are unchanged
    X2 = X.clone()
    X2[:, spec.structural_index(STRUCTURAL[0])] = 1.
    o2 = m(X2, *inp)
    torch.testing.assert_close(o2['g'], o['g'])
    torch.testing.assert_close(o2['channels'][:, :, :2], o['channels'][:, :, :2])
    # the gate is not part of the saved state
    assert not any(k in m.state_dict() for k in ('gate_ch', 'gate_int', 'enc_index'))
    # end to end: a short ensemble fit, checkpoint round trip and user queries
    ens = T.train_ensemble(samples, points, 1, cfg, T.LossConfig(), T.TrainConfig(epochs=3, batch_size=16,
                                                                                    min_epochs=1, patience=2),
                           calib_frac=0., spec_kwargs=dict(downstream=True))
    _, outs = ens.predict_samples(samples)
    assert np.isfinite(outs['coef']).all()
    import warnings
    from fermentation_insights.coeff_predictor.predict import CoefficientPredictor
    P = CoefficientPredictor(T.Ensemble.load(ens.save(str(tmp_path / 'ck'))))
    assert P.spec.downstream and P.ensemble.members[0].gate_index is not None
    new = dict(feedstock='glucose', fermentation_product='Mystery', final_product='Mystery', r=0.4,
               feedstock_price=0.4, theoretical_yield=0.9, inflection_yield=0.8,
               final_product_stoich_mass_ratio=1.0, fp_MW=100., fp_n_C=3, fp_n_H=4, fp_n_O=4, n_acid_groups=1,
               final_product_is_broth=True, water_removal_route='none', final_product_phase='aqueous_solution',
               separation_ops='centrifugation+acidulation', upgrading_ops='none', coproducts='electricity',
               neutralization_degree=1., acid_regeneration='H2SO4', neutralizing_agent='Ca(OH)2')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        o = P.predict(new)
        o_evap = P.predict(dict(new, water_removal_route='evaporate_water', final_product_is_broth=False,
                                final_product_mass_fraction=0.99))
    names = P.ensemble.channel_names()
    for c in GATED_CHANNELS:
        assert o[f'ch__{c}__c'].iloc[0] == 0 and o[f'ch__{c}__d'].iloc[0] == 0
    assert o_evap['ch__capex_downstream__c'].iloc[0] != 0 and 'capex_downstream' in names
    assert 'struct_water_removal' not in o['imputed'].iloc[0]
    with pytest.raises(ValueError, match='allowed'):
        P.predict(dict(new, separation_ops='evaporation+magic'))
    assert P.known_downstream('ethanol_corn')['water_removal_route'] == 'distill_product'
    assert set(P.downstream_vocabulary['water_removal_route']) == set(ds.ROUTES)
    # an omitted regime field that varies in training is reported unless the checkpoint does not use it
    reg = P.spec.spec_registry['regime']['HP|AcrylicAcid']
    reg['varying'] = ['separation_principle', 'neutralization', 'separation_ops']
    for c in reg['varying']:
        reg['fields'].pop(c, None)
    q = dict(feedstock='glucose', fermentation_product='HP', final_product='AcrylicAcid', r=0.4, feedstock_price=0.3)
    _, det = P._frame(q, details=True)
    assert det[0]['ambiguous'] == ['neutralization', 'separation_ops']


def test_gram_errors_and_split_metrics(broth_tables):
    samples, points, _ = broth_tables
    Y, Tt, M, _, mask = T.pad_points(samples, points)
    C = samples[['a', 'b', 'c', 'd']].to_numpy(float)
    rel, per = gram_errors(C, C, Y, Tt, mask)
    assert np.allclose(rel, 0) and np.allclose(per, 0)
    E = C.copy()
    E[:, 2] += 1.
    rel, per = gram_errors(E, C, Y, Tt, mask)
    # equals the RMS of the numerator-surface error relative to the RMS numerator
    i = 0
    m = mask[i]
    N = C[i, 0] + C[i, 1] / Y[i, m] + C[i, 2] / Tt[i, m] + C[i, 3] / (Y[i, m] * Tt[i, m])
    assert np.isclose(rel[i], np.sqrt(np.mean((1. / Tt[i, m]) ** 2)) / np.sqrt(np.mean(N ** 2)))
    assert np.allclose(per[:, 2], rel) and np.allclose(per[:, [0, 1, 3]], 0)
    pred = dict(coef=E, g=samples['g'].to_numpy(float))
    rows = {(a, b): v for a, b, v in compute_metrics(pred, samples, points)}
    sd = spec_determined(samples)
    assert sd.sum() == (samples['sample_kind'] != 'random').sum() and 0 < sd.sum() < len(sd)
    for k in ('mpsp_mape_spec', 'mpsp_mape_random', 'mpsp_mape_cfg_spec', 'mpsp_mape_cfg_random'):
        assert np.isfinite(rows[(k, 'MPSP')])
    assert np.isclose(rows[('coef_gram_rel', 'abcd')], np.nanmean(rel))
    assert rows[('n_samples_spec', 'rows')] == sd.sum()


def test_noise_floor():
    samples, points, _ = make_fake_tables(n_samples_per_config=14, products=['HP', 'EtOH'],
                                          feedstocks=['glucose'], seed=1)
    df, summ = noise_floor(samples, points, k=3, min_random=8)
    assert len(df) == 2 and summ['n_configs'] == 2
    assert (df['floor_mape'] > 0).all() and np.isfinite(df[['with_bp_mape', 'const_mape']].to_numpy()).all()
