# -*- coding: utf-8 -*-
"""FeatureSpec tests: serialization, single-dict transform, leakage guards."""
import json

import numpy as np
import pandas as pd
import pytest

from fermentation_insights.coeff_predictor.features import (
    FeatureSpec, MissingSpecError, OTHER, PRICE_COLUMNS, psat_water, register_chemical,
    resolve_chemical_descriptors, ApplicabilityDomain, _USER_CHEMICALS)


def test_json_round_trip_and_single_dict(fitted):
    spec, samples = fitted['spec'], fitted['samples']
    s = spec.to_json()
    json.loads(s)  # valid JSON (no NaN tokens)
    assert 'NaN' not in s
    spec2 = FeatureSpec.from_json(s)
    assert spec2.columns == spec.columns
    b1 = spec.transform(samples)
    b2 = spec2.transform(samples)
    for k in ('X', 'log_r', 'prices', 'w_sugar', 'g_max'):
        assert np.array_equal(getattr(b1, k), getattr(b2, k))
    rec = samples.iloc[7].to_dict()
    b3 = spec2.transform(rec)
    assert b3.X.shape == (1, spec.n_features)
    assert np.allclose(b3.X[0], b1.X[7]) and np.isclose(b3.log_r[0], b1.log_r[7])


def test_json_file_round_trip(fitted, tmp_path):
    p = tmp_path / 'fs.json'
    fitted['spec'].to_json(str(p))
    assert FeatureSpec.from_json(str(p)).columns == fitted['spec'].columns


def test_never_uses_ids_bp_or_targets(fitted):
    cols = fitted['spec'].columns
    banned = ('sample_id', 'config_id', 'biorefinery', 'seed', 'sample_kind', 'a', 'b', 'c', 'd', 'g',
              'fit_mape', 'R2', 'phi', 'F_fixed', 'y_min', 'y_max', 't_min', 't_max', 'r_mult', 'fs_price')
    for c in cols:
        base = c.split('==')[0].replace('missing__', '')
        assert base not in banned and not base.startswith(('bp_', 'ch__', 'int__'))
    # price-explicit spec keeps prices out of the encoder ...
    for c in PRICE_COLUMNS + ('p_feedstock_price', 'p_electricity_price', 'p_natural_gas_price'):
        assert c not in cols
    # ... the non-explicit spec includes them
    spec2 = FeatureSpec(price_explicit=False).fit(fitted['samples'])
    assert all(c in spec2.columns for c in PRICE_COLUMNS)


def test_unseen_category_maps_to_other(fitted):
    spec, samples = fitted['spec'], fitted['samples']
    rec = samples.iloc[0].to_dict()
    rec['feedstock_type'] = 'microalgae'
    b = spec.transform(rec)
    j = spec.columns.index(f'feedstock_type=={OTHER}')
    assert b.X[0, j] == 1.
    others = [i for i, c in enumerate(spec.columns) if c.startswith('feedstock_type==') and i != j]
    assert np.all(b.X[0, others] == 0)


def test_missing_required_raises(fitted):
    rec = fitted['samples'].iloc[0].to_dict()
    for k in ('r', 'theoretical_yield'):
        bad = dict(rec)
        bad.pop(k)
        with pytest.raises(MissingSpecError, match=k):
            fitted['spec'].transform(bad)
    # feedstock price may come from an alias
    alias = dict(rec)
    alias.pop('feedstock_price')
    b = fitted['spec'].transform(alias)
    assert np.isclose(b.prices[0, 0], rec['p_feedstock_price'])


def test_parameters_baseline_imputed_with_indicator(fitted):
    spec, samples = fitted['spec'], fitted['samples']
    assert 'p_boiler_efficiency' in spec.indicators  # missing for the EtOH adapter
    rec = samples.iloc[0].to_dict()
    rec['p_boiler_efficiency'] = np.nan
    b = spec.transform(rec)
    d = next(x for x in spec.numeric if x['name'] == 'p_boiler_efficiency')
    j = spec.columns.index('p_boiler_efficiency')
    base = samples.loc[samples.sample_kind == 'baseline', 'p_boiler_efficiency'].median()
    assert np.isclose(d['fill'], base)
    assert np.isclose(b.X[0, j], (base - d['mean']) / d['std'])
    assert b.X[0, spec.columns.index('missing__p_boiler_efficiency')] == 1.


def test_registry_descriptors_and_hook(fitted):
    samples = fitted['samples']
    df = resolve_chemical_descriptors(samples.iloc[:1].assign(fp_dor=99.))  # HP: registry wins
    assert np.isclose(df['fp_dor'].iloc[0], 4.0)
    assert np.isfinite(df['fp_log10_alpha_inf'].iloc[0])
    # unknown chemical falls back to its own fp_* columns (dor from atoms)
    new = samples[samples.ferm_product_chemical == 'NewAcid'].iloc[:1]
    d2 = resolve_chemical_descriptors(new)
    assert np.isclose(d2['fp_dor'].iloc[0], (4 * 4 + 8 - 2 * 3) / 4)
    # user hook overrides
    d3 = resolve_chemical_descriptors(new, descriptors={'NewAcid': dict(MW=104.1, n_C=4, n_H=8, n_O=3,
                                                                        n_N=0, gamma_inf=3.0)})
    assert np.isclose(d3['fp_gamma_inf'].iloc[0], 3.0)
    register_chemical('NewAcid', dict(MW=104.1, n_C=4, n_H=8, n_O=3, n_N=0, Tb=480., gamma_inf=2.0))
    try:
        assert np.isclose(resolve_chemical_descriptors(new)['fp_gamma_inf'].iloc[0], 2.0)
    finally:
        _USER_CHEMICALS.pop('NewAcid', None)


def test_derived_physics_features(fitted):
    spec, samples = fitted['spec'], fitted['samples']
    df = spec._prepare(samples.iloc[:1])
    assert np.isclose(df['phys_fp_Tb_minus_water'].iloc[0], df['fp_Tb'].iloc[0] - 373.15)
    assert np.isclose(df['phys_inflection_over_theoretical'].iloc[0], 0.8)
    assert np.isclose(df['phys_ty_over_Ymax_stoich'].iloc[0],
                      samples['theoretical_yield'].iloc[0] / df['fp_Ymax_stoich'].iloc[0])
    assert 3000 < psat_water(310.15) < 9000  # ~6.3 kPa


def test_applicability_domain(fitted):
    spec, batch = fitted['spec'], fitted['batch']
    Z = spec.full_matrix(batch)
    ad = ApplicabilityDomain().fit(Z, spec.full_columns())
    r_in, o_in = ad.score(Z[:5])
    assert np.all(r_in <= 1.0 + 1e-9) and not any(o_in)
    far = Z[:1] + 25.
    r_out, o_out = ad.score(far)
    assert r_out[0] > 1 and len(o_out[0]) > 0
    ad2 = ApplicabilityDomain.from_dict(json.loads(json.dumps(ad.to_dict())), ad.Z)
    assert np.allclose(ad2.score(far)[0], r_out)


def test_nonsugar_solids_not_computed_as_dry(fitted):
    spec, samples = fitted['spec'], fitted['samples']
    rec = samples.iloc[0].to_dict()
    rec.pop('fs_w_water')
    b = spec.transform([rec])
    assert 'fs_w_water' in b.imputed and 'phys_fs_nonsugar_solids' in b.imputed
    assert b.imputed_rows == [b.imputed] and not b.invalid.any()


def test_readme_field_reference_matches_field_units():
    """README.md section 5 is generated from FIELD_UNITS: every field has its
    row with the current unit text (regenerate the row if this fails)."""
    import os
    from fermentation_insights.coeff_predictor import features
    path = os.path.join(os.path.dirname(features.__file__), 'README.md')
    with open(path) as fh:
        text = fh.read()
    stale = [f'| `{k}` | {u.replace("|", chr(92) + "|")} |' for k, u in features.FIELD_UNITS.items()
             if f'| `{k}` | {u.replace("|", chr(92) + "|")} |' not in text]
    assert not stale, 'README field reference out of date; expected rows starting with: ' + '; '.join(stale)
