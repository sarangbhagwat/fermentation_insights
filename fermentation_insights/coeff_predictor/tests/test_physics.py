# -*- coding: utf-8 -*-
"""Tests of the physics layer: OLS fit vs. the publication's fit, channel sum
identity, phi/F estimation guards, natural-gas pool detection."""
import numpy as np
import pandas as pd
import pytest

from fermentation_insights.coeff_predictor import physics, synthetic


def _grid(n_y=9, n_t=9, y=(0.1, 0.9), t=(10., 200.)):
    Y, T = np.meshgrid(np.linspace(*y, n_y), np.linspace(*t, n_t))
    return Y.ravel(), T.ravel()


def test_equation_and_basis():
    y, t = _grid(4, 3)
    coef = np.array([0.8, 0.3, 12., 2.])
    z = physics.mpsp_equation(y, t, *coef)
    np.testing.assert_allclose(physics.basis_matrix(y, t) @ coef, z, rtol=1e-14)
    np.testing.assert_allclose(physics.mpsp_equation(y, t, *coef, g=0.5), 2 * z)
    fit = physics.fit_ols(y, t, z)
    np.testing.assert_allclose(fit.coef, coef, rtol=1e-10)
    assert fit.R2 == pytest.approx(1.0) and fit.fit_mape < 1e-12


@pytest.mark.parametrize('noise', [0.0, 0.03])
def test_ols_matches_publication_fit(noise):
    utils = pytest.importorskip('fermentation_insights.utils')
    rng = np.random.default_rng(0)
    y, t = _grid()
    z = physics.mpsp_equation(y, t, 0.9, 0.35, 15., 3.)
    z = z * (1 + noise * rng.standard_normal(len(z))) + noise * 0.2 * (y - 0.5) ** 2
    # emulate infeasible (NaN) points removed from the grid as in the paper
    keep = ~((y < 0.2) & (t > 150))
    a, b, c, d, _ = utils.fit_shifted_rect_hyperbola_two_param(
        (list(y[keep]), list(t[keep])), list(z[keep]), [1, 1, 1, 1])
    fit = physics.fit_ols(y[keep], t[keep], z[keep])
    np.testing.assert_allclose(fit.coef, [a, b, c, d], rtol=1e-6, atol=1e-8)
    # R^2 identical to the publication's get_Rsq definition
    zhat = physics.basis_matrix(y[keep], t[keep]) @ np.array([a, b, c, d])
    rss = np.sum((z[keep] - zhat) ** 2)
    tss = np.sum((z[keep] - z[keep].mean()) ** 2)
    assert fit.R2 == pytest.approx(1 - rss / tss, rel=1e-8)


def test_fit_many_is_linear():
    rng = np.random.default_rng(1)
    y, t = _grid(5, 5)
    Z = rng.normal(size=(len(y), 6))
    C, _ = physics.fit_ols_many(y, t, Z)
    Csum, _ = physics.fit_ols_many(y, t, Z.sum(1, keepdims=True))
    np.testing.assert_allclose(C.sum(1), Csum[:, 0], rtol=1e-12, atol=1e-12)


def _synthetic_sample(name='SYNB', kappa=0., seed=0):
    """Points of one synthetic sample as a canonical-column DataFrame."""
    bio = synthetic.SYNTHETIC_BIOREFINERIES[name]
    ad = synthetic.SyntheticAdapter(name, kappa=kappa)
    config = ad.configurations()[0]
    ad.load(config)
    rng = np.random.default_rng(seed)
    for p in ad.parameters():
        p.setter(p.inverse_cdf(rng.uniform()))
    from fermentation_insights.coeff_predictor.data_generation.harness import (
        collect_area_metrics, collect_tea_metrics, InfeasibleSpecification)
    rows = []
    r = bio['baseline_productivity'] * 1.3
    for y in np.linspace(bio['y_min'], bio['inflection_yield'], 5):
        for t in np.linspace(bio['t_min'], bio['t_max'], 5):
            try:
                ad.set_try(y, t, r)
            except InfeasibleSpecification:
                continue
            ad.simulate()
            res = ad.collect()
            res.update(collect_tea_metrics(ad.tea, ad.product_stream))
            res.update(collect_area_metrics(ad.unit_groups, ad.area_map))
            res.update(y=y, t=t, r=r, status='ok')
            rows.append(res)
    return pd.DataFrame(rows), ad


@pytest.mark.parametrize('name,kappa', [('SYNA', 0.), ('SYNB', 0.), ('SYNC', 0.02)])
def test_channel_sum_identity(name, kappa):
    df, _ = _synthetic_sample(name, kappa)
    d = physics.decompose_sample(df)
    total = np.sum([d.channel_coef[ch] for ch in physics.CHANNELS], axis=0)
    scale = np.max(np.abs(d.coef))
    np.testing.assert_allclose(total, d.coef, rtol=0, atol=1e-11 * scale)
    assert d.sum_check < 1e-11
    # point level: channels sum to MPSP
    M = d.point_contrib[[f'm__{c}' for c in physics.CHANNELS]].sum(axis=1)
    np.testing.assert_allclose(M, df['MPSP'], rtol=1e-12)
    # price-explicit channels = price x intensity
    for ch, q in physics.PRICE_CHANNELS.items():
        p = d.prices[q]
        np.testing.assert_allclose(d.channel_coef[ch], p * d.intensity_coef[q],
                                   rtol=1e-9, atol=1e-12 * scale)
    row = d.to_row()
    for k in physics.COEFF_KEYS:
        s = sum(row[f'ch__{c}__{k}'] for c in physics.CHANNELS)
        assert s == pytest.approx(row[k], abs=1e-11 * scale)


def test_feedstock_channel_analytic():
    # SYNA: constant capacity, constant recovery -> m_fs = p/(w*eta*g*y):
    # primed b'_fs = p/(w*eta), other feedstock coefficients 0
    df, ad = _synthetic_sample('SYNA')
    d = physics.decompose_sample(df)
    bio = ad.bio
    p = ad.values['p_feedstock_price']
    np.testing.assert_allclose(d.channel_coef['feedstock'],
                               [0, p / (bio['w_sugar'] * bio['eta']), 0, 0], atol=1e-10)
    np.testing.assert_allclose(d.intensity_coef['fs'][1], 1 / (bio['w_sugar'] * bio['eta']),
                               rtol=1e-10)
    assert d.phi == pytest.approx(synthetic.capital_charge(ad.values['p_IRR'],
                                                           ad.values['p_income_tax']), rel=1e-9)
    assert d.F_fixed == pytest.approx(synthetic.fixed_cost(bio, ad.values), rel=1e-7)
    assert d.closure_rel < 1e-12


def test_estimate_phi_F_guards():
    I = np.linspace(1e8, 2e8, 10)
    resid = 0.2 * I + 3e6
    phi, F, dg = physics.estimate_phi_F(resid, I)
    assert dg['phi_flag'] == 'ok' and phi == pytest.approx(0.2) and F == pytest.approx(3e6)
    # too few points -> prior / F=0 fallback
    phi, F, dg = physics.estimate_phi_F(resid[:2], I[:2], phi_prior=0.25)
    assert dg['phi_flag'] == 'prior' and phi == 0.25
    assert F == pytest.approx(np.mean(resid[:2] - 0.25 * I[:2]))
    phi, F, dg = physics.estimate_phi_F(resid[:2], I[:2])
    assert dg['phi_flag'] == 'F0' and F == pytest.approx(0, abs=1e-6)
    # constant installed cost -> unidentifiable
    phi, F, dg = physics.estimate_phi_F(np.full(10, 5e7), np.full(10, 1.5e8), phi_prior=0.3)
    assert dg['phi_flag'] == 'prior'
    # negative slope -> clipped to 0 (or prior)
    phi, F, dg = physics.estimate_phi_F(-0.1 * I + 5e7, I)
    assert dg['phi_flag'] == 'clipped' and phi == 0.0
    phi, F, dg = physics.estimate_phi_F(-0.1 * I + 5e7, I, phi_prior=0.2)
    assert dg['phi_flag'] == 'prior' and phi == 0.2


def test_ng_pool_detection_and_override():
    df, _ = _synthetic_sample('SYNB')   # NG is a priced feed (material pool)
    assert physics.decompose_sample(df).diagnostics['ng_pool'] == 'material'
    df2, _ = _synthetic_sample('SYNA')  # NG is a stream utility
    d2 = physics.decompose_sample(df2)
    assert d2.diagnostics['ng_pool'] == 'utility'
    assert d2.diagnostics['materials_other_min_rel'] >= 0
    assert d2.diagnostics['utilities_other_min_rel'] >= 0
    # wrong override only moves cost between the two remainder channels
    d3 = physics.decompose_sample(df2, ng_pool='material')
    np.testing.assert_allclose(d3.coef, d2.coef)
    np.testing.assert_allclose(
        d3.channel_coef['materials_other'] + d3.channel_coef['utilities_other'],
        d2.channel_coef['materials_other'] + d2.channel_coef['utilities_other'], rtol=1e-9,
        atol=1e-12)


def test_missing_cost_columns_fall_into_remainder():
    df, _ = _synthetic_sample('SYNA')
    full = physics.decompose_sample(df)
    d = physics.decompose_sample(df.drop(columns=['natural_gas_flow', 'base_cost']))
    np.testing.assert_allclose(d.coef, full.coef, rtol=1e-12)
    assert d.sum_check < 1e-11
    assert np.allclose(d.channel_coef['natural_gas'], 0)
    assert 'natural_gas_flow' in d.diagnostics['imputed']
    assert d.phi == pytest.approx(full.phi, rel=1e-9)


def test_capital_basis_fci_vs_installed():
    # synthetic FCI = 1.6*(I + 0.175*ISBL): FCI basis is exact, installed-cost
    # basis leaves a (small) closure; without FCI the loader falls back.
    df, ad = _synthetic_sample('SYNB')
    d_fci = physics.decompose_sample(df)
    assert d_fci.diagnostics['capital_basis'] == 'FCI' and d_fci.closure_rel < 1e-12
    d_I = physics.decompose_sample(df, capital_basis='installed_cost')
    assert d_I.closure_rel > 1e-6 and d_I.sum_check < 1e-11
    np.testing.assert_allclose(d_I.coef, d_fci.coef, rtol=1e-12)
    d_fb = physics.decompose_sample(df.drop(columns=['FCI']))
    assert d_fb.diagnostics['capital_basis'] == 'installed_cost'
    assert d_fb.phi == pytest.approx(d_I.phi, rel=1e-12)
    # capital channels sum to phi*FCI at every point
    cap = sum(d_fci.point_contrib[f'm__{c}'] for c in physics.CAPEX_CHANNEL_AREAS)
    h, Ff = df['operating_hours'], df['F_final_product']
    np.testing.assert_allclose(cap * Ff * h, d_fci.phi * df['FCI'], rtol=1e-10)
