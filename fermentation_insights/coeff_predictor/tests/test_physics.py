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


def _synthetic_sample(name='SYNB', kappa=0., seed=0, v2=False):
    """Points of one synthetic sample as a canonical-column DataFrame
    (``v2``: with the harness v2 columns: NG per pool, labor, OSBL installed
    cost per area, per-point fs_*)."""
    bio = synthetic.SYNTHETIC_BIOREFINERIES[name]
    ad = synthetic.SyntheticAdapter(name, kappa=kappa)
    config = ad.configurations()[0]
    ad.load(config)
    rng = np.random.default_rng(seed)
    for p in ad.parameters():
        p.setter(p.inverse_cdf(rng.uniform()))
    from fermentation_insights.coeff_predictor.data_generation.harness import (
        collect_area_metrics, collect_tea_metrics, collect_natural_gas, collect_labor,
        describe_feedstock, InfeasibleSpecification)
    rows = []
    r = bio['baseline_productivity'] * 1.3
    with synthetic.biosteam_stub():
        for y in np.linspace(bio['y_min'], bio['inflection_yield'], 5):
            for t in np.linspace(bio['t_min'], bio['t_max'], 5):
                try:
                    ad.set_try(y, t, r)
                except InfeasibleSpecification:
                    continue
                ad.simulate()
                res = ad.collect()
                res.update(collect_tea_metrics(ad.tea, ad.product_stream))
                if v2:
                    res.update(collect_area_metrics(ad.unit_groups, ad.area_map,
                                                    set(getattr(ad.tea, 'OSBL_units', None) or ())))
                    res.update(collect_natural_gas(ad.system))
                    res.update(collect_labor(ad.tea))
                    res.update(describe_feedstock(ad.feedstock_stream))
                else:
                    res.update(collect_area_metrics(ad.unit_groups, ad.area_map))
                    res = {k: v for k, v in res.items() if not k.endswith('__installed_cost_osbl')}
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


def _ng_exact(d, df):
    """price x intensity reproduces the NG channel at every point and for the
    coefficients."""
    pc = d.point_contrib
    np.testing.assert_allclose(d.prices['ng'] * pc['q__ng'], pc['m__natural_gas'], rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(d.channel_coef['natural_gas'], d.prices['ng'] * d.intensity_coef['ng'],
                               rtol=1e-9, atol=1e-14)


def test_ng_v2_pools_override_unreliable_adapter_columns():
    # SYNA: NG is a stream utility (charged at biosteam.stream_prices) and the
    # dryer burns unpriced methane (material pool, cost 0) -- like TAL.
    df, ad = _synthetic_sample('SYNA', v2=True)
    ref = physics.decompose_sample(_synthetic_sample('SYNA')[0])      # v1, correct adapter columns
    assert (df['ng_flow_material'] > 0).all() and (df['ng_cost_material'] == 0).all()
    np.testing.assert_allclose(df['ng_cost_utility'], df['natural_gas_flow'] * df['natural_gas_price'])
    # adapter columns as unreliable as TAL's: unpriced dryer gas counted at a wrong price
    bad = df.assign(natural_gas_flow=df['natural_gas_flow'] + df['ng_flow_material'],
                    natural_gas_price=3 * df['natural_gas_price'])
    assert not physics.decompose_sample(bad, ng_pool='utility').diagnostics['ng_pool_admissible']
    d = physics.decompose_sample(bad)
    dg = d.diagnostics
    assert dg['ng_pool'] == 'split' and dg['ng_pool_admissible'] and dg['ng_price_source'] == 'pools'
    assert dg['ng_mixed_pool'] == '' and dg['ng_unpriced_flow_share'] > 0
    assert dg['ng_utility_cost_share'] == pytest.approx(1.0)
    # same as v1 with correct adapter columns (capital channels differ only by
    # the v2 FCI-contribution weights)
    for ch in physics.OPEX_CHANNELS + ('fixed', 'closure'):
        np.testing.assert_allclose(d.channel_coef[ch], ref.channel_coef[ch], rtol=1e-9,
                                   atol=1e-12 * np.max(np.abs(d.coef)))
    cap = lambda x: sum(x.channel_coef[c] for c in physics.CAPEX_CHANNEL_AREAS)
    np.testing.assert_allclose(cap(d), cap(ref), rtol=1e-9, atol=1e-12 * np.max(np.abs(d.coef)))
    assert d.diagnostics['capex_weights'] == 'tdc' and ref.diagnostics['capex_weights'] == 'installed'
    # priced NG only (the utility gas), at the price actually charged
    assert d.prices['ng'] == pytest.approx(ad.values['p_natural_gas_price'], rel=1e-12)
    np.testing.assert_allclose(d.intensity_coef['ng'], ref.intensity_coef['ng'], rtol=1e-9, atol=1e-14)
    np.testing.assert_allclose(d.point_contrib['q__ng'], df['natural_gas_flow'] / df['F_final_product'],
                               rtol=1e-12)
    _ng_exact(d, df)
    # explicit legacy override still possible; 'split' without the columns fails
    assert physics.decompose_sample(bad, ng_pool='utility').diagnostics['ng_pool'] == 'utility'
    with pytest.raises(ValueError):
        physics.decompose_sample(bad.drop(columns=['ng_cost_material']), ng_pool='split')
    # incomplete v2 columns -> v1 fallback
    part = df.assign(ng_cost_utility=np.where(np.arange(len(df)) == 0, np.nan, df['ng_cost_utility']))
    assert physics.decompose_sample(part).diagnostics['ng_pool'] == 'utility'


def _ng_frame(n=12, seed=3):
    """Minimal per-point NG data (pools only) for natural_gas_pools tests."""
    rng = np.random.default_rng(seed)
    return pd.DataFrame({c: np.zeros(n) for c in physics.NG_POOL_COLUMNS}), rng


def test_natural_gas_pools_mixed_and_two_prices():
    df, rng = _ng_frame()
    n = len(df)
    # material pool mixes a priced CH4 feed (0.25 USD/kg) and an unpriced one;
    # one point is free of unpriced gas -> the priced flow is exact
    Fp, U = rng.uniform(100, 1000, n), rng.uniform(1e3, 5e4, n)
    U[4] = 0.
    df['ng_flow_material'], df['ng_cost_material'] = Fp + U, 0.25 * Fp
    ng = physics.natural_gas_pools(df)
    assert ng['diagnostics']['ng_mixed_pool'] == 'material'
    assert ng['price'] == pytest.approx(0.25) and ng['diagnostics']['ng_price_material'] == pytest.approx(0.25)
    np.testing.assert_allclose(ng['q_priced'], Fp, rtol=1e-12)
    np.testing.assert_allclose(ng['price'] * ng['q_equiv'], ng['cost'], rtol=1e-12)
    assert ng['diagnostics']['ng_unpriced_flow_share'] == pytest.approx(U.sum() / (U + Fp).sum())
    # unpriced gas at every point: price is a lower bound, price x q still exact
    U[4] = 10.
    df['ng_flow_material'] = Fp + U
    ng = physics.natural_gas_pools(df)
    assert ng['price'] < 0.25 and np.all(ng['q_priced'] >= Fp)
    np.testing.assert_allclose(ng['price'] * ng['q_equiv'], ng['cost'], rtol=1e-12)
    # two priced pools at different prices: effective price = cost / priced flow
    Fu = rng.uniform(100, 1000, n)
    df['ng_flow_utility'], df['ng_cost_utility'] = Fu, 0.20 * Fu
    df['ng_flow_material'], df['ng_cost_material'] = Fp, 0.30 * Fp
    ng = physics.natural_gas_pools(df)
    assert ng['diagnostics']['ng_mixed_pool'] == ''
    assert ng['price'] == pytest.approx((0.2 * Fu.sum() + 0.3 * Fp.sum()) / (Fu.sum() + Fp.sum()))
    np.testing.assert_allclose(ng['q_priced'], Fu + Fp, rtol=1e-12)
    np.testing.assert_allclose(ng['price'] * ng['q_equiv'], ng['cost'], rtol=1e-12)
    assert ng['diagnostics']['ng_q_equiv_rel_dev'] > 0          # cost-equivalent != physical flow
    assert ng['diagnostics']['ng_utility_cost_share'] == pytest.approx(0.2 * Fu.sum() / ng['cost'].sum())
    # no priced NG at all: channel and intensity 0, adapter price as fallback
    df[:] = 0.
    df['ng_flow_material'] = U
    ng = physics.natural_gas_pools(df, fallback_price=0.276)
    assert ng['price'] == 0.276 and ng['diagnostics']['ng_price_source'] == 'adapter'
    assert not ng['cost'].any() and not ng['q_equiv'].any()
    assert ng['diagnostics']['ng_unpriced_flow_share'] == 1.0
    assert np.isnan(physics.natural_gas_pools(df)['price'])
    assert physics.natural_gas_pools(df.assign(ng_cost_utility=np.nan)) is None


def test_ng_v2_no_priced_gas_in_decomposition():
    df, _ = _synthetic_sample('SYNA', v2=True)
    # all gas unpriced (e.g. HP, oxalic): utility pool cost moved out of the TEA
    df2 = df.assign(utility_cost=df['utility_cost'] - df['ng_cost_utility'] * df['operating_hours'],
                    MPSP=df['MPSP'] - df['ng_cost_utility'] / df['F_final_product'],
                    ng_cost_utility=0.)
    d = physics.decompose_sample(df2)
    assert d.diagnostics['ng_price_source'] == 'adapter'
    assert np.allclose(d.channel_coef['natural_gas'], 0) and np.allclose(d.intensity_coef['ng'], 0)
    assert d.prices['ng'] == pytest.approx(float(np.median(df['natural_gas_price'])))
    assert d.diagnostics['ng_pool_admissible'] and d.closure_rel < 1e-12


@pytest.mark.parametrize('name,source,mult', [('SYNA', 'labor', 1.9), ('SYNB', 'labor', 1.9),
                                              ('SYNC', 'labor_foc', 1.6)])
def test_fixed_channel_from_labor(name, source, mult):
    df, ad = _synthetic_sample(name, v2=True)
    d = physics.decompose_sample(df)
    dg = d.diagnostics
    F_true = synthetic.fixed_cost(ad.bio, ad.values)
    assert dg['fixed_source'] == source and dg['labor_multiplier'] == pytest.approx(mult, rel=1e-9)
    assert d.F_fixed == pytest.approx(F_true, rel=1e-9)
    assert d.F_fixed == pytest.approx(mult * dg['labor_cost'], rel=1e-9)
    assert d.phi == pytest.approx(synthetic.capital_charge(ad.values['p_IRR'], ad.values['p_income_tax']),
                                  rel=1e-9)
    assert not dg['F_negative'] and not dg['fixed_exceeds_foc'] and d.closure_rel < 1e-12
    # robust to DCF weights != 1 (the legacy free intercept is not)
    E, _ = physics.opex_channels(df)
    extra = 0.0067 * (E['feedstock'] + E['materials_other'] + E['utilities_other'])
    df2 = df.assign(MPSP=df['MPSP'] + extra / (df['F_final_product'] * df['operating_hours']))
    assert physics.decompose_sample(df2).F_fixed == pytest.approx(F_true, rel=1e-9)
    # recorded labor above the TEA's FOC is inconsistent -> flagged
    hi = df.assign(labor_cost=10 * df['FOC'])
    assert physics.decompose_sample(hi).diagnostics['fixed_exceeds_foc'] or source == 'labor_foc'
    # without labor data: legacy FOC intercept (exact here: FOC affine in FCI)
    d0 = physics.decompose_sample(df.drop(columns=['labor_cost', 'labor_burden']))
    assert d0.diagnostics['fixed_source'] == 'foc_intercept'
    assert d0.F_fixed == pytest.approx(F_true, rel=1e-6)


def test_labor_fixed_cost_varying_labor():
    # labor scaled with plant throughput (e.g. lactic) and no labor_burden:
    # FOC = 0.02*K + 1.6*labor_i is recovered point by point
    K = np.linspace(1e8, 3e8, 8)
    labor = np.linspace(2e6, 2.6e6, 8)[::-1]
    df = pd.DataFrame({'labor_cost': labor})
    F, src, dg = physics.labor_fixed_cost(df, K, 0.02 * K + 1.6 * labor)
    assert src == 'labor_foc' and dg['labor_multiplier'] == pytest.approx(1.6)
    np.testing.assert_allclose(F, 1.6 * labor, rtol=1e-10)
    # implausible multiplier (< 1, i.e. negative burden) -> no labor-based F
    assert physics.labor_fixed_cost(df, K, 0.02 * K + 0.5 * labor)[0] is None
    assert physics.labor_fixed_cost(df.assign(labor_cost=np.nan), K, K)[0] is None


def test_foc_maintenance_allocated_by_isbl():
    # cellulosic-TEA FOC: labor*(1+burden) + 0.007*FCI + 0.03*ISBL installed
    df, ad = _synthetic_sample('SYNB', v2=True)
    osbl = sum(df[f'{a}__installed_cost_osbl'] for a in physics.AREAS)
    S = df['installed_cost'] - osbl
    F = df['labor_cost'] * (1 + df['labor_burden'])
    foc_old = df['FOC']
    foc = F + 0.007 * df['FCI'] + 0.03 * S
    denom = df['F_final_product'] * df['operating_hours']
    df2 = df.assign(FOC=foc, MPSP=df['MPSP'] + (foc - foc_old) / denom)
    d = physics.decompose_sample(df2)
    dg = d.diagnostics
    assert dg['foc_alloc'] == 'isbl' and dg['foc_isbl_rate'] == pytest.approx(0.03, rel=1e-9)
    assert dg['capex_weights'] == 'tdc' and d.closure_rel < 1e-12 and d.sum_check < 1e-11
    phi_K = dg['phi_capital']
    assert phi_K == pytest.approx(synthetic.capital_charge(ad.values['p_IRR'], ad.values['p_income_tax'])
                                  - 0.03, rel=1e-9)
    pc = d.point_contrib
    w = {a: df[f'{a}__installed_cost'] * (1 if a in synthetic.OSBL_AREAS else 1.175) for a in physics.AREAS}
    W = sum(w.values())
    base = (phi_K + 0.007) * df['FCI'] / W
    # OSBL area: capital + insurance only; ISBL area: + maintenance on its ISBL cost
    np.testing.assert_allclose(pc['m__capex_wastewater'] * denom, base * w['wastewater'], rtol=1e-9)
    np.testing.assert_allclose(pc['m__capex_fermentation'] * denom,
                               base * w['fermentation'] + 0.03 * df['fermentation__installed_cost'], rtol=1e-9)
    cap = sum(pc[f'm__{c}'] for c in physics.CAPEX_CHANNEL_AREAS) * denom
    np.testing.assert_allclose(cap, phi_K * df['FCI'] + foc - F, rtol=1e-10)


def test_fixed_from_foc_is_robust_to_dcf_weights():
    # VOC carried with a DCF weight != 1 (startup fractions): the legacy free
    # intercept absorbs it, the FOC-based fixed channel stays at the labor cost
    df, ad = _synthetic_sample('SYNB')
    E, info = physics.opex_channels(df)
    extra = 0.0067 * (E['feedstock'] + E['materials_other'] + E['utilities_other'])
    df2 = df.assign(MPSP=df['MPSP'] + extra / (df['F_final_product'] * df['operating_hours']))
    F_true = synthetic.fixed_cost(ad.bio, ad.values)
    d = physics.decompose_sample(df2)
    assert d.diagnostics['fixed_basis'] == 'FOC'
    assert d.F_fixed == pytest.approx(F_true, rel=1e-9) and d.closure_rel > 0
    d_old = physics.decompose_sample(df2, fixed_basis='intercept')
    assert abs(d_old.F_fixed / F_true - 1) > 1e-3
    assert d.sum_check < 1e-11


def test_capex_weights_tdc_contributions():
    # with the OSBL part of each area recorded, ISBL areas carry (1 + f) more
    # capital per installed USD; synthetic FCI = 1.6*(I + 0.175*ISBL)
    df, _ = _synthetic_sample('SYNB')
    osbl = {a: (df[f'{a}__installed_cost'] if a in ('wastewater', 'chp', 'facilities') else 0. * df['y'])
            for a in physics.AREAS}
    df2 = df.assign(**{f'{a}__installed_cost_osbl': v for a, v in osbl.items()})
    d = physics.decompose_sample(df2)
    assert d.diagnostics['capex_weights'] == 'tdc'
    assert d.diagnostics['isbl_factor'] == pytest.approx(0.175, rel=1e-6)
    pc = d.point_contrib
    r_isbl = pc['m__capex_fermentation'] / df['fermentation__installed_cost']
    r_osbl = pc['m__capex_wastewater'] / df['wastewater__installed_cost']
    np.testing.assert_allclose(r_isbl / r_osbl, 1.175, rtol=1e-6)
    cap = sum(pc[f'm__{c}'] for c in physics.CAPEX_CHANNEL_AREAS)
    np.testing.assert_allclose(cap * df['F_final_product'] * df['operating_hours'], d.phi * df['FCI'], rtol=1e-9)
    assert physics.decompose_sample(df).diagnostics['capex_weights'] == 'installed'
