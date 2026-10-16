# -*- coding: utf-8 -*-
"""
Small synthetic tables in the processed-data format of DESIGN.md section 7
(``samples``, ``points``, ``configs``) for fast tests of the model code.

The generator is physically plausible but simple: channel coefficients are
built directly (price-explicit channels as price x intensity), totals are
their sums, MPSP points follow the unifying equation with 0.5 % noise,
fermenter capital decreases with productivity, and the feedstock intensity
``b'`` equals ``1/(w_sugar * eta)`` with ``eta < 1``.
"""
import numpy as np
import pandas as pd

K = ('a', 'b', 'c', 'd')
CHANNELS = ('feedstock', 'electricity', 'natural_gas', 'neutralization',
            'materials_other', 'utilities_other', 'capex_upstream',
            'capex_fermentation', 'capex_downstream', 'capex_wastewater',
            'capex_utilities', 'fixed', 'coproducts', 'closure')

FEEDSTOCKS = {
    #  name: (type, w_sugar_eq, w_water, w_lignin, price, capacity kg/h, eta)
    'glucose': ('pure_sugar', 0.95, 0.05, 0.0, 0.40, 3.0e4, 0.97),
    'corn': ('starch', 0.62, 0.15, 0.0, 0.20, 1.0e5, 0.93),
    'sugarcane': ('sugar_juice', 0.14, 0.70, 0.05, 0.035, 4.0e5, 0.90),
    'cornstover': ('lignocellulose', 0.50, 0.20, 0.15, 0.08, 1.0e5, 0.82),
}

PRODUCTS = {
    # label: (fp chem, final, final chem, ty, stoich ratio, neutralization,
    #         n_acid, pKa1, upgrading, aerobic, solubility, separation,
    #         r0, g_frac)
    'HP': ('HP', 'AcrylicAcid', 'AcrylicAcid', 1.0, 0.80, True, 1, 4.51, True, True,
           1000., 'evaporation+crystallization', 0.5, 0.90),
    'LA': ('LacticAcid', 'LA', 'LacticAcid', 1.0, 1.0, True, 1, 3.86, False, False,
           1000., 'esterification', 0.9, 0.92),
    'SA': ('SuccinicAcid', 'SA', 'SuccinicAcid', 1.12, 1.0, True, 2, 4.2, False, False,
           83., 'crystallization', 0.6, 0.88),
    'EtOH': ('Ethanol', 'EtOH', 'Ethanol', 0.511, 1.0, False, 0, np.nan, False, False,
             1000., 'distillation', 1.5, 0.97),
    'TAL': ('TAL', 'SorbicAcid', 'SorbicAcid', 0.467, 0.889, False, 0, np.nan, True, True,
            8.6, 'crystallization', 0.2, 0.85),
    'NewAcid': ('NewAcid', 'NewAcid', 'NewAcid', 0.58, 1.0, True, 1, 4.6, False, True,
                300., 'crystallization', 0.4, 0.90),
}

# fp_* columns as the harness would record them (registry overrides known IDs)
_FP_COLS = {
    'HP': dict(MW=90.08, n_C=3, n_H=6, n_O=3, n_N=0, Tb=453., Tm=345., Psat_ferm=180.,
               LHV=1.22e6, Hvap_298=5.6e4),
    'LacticAcid': dict(MW=90.08, n_C=3, n_H=6, n_O=3, n_N=0, Tb=505., Tm=290., Psat_ferm=20.,
                       LHV=1.25e6, Hvap_298=6.5e4),
    'SuccinicAcid': dict(MW=118.09, n_C=4, n_H=6, n_O=4, n_N=0, Tb=508., Tm=461., Psat_ferm=1e-3,
                         LHV=1.36e6, Hvap_298=8.0e4),
    'Ethanol': dict(MW=46.07, n_C=2, n_H=6, n_O=1, n_N=0, Tb=351.4, Tm=159., Psat_ferm=1.4e4,
                    LHV=1.24e6, Hvap_298=4.2e4),
    'TAL': dict(MW=126.11, n_C=6, n_H=6, n_O=3, n_N=0, Tb=560., Tm=461., Psat_ferm=0.5,
                LHV=2.6e6, Hvap_298=7.0e4),
    'NewAcid': dict(MW=104.1, n_C=4, n_H=8, n_O=3, n_N=0, Tb=480., Tm=320., Psat_ferm=30.,
                    LHV=1.8e6, Hvap_298=6.8e4),
}
_FP_COLS.update({'AcrylicAcid': dict(MW=72.06, n_C=3, n_H=4, n_O=2, n_N=0, Tb=414., Tm=286.,
                                     Psat_ferm=1.0e3, LHV=1.3e6, Hvap_298=5.3e4),
                 'SorbicAcid': dict(MW=112.13, n_C=6, n_H=8, n_O=2, n_N=0, Tb=501., Tm=407.,
                                    Psat_ferm=2., LHV=3.0e6, Hvap_298=7.0e4)})


def _channels(fs, prod, params, r, rng):
    ftype, w, w_water, w_lig, _, cap, eta = FEEDSTOCKS[fs]
    (_, _, _, ty, ratio, neut, n_acid, _, upg, aer, sol, _, _, _) = PRODUCTS[prod]
    lig = float(ftype == 'lignocellulose')
    phi = 0.12 + 0.9 * (params['p_IRR'] - 0.10)
    capf = (cap / 1e5) ** -0.3
    q = {
        'fs': np.array([0.002, 1. / (w * eta), 0., 0.]),
        'el': np.array([0.3 + 0.6 * aer - 1.2 * lig, 0.1 - 0.3 * lig, 3.0, 0.]),
        'ng': np.array([0.05, 0.02, 2.0 * (1 - 0.8 * lig), 0.]),
    }
    ch = {
        'neutralization': np.array([0.25 * n_acid * neut, 0., 0., 0.]),
        'materials_other': np.array([0.05 + 0.1 * lig, 0.03, 0., 0.]),
        'utilities_other': np.array([0.02, 0., 0.8, 0.]),
        'capex_upstream': phi * capf * np.array([0.1, 0.25 * (1 + lig), 0., 0.]),
        'capex_fermentation': phi * capf * np.array([0.8 / r * (1 + aer), 0.05, 0.5, 0.]),
        'capex_downstream': phi * capf * np.array([0.3 + 0.2 * upg, 0., 6.0 * (100. / sol) ** 0.2, 0.1]),
        'capex_wastewater': phi * capf * np.array([0.05, 0.05, 1.0, 0.]),
        'capex_utilities': phi * capf * np.array([0.1, 0.1, 0.5, 0.]),
        'fixed': np.array([0., 0.08 * (cap / 1e5) ** -0.5, 0., 0.]),
        'coproducts': np.array([0., -0.02 * lig, 0., 0.]),
        'closure': rng.normal(0., 0.003, 4),
    }
    for key, pcol in (('fs', 'feedstock_price'), ('el', 'electricity_price'), ('ng', 'natural_gas_price')):
        ch[{'fs': 'feedstock', 'el': 'electricity', 'ng': 'natural_gas'}[key]] = params[pcol] * q[key]
    return q, ch, phi


def make_fake_tables(n_samples_per_config=6, products=None, feedstocks=None, seed=0,
                     n_y=5, n_t=5, noise=0.005):
    """Return ``(samples, points, configs)`` DataFrames in the section 7 format."""
    rng = np.random.RandomState(seed)
    products = list(products or PRODUCTS)
    feedstocks = list(feedstocks or FEEDSTOCKS)
    srows, prows = [], []
    for prod in products:
        (fpc, final, ppc, ty, ratio, neut, n_acid, pKa, upg, aer, sol, sep, r0, gfrac) = PRODUCTS[prod]
        for fs in feedstocks:
            ftype, w, w_water, w_lig, fprice, cap, eta = FEEDSTOCKS[fs]
            cid = f'{prod}_{fs}'
            base = dict(feedstock_price=fprice, electricity_price=0.07, natural_gas_price=0.25,
                        p_IRR=0.10, p_income_tax=0.21, p_operating_days=330.,
                        p_boiler_efficiency=0.8, p_feedstock_capacity=cap)
            bp0 = rng.uniform(0, 1)
            for i in range(n_samples_per_config):
                if i == 0:
                    kind, params, rm = 'baseline', dict(base), 1.
                elif i in (1, 2):
                    kind, params, rm = 'productivity', dict(base), (0.5, 2.0)[i - 1]
                else:
                    kind, rm = 'random', float(np.exp(rng.uniform(np.log(0.4), np.log(2.5))))
                    params = dict(base)
                    params.update(feedstock_price=fprice * rng.uniform(0.6, 1.5),
                                  electricity_price=rng.uniform(0.03, 0.11),
                                  natural_gas_price=rng.uniform(0.15, 0.40),
                                  p_IRR=rng.uniform(0.07, 0.13), p_income_tax=rng.uniform(0.15, 0.3),
                                  p_operating_days=rng.uniform(300, 350),
                                  p_boiler_efficiency=rng.uniform(0.7, 0.85),
                                  p_feedstock_capacity=cap * rng.uniform(0.7, 1.3))
                if prod == 'EtOH':
                    params['p_boiler_efficiency'] = np.nan  # not mapped by this adapter
                r = r0 * rm
                q, ch, phi = _channels(fs, prod, params, r, rng)
                g = ratio * float(np.clip(gfrac + rng.normal(0, 0.01), 0.5, 0.99))
                tot = sum(ch[c] for c in CHANNELS)
                y_min, y_max = 0.1 * ty, 0.8 * ty
                t_min, t_max = 10., min(150., 0.5 * sol + 10.)
                ys, ts = np.linspace(y_min, y_max, n_y), np.linspace(t_min, t_max, n_t)
                Y, T = np.meshgrid(ys, ts)
                Y, T = Y.ravel(), T.ravel()
                keep = T <= 300. * Y
                Y, T = Y[keep], T[keep]
                m = (tot[0] + tot[1] / Y + tot[2] / T + tot[3] / (Y * T)) / g
                M = m * (1 + rng.normal(0, noise, len(m)))
                sid = f'{cid}_{i}'
                for yy, tt, mm in zip(Y, T, M):
                    prows.append(dict(sample_id=sid, y=yy, t=tt, MPSP=mm,
                                      recovery=g * (1 + rng.normal(0, 0.003))))
                row = dict(sample_id=sid, config_id=cid, biorefinery=prod, sample_kind=kind, seed=i,
                           feedstock=fs, feedstock_type=ftype, fermentation_product=prod,
                           final_product=final, ferm_product_chemical=fpc, final_product_chemical=ppc,
                           neutralization=neut, neutralizing_agent='Ca(OH)2' if neut else '',
                           n_acid_groups=n_acid, pKa1=pKa, final_product_stoich_mass_ratio=ratio,
                           upgrading=upg, separation_variant='', separation_principle=sep,
                           aerobic=aer, fermentation_mode='batch' if prod != 'EtOH' else 'fed-batch',
                           T_fermentation=305.15, product_water_solubility=sol,
                           theoretical_yield=ty, inflection_yield=y_max, y_min=y_min, y_max=y_max,
                           t_min=t_min, t_max=t_max, baseline_productivity=r0, r=r, r_mult=rm,
                           fs_w_glucose=w * 0.5, fs_w_glucan=w * 0.45, fs_w_water=w_water,
                           fs_w_lignin=w_lig, fs_w_fermentable_sugar_eq=w, fs_price=fprice,
                           fs_F_mass=cap, bp_strain_factor=bp0 if kind != 'random' else rng.uniform(0, 1),
                           a=tot[0], b=tot[1], c=tot[2], d=tot[3], g=g, R2=0.999,
                           n_points=int(len(M)),
                           fit_mape=float(np.mean(np.abs(m - M) / np.abs(M))),
                           phi=phi, F_fixed=1e6, closure_rel=0.002)
                for k, v in _FP_COLS[fpc].items():
                    row[f'fp_{k}'] = v
                for k, v in _FP_COLS.get(ppc, _FP_COLS[fpc]).items():
                    row[f'pp_{k}'] = v
                row.update({k: v for k, v in params.items()})
                row['p_feedstock_price'] = params['feedstock_price']
                row['p_electricity_price'] = params['electricity_price']
                row['p_natural_gas_price'] = params['natural_gas_price']
                for c in CHANNELS:
                    for j, k in enumerate(K):
                        row[f'ch__{c}__{k}'] = ch[c][j]
                for key in ('fs', 'el', 'ng'):
                    for j, k in enumerate(K):
                        row[f'int__{key}__{k}'] = q[key][j]
                srows.append(row)
    samples = pd.DataFrame(srows)
    points = pd.DataFrame(prows)
    spec_cols = [c for c in samples.columns if c in (
        'config_id', 'biorefinery', 'feedstock', 'feedstock_type', 'fermentation_product',
        'final_product', 'ferm_product_chemical', 'final_product_chemical', 'neutralization',
        'neutralizing_agent', 'n_acid_groups', 'pKa1', 'final_product_stoich_mass_ratio',
        'upgrading', 'separation_variant', 'separation_principle', 'aerobic',
        'fermentation_mode', 'T_fermentation', 'product_water_solubility',
        'theoretical_yield', 'inflection_yield') or c.startswith(('fp_', 'pp_', 'fs_'))]
    configs = samples[spec_cols].drop_duplicates('config_id').reset_index(drop=True)
    return samples, points, configs
