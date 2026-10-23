#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Quickstart of the coefficient predictor: the examples of
``coeff_predictor/README.md``, end to end.

    PYTHONPATH=/path/to/fermentation_insights python \\
        fermentation_insights/coeff_predictor/examples/quickstart.py [CHECKPOINT] [--no-thermosteam]

``CHECKPOINT`` is a directory written by ``train.Ensemble.save`` (default:
``models/release`` of the package, produced by ``experiments.py``; a relative
path is looked up in the current directory, then in the package directory).
``--no-thermosteam`` uses the manual descriptor dict of the new product instead
of computing it with thermosteam.

Examples
--------
1. a known configuration, specified by names;
2. a new feedstock (custom composition and price);
3. a new fermentation product (descriptors from thermosteam or a manual dict,
   plus its regime and downstream fields);
4. the MPSP surface and the relative impact of yield vs. titer;
5. uncertainty intervals and the applicability domain;
6. the cost-channel breakdown and the exact price linearity;
7. a user-supplied recovery g.
"""
import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd

from fermentation_insights.coeff_predictor.predict import CoefficientPredictor, ApplicabilityWarning
from fermentation_insights.coeff_predictor.features import W_SUGAR, describe_new_chemical

PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_CHECKPOINT = os.path.join('models', 'release')
COEF = ['a', 'b', 'c', 'd', 'g']

#: hydrolysis mass gain to glucose-eq. sugar (as the data-generation harness,
#: ``data_generation/harness.describe_feedstock``); monomers count 1:1
SUGAR_GAIN = {'glucose': 1., 'xylose': 1., 'arabinose': 1., 'galactose': 1., 'mannose': 1., 'fructose': 1.,
              'glucan': 180.156 / 162.141, 'starch': 180.156 / 162.141, 'galactan': 180.156 / 162.141,
              'mannan': 180.156 / 162.141, 'xylan': 150.13 / 132.115, 'arabinan': 150.13 / 132.115,
              'sucrose': 2 * 180.156 / 342.297}

# %% Example specifications

#: 1. a known configuration: names fill the composition, yields, stoichiometric
#: ratio, the regime fields that are unique for HP -> acrylic acid in training
#: and the contextual baselines; ``config_id`` fills the downstream fields
#: (downstream.TABLE); the regime fields that differ between the HP designs
#: (neutral vs. low pH) are given explicitly
KNOWN = dict(
    config_id='HP_neutral_cornstover',
    feedstock='cornstover', fermentation_product='HP', final_product='AcrylicAcid',
    neutralization=True, neutralizing_agent='Ca(OH)2',
    r=0.28,                       # g/L/h
    feedstock_price=0.0517,       # USD/kg WET corn stover
)

#: 2. wheat straw (illustrative dry-basis composition, 15 % moisture)
WHEAT_STRAW_DRY = dict(glucan=0.354, xylan=0.197, arabinan=0.024, galactan=0.008, mannan=0.003,
                       lignin=0.170, ash=0.065, protein=0.030, extractives=0.120)
WHEAT_STRAW_MOISTURE = 0.15

#: 3. itaconic acid (C5H6O4): manual descriptors (keys of registry/chemicals.json;
#: dor / Ymax_* are completed from the atom counts and MW). Values from
#: thermosteam (descriptors.chemical_descriptors('Itaconic acid')).
ITACONIC = dict(MW=130.099, Tb=601.0, Tm=444.15, Hvap_298=1.1407e5, Psat_ferm=1.98e-3, LHV=1.8517e6,
                n_C=5, n_H=6, n_O=4, n_N=0, gamma_inf=7.78, log10_alpha_inf=-5.61)
NEW_PRODUCT = dict(
    feedstock='glucose', feedstock_price=0.28, r=0.6,
    fermentation_product='ItaconicAcid', final_product='ItaconicAcid',
    ferm_product_chemical='ItaconicAcid', final_product_chemical='ItaconicAcid',
    theoretical_yield=0.722,      # g/g glucose (1 mol/mol)
    inflection_yield=0.58,        # g/g; end of the fitted yield domain
    final_product_stoich_mass_ratio=1.0, upgrading=False,
    n_acid_groups=2, pKa1=3.85, product_water_solubility=83.,   # g/L
    neutralization=False, neutralizing_agent='none',
    aerobic=True, fermentation_mode='fed-batch', T_fermentation=310.15,
    # final product and downstream process (DESIGN.md section 5)
    final_product_mass_fraction=0.99, final_product_is_broth=False, final_product_phase='solid',
    water_removal_route='crystallize',
    separation_ops='centrifugation+evaporation+crystallization+filtration+drying',
    upgrading_ops='none', coproducts='electricity', neutralization_degree=0., acid_regeneration='none',
)


# %% Helpers

def section(title):
    print('\n' + '=' * 78 + f'\n{title}\n' + '=' * 78)


def fmt(df):
    return df.to_string(float_format=lambda v: f'{v:.4g}')


def fermentable_sugar_eq(comp):
    """Glucose-eq. fermentable sugar mass fraction of an ``fs_w_*`` composition."""
    return sum(comp.get(f'fs_w_{k}', 0.) * gain for k, gain in SUGAR_GAIN.items())


def wet_composition(P, dry, moisture):
    """``fs_w_*`` fields (wet basis) of a feedstock from its dry-basis
    composition: every composition column of the checkpoint is set (absent
    components 0, so nothing is imputed), the rest goes to ``fs_w_other``."""
    comp = {d['name']: 0. for d in P.spec.numeric if d['name'].startswith('fs_w_')}
    comp.update({f'fs_w_{k}': v * (1. - moisture) for k, v in dry.items()})
    comp['fs_w_water'] = moisture
    comp['fs_w_other'] = max(0., 1. - sum(v for k, v in comp.items() if k not in ('fs_w_other', W_SUGAR)))
    comp[W_SUGAR] = fermentable_sugar_eq(comp)
    return comp


def mpsp_eq(p, y, t, g=None):
    """MPSP [USD/kg final product] from predicted coefficients (one row)."""
    a, b, c, d = (float(p[k]) for k in 'abcd')
    return (a + b / y + c / t + d / (y * t)) / float(p['g'] if g is None else g)


def relative_impact(p, y, t, dy=0.01, dt=1.):
    """MPSP decrease from +dy yield [g/g] relative to that from +dt titer [g/L]
    (publication, MPSP_YT_fit_only.py): (dMPSP/dy * dy) / (dMPSP/dt * dt); g cancels."""
    b, c, d = (float(p[k]) for k in 'bcd')
    return (dy * (b + d / t) / y ** 2) / (dt * (c + d / y) / t ** 2)


def resolve_checkpoint(path):
    for p in (path, os.path.join(PKG_DIR, path)):
        if os.path.isfile(os.path.join(p, 'meta.json')):
            return p
    sys.exit(f'no checkpoint at {path!r} (nor {os.path.join(PKG_DIR, path)!r}); train one with '
             'experiments.py (models/release) or pass a checkpoint directory')


# %% Examples

def example_known(P):
    section('1. Known configuration by name (HP -> acrylic acid, neutral, corn stover)')
    p = P.predict(KNOWN)
    print(fmt(p[COEF + ['g_source', 'in_domain']]))
    filled = p['filled_from_registry'].iloc[0].split(', ')
    print(f'filled from the registry ({len(filled)} fields):', ', '.join(filled[:12]), '...')
    print('contextual defaults:', p['contextual_defaults'].iloc[0] or '-')
    print('imputed:', p['imputed'].iloc[0] or '-')
    A = P.unprimed(KNOWN, warn=False)
    print("MPSP-level coefficients A, B, C, D = a'/g, ...:\n" + fmt(A))
    return p


def example_new_feedstock(P, p_known):
    section('2. New feedstock: wheat straw (custom composition and price)')
    comp = wet_composition(P, WHEAT_STRAW_DRY, WHEAT_STRAW_MOISTURE)
    spec = dict(KNOWN, feedstock_type='lignocellulose', feedstock_price=0.0595, p_feedstock_capacity=1.0e5,
                **comp)
    spec.pop('feedstock')          # a new feedstock needs no name (it is not a model input)
    # config_id with a composition of its own and no feedstock name: the
    # configuration's process, product and downstream design with a NEW
    # feedstock (its composition is never filled from corn stover)
    print(f'fs_w_fermentable_sugar_eq = {comp[W_SUGAR]:.4f} kg glucose-eq. / kg wet straw')
    p = P.predict(spec)
    out = pd.concat([p_known[COEF + ['ch__feedstock__b']], p[COEF + ['ch__feedstock__b']]])
    out.index = ['corn stover', 'wheat straw']
    print(fmt(out))
    print(f"wheat straw: intervals_calibrated={bool(p['intervals_calibrated'].iloc[0])} (a feedstock not in "
          f"training), in_domain={bool(p['in_domain'].iloc[0])} (outside range: "
          f"{p['outside_range'].iloc[0] or '-'})")
    return spec


def example_new_product(P, use_thermosteam=True):
    section('3. New fermentation product: itaconic acid from glucose')
    desc = ITACONIC
    if use_thermosteam:
        try:
            # thermosteam (optional): computes and registers the descriptors for this session
            desc = describe_new_chemical('ItaconicAcid', 'Itaconic acid')
            print('descriptors from thermosteam:',
                  {k: round(desc[k], 4) for k in ('MW', 'Tb', 'dor', 'Ymax_stoich', 'log10_alpha_inf')})
        except Exception as e:     # ImportError without thermosteam
            print(f'thermosteam unavailable ({type(e).__name__}): manual descriptor dict')
    descriptors = {'ItaconicAcid': desc}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        p = P.predict(NEW_PRODUCT, descriptors=descriptors)
    print(fmt(p[COEF + ['g_source', 'in_domain', 'ad_ratio']]))
    for w in caught:   # ApplicabilityWarning: outside the training domain (extrapolation)
        tag = 'applicability' if issubclass(w.category, ApplicabilityWarning) else 'note'
        print(f'  [{tag}] {str(w.message)[:300]}')
    return p, descriptors


def example_surface(P, p_known):
    section('4. MPSP surface and relative impact of yield vs. titer (from the predicted coefficients)')
    reg = P.spec.spec_registry
    infl = reg['fermentation_product']['HP']['inflection_yield']
    dom = reg['domain'].get('HP|AcrylicAcid') or reg['domain']['__all__']
    ys = np.round(np.linspace(0.25 * infl, infl, 5), 3)
    ts = np.round(np.geomspace(max(dom['t_min'], 10.), min(dom['t_max'], 200.), 5), 1)
    Y, T = np.meshgrid(ys, ts)
    p = p_known.iloc[0]
    M = pd.DataFrame(mpsp_eq(p, Y, T), index=pd.Index(ts, name='t [g/L]'), columns=pd.Index(ys, name='y [g/g]'))
    print('MPSP [USD/kg acrylic acid]:\n' + fmt(M))
    RI = pd.DataFrame(relative_impact(p, Y, T), index=M.index, columns=M.columns)
    print('relative impact of +0.01 g/g yield vs. +1 g/L titer (> 1: yield matters more):\n' + fmt(RI))
    # the same surface from the ensemble (mean of the members' MPSP) with intervals
    m = P.mpsp(KNOWN, y=Y.ravel(), t=T.ravel(), warn=False)
    cols = [c for c in ('y', 't', 'MPSP', 'MPSP_std', 'MPSP_lo', 'MPSP_hi') if c in m]
    print('ensemble MPSP at the corners:\n' + fmt(m.iloc[[0, len(ys) - 1, -len(ys), -1]][cols]))


def example_uncertainty(P, p_new):
    section('5. Uncertainty intervals and the applicability domain')
    conf = P.ensemble.conformal or {}
    print(f"conformal level {1 - conf.get('alpha', float('nan')):.2f}; calibration: {conf.get('n')} samples of "
          f"{conf.get('n_calib_groups')} groups ({conf.get('calibration_grouping')}); coverage guaranteed: "
          f"{conf.get('coverage_guaranteed')}")
    p = P.predict(KNOWN, warn=False)
    cols = [f'{k}{s}' for k in COEF for s in ('', '_std', '_lo', '_hi') if f'{k}{s}' in p]
    print(fmt(p[cols]))
    for name, q in (('known HP configuration', p), ('itaconic acid', p_new)):
        r = q.iloc[0]
        print(f"{name}: in_domain={r['in_domain']}, ad_ratio={r['ad_ratio']:.2f} (> 1: farther from the "
              f"training data than 99 % of its rows), outside_range={r['outside_range'] or '-'}")


def example_channels(P, p_known):
    section('6. Channel breakdown (primed basis, USD per kg fermentation product) and price linearity')
    p = p_known.iloc[0]
    names = P.ensemble.channel_names()
    tab = pd.DataFrame({k: [p[f'ch__{c}__{k}'] for c in names] for k in 'abcd'}, index=names)
    y, t = 0.5, 100.
    tab[f'MPSP share at y={y}, t={t}'] = (tab['a'] + tab['b'] / y + tab['c'] / t + tab['d'] / (y * t)) / p['g']
    tab.loc['total'] = tab.sum()
    print(fmt(tab))
    print(f"sum of channels = total: a' {tab.loc['total', 'a']:.6g} vs {p['a']:.6g}; "
          f"MPSP {tab.iloc[-1, -1]:.6g} vs {mpsp_eq(p, y, t):.6g} USD/kg")
    print(f"feedstock channel = price x intensity: b'_fs {p['ch__feedstock__b']:.6g} = "
          f"{KNOWN['feedstock_price']} x {p['int__fs__b']:.6g}")
    p2 = P.predict(dict(KNOWN, feedstock_price=2 * KNOWN['feedstock_price']), warn=False).iloc[0]
    d = {c: p2[f'ch__{c}__b'] - p[f'ch__{c}__b'] for c in names}
    print("doubling the feedstock price changes b' by", f"{p2['b'] - p['b']:.6g}:",
          ', '.join(f'{c} {v:+.3g}' for c, v in d.items() if abs(v) > 1e-9) or '-')


def example_recovery(P, descriptors):
    section('7. User-supplied recovery g (keeps the predicted a\', b\', c\', d\')')
    y, t = 0.5, 80.
    p0 = P.predict(NEW_PRODUCT, descriptors=descriptors, warn=False)
    p1 = P.predict(NEW_PRODUCT, descriptors=descriptors, warn=False, recovery=0.85)
    print(fmt(pd.concat([p0, p1])[['a', 'b', 'c', 'd', 'g', 'g_predicted', 'g_source']].set_index('g_source')))
    m0 = P.mpsp(NEW_PRODUCT, y=y, t=t, descriptors=descriptors, warn=False)
    m1 = P.mpsp(NEW_PRODUCT, y=y, t=t, descriptors=descriptors, warn=False, recovery=0.85)
    print(f"MPSP at y={y}, t={t}: predicted g {m0['MPSP'].iloc[0]:.4g} USD/kg, with g = 0.85: "
          f"{m1['MPSP'].iloc[0]:.4g} USD/kg")
    try:
        P.predict(NEW_PRODUCT, descriptors=descriptors, warn=False, recovery=1.2)
    except ValueError as e:
        print(f'recovery=1.2 -> ValueError: {e}')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('checkpoint', nargs='?', default=DEFAULT_CHECKPOINT,
                    help='checkpoint directory (default: models/release)')
    ap.add_argument('--no-thermosteam', action='store_true', help='use the manual descriptors of the new product')
    a = ap.parse_args(argv)
    warnings.formatwarning = lambda msg, cat, *args, **kw: f'  [{cat.__name__}] {msg}\n'
    path = resolve_checkpoint(a.checkpoint)
    P = CoefficientPredictor.load(path)
    print(f'checkpoint {path}: {len(P.ensemble.members)} members; downstream features: {P.spec.downstream}')
    p_known = example_known(P)
    example_new_feedstock(P, p_known)
    p_new, descriptors = example_new_product(P, not a.no_thermosteam)
    example_surface(P, p_known)
    example_uncertainty(P, p_new)
    example_channels(P, p_known)
    example_recovery(P, descriptors)
    return 0


if __name__ == '__main__':
    sys.exit(main())
