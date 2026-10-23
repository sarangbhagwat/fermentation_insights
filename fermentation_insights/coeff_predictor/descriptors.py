# -*- coding: utf-8 -*-
"""
Thermodynamic and stoichiometric descriptors of fermentation products and
final products.

Descriptors are computed from the thermosteam chemical database (the same
source at training and inference time) and stored in a JSON registry
(``registry/chemicals.json``) so that the predictor itself does not need
thermosteam. For a new chemical, call :func:`chemical_descriptors` (requires
thermosteam) or supply the descriptor values manually.

Descriptors
-----------
MW, Tb, Tm [K], Hvap_298 [J/mol], Psat_ferm [Pa] (at T_ferm), LHV [J/mol],
n_C, n_H, n_O, n_N,
dor : degree of reduction per carbon, (4C + H - 2O - 3N)/C (NH3 N-source),
Ymax_carbon, Ymax_electron : stoichiometric maximum yields on glucose [g/g]
    from carbon (no CO2 fixation) and available-electron balances,
Ymax_stoich : min of the two,
gamma_inf : infinite-dilution activity coefficient in water at T_ferm
    (Dortmund-modified UNIFAC; falls back to original UNIFAC),
log10_alpha_inf : log10 relative volatility vs. water at infinite dilution,
    alpha = gamma_inf * Psat_i / Psat_water.
"""
import json
import math
import os

__all__ = ('chemical_descriptors', 'load_registry', 'build_registry',
           'KNOWN_CHEMICALS', 'REGISTRY_PATH')

REGISTRY_PATH = os.path.join(os.path.dirname(__file__), 'registry', 'chemicals.json')

#: Known fermentation/final products: ID -> thermosteam search identifier.
KNOWN_CHEMICALS = {
    'Ethanol': 'Ethanol',
    'LacticAcid': 'Lactic acid',
    'HP': '3-Hydroxypropionic acid',
    'AcrylicAcid': 'Acrylic acid',
    'SuccinicAcid': 'Succinic acid',
    'TAL': '675-10-5',  # triacetic acid lactone (4-hydroxy-6-methyl-2-pyrone)
    'SorbicAcid': 'Sorbic acid',
    'PotassiumSorbate': '24634-61-5',
    'OxalicAcid': 'Oxalic acid',
    'Isobutanol': 'Isobutanol',
    'SodiumHP': None,  # salts: descriptors supplied manually if needed
    'Glucose': 'Glucose',
}

_MW_GLUCOSE = 180.156
_T_FERM = 310.15


def _f(x):
    try:
        x = float(x)
        return x if math.isfinite(x) else float('nan')
    except Exception:
        return float('nan')


def stoichiometric_descriptors(n_C, n_H, n_O, n_N, MW):
    """Degree of reduction and stoichiometric maximum yields on glucose."""
    nan = float('nan')
    if not n_C:
        return dict(dor=nan, Ymax_carbon=nan, Ymax_electron=nan, Ymax_stoich=nan)
    dor = (4 * n_C + n_H - 2 * n_O - 3 * n_N) / n_C
    Y_c = 6. / n_C * MW / _MW_GLUCOSE  # mol C conserved, no CO2 fixation
    Y_e = 24. / (dor * n_C) * MW / _MW_GLUCOSE if dor > 0 else nan
    return dict(dor=dor, Ymax_carbon=Y_c, Ymax_electron=Y_e,
                Ymax_stoich=min(Y_c, Y_e) if dor > 0 else Y_c)


def chemical_descriptors(chemical, T_ferm=_T_FERM):
    """Return the descriptor dict of a thermosteam Chemical (or search
    identifier). Requires thermosteam."""
    import numpy as np
    import thermosteam as tmo
    if isinstance(chemical, str):
        chemical = tmo.Chemical(chemical.replace(' ', '_'), search_ID=chemical)
    c = chemical
    atoms = {}
    try:
        atoms = dict(c.atoms)
    except Exception:
        pass
    d = dict(
        MW=_f(c.MW), Tb=_f(c.Tb), Tm=_f(c.Tm),
        Hvap_298=_f(_try(lambda: c.Hvap(298.15))),
        Psat_ferm=_f(_try(lambda: c.Psat(T_ferm))),
        LHV=_f(_try(lambda: c.LHV)),
        n_C=float(atoms.get('C', 0)), n_H=float(atoms.get('H', 0)),
        n_O=float(atoms.get('O', 0)), n_N=float(atoms.get('N', 0)),
        formula=_try(lambda: c.formula) or '',
    )
    d.update(stoichiometric_descriptors(d['n_C'], d['n_H'], d['n_O'], d['n_N'], d['MW']))
    # infinite-dilution activity coefficient & relative volatility vs water
    gamma = float('nan')
    try:
        water = tmo.Chemical('Water')
        chems = tmo.Chemicals([water, c.copy(c.ID + '_d') if c.ID == 'Water' else c])
        from thermosteam.equilibrium import (DortmundActivityCoefficients,
                                             UNIFACActivityCoefficients)
        x = np.array([1 - 1e-6, 1e-6])
        for cls in (DortmundActivityCoefficients, UNIFACActivityCoefficients):
            try:
                gamma = float(cls(chems)(x, T_ferm)[1])
                if math.isfinite(gamma) and gamma > 0:
                    break
            except Exception:
                gamma = float('nan')
        Pw = water.Psat(T_ferm)
        d['gamma_inf'] = gamma
        d['log10_alpha_inf'] = (math.log10(gamma * d['Psat_ferm'] / Pw)
                                if gamma > 0 and d['Psat_ferm'] > 0 else float('nan'))
    except Exception:
        d['gamma_inf'] = float('nan')
        d['log10_alpha_inf'] = float('nan')
    return d


def _try(f):
    try:
        return f()
    except Exception:
        return None


def build_registry(path=REGISTRY_PATH, extra=None):
    """Compute descriptors for KNOWN_CHEMICALS (+ extra {ID: search_ID}) and
    write the JSON registry."""
    out = {}
    items = dict(KNOWN_CHEMICALS)
    if extra:
        items.update(extra)
    for ID, sid in items.items():
        if sid is None:
            continue
        try:
            d = chemical_descriptors(sid)
            d['search_ID'] = sid
            out[ID] = d
        except Exception as e:  # pragma: no cover
            print(f'could not describe {ID} ({sid}): {e}')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=True,
                  default=lambda v: None)
    return out


def load_registry(path=REGISTRY_PATH):
    """Load the descriptor registry ({ID: descriptors})."""
    with open(path) as fh:
        reg = json.load(fh)
    for d in reg.values():
        for k, v in d.items():
            if v is None:
                d[k] = float('nan')
    return reg


if __name__ == '__main__':
    reg = build_registry()
    for k, v in reg.items():
        print(f"{k:18s} MW={v['MW']:.1f} dor={v['dor']:.2f} Ymax={v['Ymax_stoich']:.3f} "
              f"Tb={v['Tb']:.0f} gamma_inf={v['gamma_inf']:.2f} "
              f"log10a={v['log10_alpha_inf']:.2f}")
