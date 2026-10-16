#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Feature construction (DESIGN.md section 5): specifications of the feedstock,
fermentation product, final product and fermentation regime -> standardized
feature vector.

A :class:`FeatureSpec` is fit on training samples only (imputation values,
log-transform decisions, standardization statistics, category levels) and is
JSON-serializable, so that exactly the same transformation is applied at
inference time to a single user-specified dict.

Feature groups
--------------
* specifications: ``n_acid_groups``, ``pKa1``, ``final_product_stoich_mass_ratio``,
  ``T_fermentation``, ``theoretical_yield``, ``inflection_yield`` and the
  boolean flags ``neutralization``, ``upgrading``, ``aerobic``;
* chemical descriptors of the fermentation product (``fp_*``) and final
  product (``pp_*``), taken from the descriptor registry
  (``registry/chemicals.json``, see :mod:`descriptors`) keyed by
  ``ferm_product_chemical`` / ``final_product_chemical`` (falling back to the
  product labels), else from the sample's own ``fp_*``/``pp_*`` columns;
  user-supplied descriptor dicts take precedence (``descriptors=`` hook);
* feedstock composition: every ``fs_w_*`` mass fraction (wet basis);
* derived physics features (``phys_*``): ``Tb - 373.15``, log10 of the
  infinite-dilution activity coefficient, log10 relative volatility vs. water,
  log10 vapor-pressure ratio to water at 310.15 K, theoretical yield relative
  to the stoichiometric (and carbon) maxima, inflection/theoretical yield,
  mass-based LHV and heat of vaporization, MW ratio final/fermentation
  product, log10 solubility, log10 capacity, non-sugar solids fraction;
* canonical contextual parameters ``p_*`` (baseline-imputed + missing
  indicators);
* one-hot categoricals (``feedstock_type``, ``neutralizing_agent``,
  ``fermentation_mode``, ``separation_principle``) with an explicit
  ``__other__`` level for unseen values;
* productivity enters the network through ``log r`` and ``1/r``, which the
  model computes itself from ``log_r`` (so that ``dMPSP/dr`` is available by
  autograd); their standardization statistics are stored here.

Never used as features: ids (``sample_id``, ``config_id``, ``biorefinery``,
``sample_kind``, ``seed``), biorefinery-specific parameters ``bp_*``, product
and feedstock names, the TY domain bounds, and all targets. Prices used by the
price-explicit channels (``feedstock_price``, ``electricity_price``,
``natural_gas_price``) are returned separately in natural units and are only
fed to the encoder when ``price_explicit=False``.
"""
import json
import math
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

__all__ = (
    'FeatureSpec', 'FeatureBatch', 'MissingSpecError', 'ApplicabilityDomain',
    'ApplicabilityWarning', 'resolve_chemical_descriptors', 'derive_features',
    'build_spec_registry', 'apply_spec_registry', 'register_chemical',
    'describe_new_chemical', 'psat_water', 'PRICE_COLUMNS', 'PRICE_KEYS',
    'REQUIRED_FIELDS',
)

try:
    from .physics import PRICE_COLUMNS as _PC, INTENSITY_KEYS as _IK
    PRICE_KEYS = tuple(_IK)
    PRICE_COLUMNS = tuple(_PC[k] for k in PRICE_KEYS)
except Exception:  # pragma: no cover
    PRICE_KEYS = ('fs', 'el', 'ng')
    PRICE_COLUMNS = ('feedstock_price', 'electricity_price', 'natural_gas_price')

try:
    from .data_generation.schema import CANONICAL_PARAMETERS
except Exception:  # pragma: no cover
    CANONICAL_PARAMETERS = (
        'p_feedstock_price', 'p_feedstock_capacity', 'p_electricity_price',
        'p_natural_gas_price', 'p_operating_days', 'p_IRR', 'p_income_tax',
        'p_boiler_efficiency', 'p_turbogenerator_efficiency',
        'p_cell_mass_yield', 'p_base_price', 'p_acid_price',
        'p_enzyme_loading', 'p_saccharification_yield')

#: Aliases from which a missing price column is filled.
PRICE_ALIASES = {
    'feedstock_price': ('p_feedstock_price', 'fs_price'),
    'electricity_price': ('p_electricity_price',),
    'natural_gas_price': ('p_natural_gas_price',),
}
#: Columns that duplicate prices or capacity (never used directly).
EXCLUDED_ALWAYS = ('p_feedstock_price', 'p_electricity_price', 'p_natural_gas_price',
                   'fs_price', 'fs_F_mass', 'p_feedstock_capacity')

W_SUGAR = 'fs_w_fermentable_sugar_eq'
G_MAX = 'final_product_stoich_mass_ratio'

#: Fields that cannot be imputed (prediction raises MissingSpecError).
REQUIRED_FIELDS = ('r', 'theoretical_yield', 'inflection_yield', G_MAX, W_SUGAR,
                   'feedstock_price')

SPEC_NUMERIC = ('n_acid_groups', 'pKa1', G_MAX, 'T_fermentation',
                'theoretical_yield', 'inflection_yield')
SPEC_BOOLEAN = ('neutralization', 'upgrading', 'aerobic')
SPEC_CATEGORICAL = ('feedstock_type', 'neutralizing_agent', 'fermentation_mode',
                    'separation_principle')

#: Chemical descriptors (registry keys) resolved for fp_ and pp_.
CHEM_DESCRIPTORS = ('MW', 'Tb', 'Tm', 'Hvap_298', 'Psat_ferm', 'LHV', 'n_C', 'n_H',
                    'n_O', 'n_N', 'dor', 'Ymax_carbon', 'Ymax_electron',
                    'Ymax_stoich', 'gamma_inf', 'log10_alpha_inf')
#: Raw descriptors used as features directly (others enter via phys_*).
CHEM_FEATURES = ('n_C', 'n_H', 'n_O', 'n_N', 'dor', 'Ymax_carbon', 'Ymax_electron',
                 'Ymax_stoich', 'Tm', 'log10_alpha_inf')
CHEM_KEYS = {'fp': ('ferm_product_chemical', 'fermentation_product'),
             'pp': ('final_product_chemical', 'final_product')}

OTHER = '__other__'
_T_REF = 310.15
_MW_GLUCOSE = 180.156

#: Session-level registry additions (see :func:`register_chemical`).
_USER_CHEMICALS = {}


class MissingSpecError(ValueError):
    """Raised when a required specification field is missing."""


class ApplicabilityWarning(UserWarning):
    """Issued when a query lies outside the training applicability domain."""


# %% Descriptor registry

def _package_registry():
    try:
        from .descriptors import load_registry
        return load_registry()
    except Exception:
        return {}


def register_chemical(ID, descriptors):
    """Register descriptors of a new chemical for this Python session (used by
    every FeatureSpec.transform call). ``descriptors`` follows the keys of
    ``registry/chemicals.json``; missing ``dor``/``Ymax_*`` are computed from
    the atom counts and MW."""
    _USER_CHEMICALS[str(ID)] = _complete_descriptors(dict(descriptors))


def describe_new_chemical(ID, search_ID=None, register=True):
    """Compute descriptors of a new chemical with thermosteam (if installed)
    via ``descriptors.chemical_descriptors`` and optionally register them."""
    from .descriptors import chemical_descriptors
    d = chemical_descriptors(search_ID or ID)
    if register:
        register_chemical(ID, d)
    return d


def _num(v):
    try:
        v = float(v)
        return v if math.isfinite(v) else float('nan')
    except Exception:
        return float('nan')


def _complete_descriptors(d):
    """Fill dor / Ymax_* from atom counts if absent."""
    try:
        from .descriptors import stoichiometric_descriptors
    except Exception:  # pragma: no cover
        return d
    need = [k for k in ('dor', 'Ymax_carbon', 'Ymax_electron', 'Ymax_stoich')
            if not np.isfinite(_num(d.get(k)))]
    atoms = [_num(d.get(k, 0.)) for k in ('n_C', 'n_H', 'n_O', 'n_N')]
    MW = _num(d.get('MW'))
    if need and np.isfinite(MW) and np.isfinite(atoms[0]) and atoms[0] > 0:
        atoms = [0. if not np.isfinite(a) else a for a in atoms]
        s = stoichiometric_descriptors(*atoms, MW)
        for k in need:
            d[k] = s[k]
    return d


def resolve_chemical_descriptors(df, registry=None, descriptors=None):
    """Return a copy of ``df`` whose ``fp_<k>``/``pp_<k>`` columns (k in
    CHEM_DESCRIPTORS) are taken from (in order of precedence) the
    ``descriptors`` hook ({chemical_ID: dict}), session registrations
    (:func:`register_chemical`), ``registry`` and the rows' own columns."""
    df = df.copy()
    reg = {}
    reg.update(registry if registry is not None else _package_registry())
    reg.update(_USER_CHEMICALS)
    if descriptors:
        reg.update({k: _complete_descriptors(dict(v)) for k, v in descriptors.items()})
    for pre, keys in CHEM_KEYS.items():
        key = pd.Series([None] * len(df), index=df.index, dtype=object)
        for kc in keys:
            if kc in df:
                v = df[kc].where(df[kc].notna() & (df[kc].astype(str) != ''))
                key = key.where(key.notna(), v)
        cols = {k: (pd.to_numeric(df[f'{pre}_{k}'], errors='coerce').to_numpy(float).copy()
                    if f'{pre}_{k}' in df else np.full(len(df), np.nan))
                for k in CHEM_DESCRIPTORS}
        keyv = key.to_numpy(object)
        for ID in pd.unique(key.dropna()):
            d = reg.get(str(ID))
            if d is None:
                continue
            m = keyv == ID
            for k in CHEM_DESCRIPTORS:
                v = _num(d.get(k))
                if np.isfinite(v):
                    cols[k][m] = v
        # complete stoichiometric descriptors from atoms where missing
        try:
            from .descriptors import stoichiometric_descriptors
            miss = ~np.isfinite(cols['dor']) & np.isfinite(cols['MW']) & (np.nan_to_num(cols['n_C']) > 0)
            for i in np.flatnonzero(miss):
                s = stoichiometric_descriptors(*(np.nan_to_num(cols[a][i]) for a in ('n_C', 'n_H', 'n_O', 'n_N')),
                                               cols['MW'][i])
                for k, v in s.items():
                    if not np.isfinite(cols[k][i]):
                        cols[k][i] = v
        except Exception:  # pragma: no cover
            pass
        for k, v in cols.items():
            df[f'{pre}_{k}'] = v
    return df


def used_registry_entries(df, registry=None):
    """Registry entries referenced by ``df`` (stored with the FeatureSpec)."""
    reg = registry if registry is not None else _package_registry()
    out = {}
    for keys in CHEM_KEYS.values():
        for kc in keys:
            if kc in df:
                for ID in pd.unique(df[kc].dropna().astype(str)):
                    if ID in reg:
                        out[ID] = {k: (None if isinstance(v, float) and not math.isfinite(v) else v)
                                   for k, v in reg[ID].items()}
    return out


# %% Derived physics features

def psat_water(T):
    """Water vapor pressure [Pa] (Antoine, 1-100 C)."""
    TC = np.asarray(T, float) - 273.15
    return 133.322 * 10 ** (8.07131 - 1730.63 / (233.426 + TC))


def _col(df, name):
    if name in df:
        return pd.to_numeric(df[name], errors='coerce').astype(float)
    return pd.Series(np.nan, index=df.index, dtype=float)


def _log10_pos(x):
    x = np.asarray(x, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(x > 0, np.log10(np.where(x > 0, x, 1.)), np.nan)


def _ratio(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(np.abs(b) > 0, a / np.where(np.abs(b) > 0, b, 1.), np.nan)


def derive_features(df):
    """Derived physics features (``phys_*``) from resolved descriptors."""
    out = {}
    Pw = float(psat_water(_T_REF))
    for p in ('fp', 'pp'):
        MW = _col(df, f'{p}_MW')
        out[f'phys_{p}_log10_MW'] = _log10_pos(MW)
        out[f'phys_{p}_Tb_minus_water'] = _col(df, f'{p}_Tb').to_numpy() - 373.15
        out[f'phys_{p}_log10_gamma_inf'] = _log10_pos(_col(df, f'{p}_gamma_inf'))
        out[f'phys_{p}_log10_Psat_ratio'] = _log10_pos(_col(df, f'{p}_Psat_ferm') / Pw)
        out[f'phys_{p}_LHV_mass'] = _ratio(_col(df, f'{p}_LHV'), MW) / 1e3  # MJ/kg
        out[f'phys_{p}_Hvap_mass'] = _ratio(_col(df, f'{p}_Hvap_298'), MW) / 1e3  # kJ/g
    ty = _col(df, 'theoretical_yield')
    out['phys_ty_over_Ymax_stoich'] = _ratio(ty, _col(df, 'fp_Ymax_stoich'))
    out['phys_ty_over_Ymax_carbon'] = _ratio(ty, _col(df, 'fp_Ymax_carbon'))
    out['phys_inflection_over_theoretical'] = _ratio(_col(df, 'inflection_yield'), ty)
    out['phys_pp_fp_MW_ratio'] = _ratio(_col(df, 'pp_MW'), _col(df, 'fp_MW'))
    out['phys_log10_solubility'] = _log10_pos(_col(df, 'product_water_solubility'))
    cap = _col(df, 'p_feedstock_capacity')
    cap = cap.where(cap.notna(), _col(df, 'fs_F_mass'))
    out['phys_log10_capacity'] = _log10_pos(cap)
    out['phys_fs_nonsugar_solids'] = (1. - _col(df, W_SUGAR) - _col(df, 'fs_w_water').fillna(0.)).to_numpy()
    return pd.DataFrame(out, index=df.index)


# %% Spec-level registry (fill omitted fields by name)

_FS_REG_COLS = ('feedstock_type',)
_FP_REG_COLS = ('ferm_product_chemical', 'theoretical_yield', 'inflection_yield',
                'n_acid_groups', 'pKa1', 'product_water_solubility')
_PAIR_REG_COLS = ('final_product_chemical', G_MAX, 'upgrading')


def _mode_or_median(s):
    s = s.dropna()
    if not len(s):
        return None
    if pd.api.types.is_bool_dtype(s):
        return bool(s.mode().iloc[0])
    if pd.api.types.is_numeric_dtype(s):
        return float(np.median(s.astype(float)))
    return str(s.astype(str).mode().iloc[0])


def build_spec_registry(samples):
    """Per-name defaults learned from training samples so that users can
    specify known feedstocks/products by name: {'feedstock': {name: {...}},
    'fermentation_product': {...}, 'pair': {'fp|pp': {...}}}."""
    reg = {'feedstock': {}, 'fermentation_product': {}, 'pair': {}}
    if 'feedstock' in samples:
        cols = [c for c in samples if c.startswith('fs_w_')] + list(_FS_REG_COLS)
        for name, grp in samples.groupby('feedstock'):
            reg['feedstock'][str(name)] = {c: _mode_or_median(grp[c]) for c in cols if c in grp}
    if 'fermentation_product' in samples:
        for name, grp in samples.groupby('fermentation_product'):
            reg['fermentation_product'][str(name)] = {c: _mode_or_median(grp[c]) for c in _FP_REG_COLS if c in grp}
    if 'fermentation_product' in samples and 'final_product' in samples:
        for (a, b), grp in samples.groupby(['fermentation_product', 'final_product']):
            reg['pair'][f'{a}|{b}'] = {c: _mode_or_median(grp[c]) for c in _PAIR_REG_COLS if c in grp}
    return reg


def _missing(v):
    return v is None or (isinstance(v, float) and not math.isfinite(v)) or (isinstance(v, str) and v == '')


def apply_spec_registry(spec, reg):
    """Return a copy of the dict ``spec`` with omitted fields filled from the
    spec registry (user values always take precedence)."""
    spec = dict(spec)
    filled = []
    entries = []
    if spec.get('feedstock') is not None:
        entries.append(reg.get('feedstock', {}).get(str(spec['feedstock']), {}))
    if spec.get('fermentation_product') is not None:
        entries.append(reg.get('fermentation_product', {}).get(str(spec['fermentation_product']), {}))
        if spec.get('final_product') is not None:
            entries.append(reg.get('pair', {}).get(f"{spec['fermentation_product']}|{spec['final_product']}", {}))
    for e in entries:
        for k, v in e.items():
            if _missing(spec.get(k)) and not _missing(v):
                spec[k] = v
                filled.append(k)
    return spec, filled


# %% Feature spec

def _to_float(s):
    """Coerce booleans/strings/numbers to float (NaN if not convertible)."""
    if pd.api.types.is_bool_dtype(s):
        return s.astype(float)
    if s.dtype == object:
        m = {'true': 1., 'false': 0., 'yes': 1., 'no': 0.}
        s = s.map(lambda v: m.get(v.strip().lower(), v) if isinstance(v, str) else v)
        s = s.map(lambda v: float(v) if isinstance(v, (bool, np.bool_)) else v)
    return pd.to_numeric(s, errors='coerce').astype(float)


def _cat(s):
    def f(v):
        if v is None or (isinstance(v, float) and not math.isfinite(v)):
            return 'none'
        v = str(v).strip()
        return v if v and v.lower() != 'nan' else 'none'
    return s.map(f)


@dataclass
class FeatureBatch:
    """Transformed inputs. ``X`` excludes productivity (the model computes the
    r-features from ``log_r``); ``prices`` are in natural units, columns in the
    order of ``PRICE_COLUMNS``."""
    X: np.ndarray
    log_r: np.ndarray
    prices: np.ndarray
    w_sugar: np.ndarray
    g_max: np.ndarray
    index: object = None
    imputed: list = field(default_factory=list)

    def __len__(self):
        return len(self.log_r)

    def subset(self, idx):
        idx = np.asarray(idx)
        return FeatureBatch(self.X[idx], self.log_r[idx], self.prices[idx],
                            self.w_sugar[idx], self.g_max[idx],
                            None if self.index is None else np.asarray(self.index)[idx],
                            list(self.imputed))


class FeatureSpec:
    """Feature transformation fit on training samples only.

    Parameters
    ----------
    price_explicit : bool
        If True (physics model), the prices used by price-explicit channels
        are excluded from the encoder features so that predictions are exactly
        linear in them. If False, they are ordinary (standardized) features.
    use_names : bool
        One-hot encode feedstock / product names too (off by default: forces
        generalization through descriptors).
    clip : float
        Standardized continuous features are clipped to [-clip, clip].
    exclude : iterable of str
        Additional columns never to use.
    """

    VERSION = 1

    def __init__(self, price_explicit=True, use_names=False, clip=10., exclude=(),
                 log_ratio=100., min_level_count=1):
        self.price_explicit = bool(price_explicit)
        self.use_names = bool(use_names)
        self.clip = float(clip)
        self.exclude = tuple(exclude)
        self.log_ratio = float(log_ratio)
        self.min_level_count = int(min_level_count)
        self.numeric = []        # list of dicts: name, log, fill, mean, std, kind
        self.indicators = []     # numeric names with a missing indicator
        self.categoricals = {}   # name -> levels (incl. OTHER)
        self.r_stats = {}
        self.price_stats = {}
        self.price_fill = {}
        self.registry = {}       # chemical descriptor entries used in training
        self.spec_registry = {}  # per-name defaults
        self.fitted = False

    # ---- helpers ----
    def _prepare(self, df, descriptors=None):
        df = resolve_chemical_descriptors(df, self._chem_registry(), descriptors)
        der = derive_features(df)
        df = pd.concat([df.drop(columns=[c for c in der if c in df]), der], axis=1)
        for c, aliases in PRICE_ALIASES.items():
            v = _col(df, c)
            for a in aliases:
                v = v.where(v.notna(), _col(df, a))
            df[c] = v
        return df

    def _chem_registry(self):
        reg = _package_registry()
        reg.update({k: {kk: (float('nan') if vv is None else vv) for kk, vv in v.items()}
                    for k, v in self.registry.items()})
        return reg

    def _numeric_candidates(self, df):
        cands = []
        for c in SPEC_NUMERIC:
            cands.append((c, 'spec'))
        for c in SPEC_BOOLEAN:
            cands.append((c, 'flag'))
        for p in ('fp', 'pp'):
            for k in CHEM_FEATURES:
                cands.append((f'{p}_{k}', 'descriptor'))
        for c in df.columns:
            if c.startswith('fs_w_'):
                cands.append((c, 'feedstock'))
        for c in df.columns:
            if c.startswith('phys_'):
                # capacity is contextual (baseline-imputed like p_* parameters)
                cands.append((c, 'param' if c == 'phys_log10_capacity' else 'derived'))
        params = list(CANONICAL_PARAMETERS) + sorted(
            c for c in df.columns if c.startswith('p_') and c not in CANONICAL_PARAMETERS)
        for c in params:
            cands.append((c, 'param'))
        if not self.price_explicit:
            for c in PRICE_COLUMNS:
                cands.append((c, 'price'))
        seen, out = set(), []
        for c, kind in cands:
            if c in seen or c in EXCLUDED_ALWAYS or c in self.exclude or c.startswith('bp_'):
                continue
            seen.add(c)
            out.append((c, kind))
        return out

    # ---- fit ----
    def fit(self, samples, descriptors=None):
        """Fit on training samples (DataFrame in the ``samples`` format)."""
        self.registry = used_registry_entries(samples, self._chem_registry())
        if descriptors:
            self.registry.update({k: dict(v) for k, v in descriptors.items()})
        df = self._prepare(samples, descriptors)
        n = len(df)
        base = (df['sample_kind'].astype(str) == 'baseline').to_numpy() if 'sample_kind' in df else np.zeros(n, bool)
        self.numeric, self.indicators = [], []
        for c, kind in self._numeric_candidates(df):
            if c not in df:
                continue
            v = _to_float(df[c]).to_numpy()
            fin = np.isfinite(v)
            if not fin.any():
                continue
            use_log = False
            if kind not in ('flag',) and not c.startswith('phys_') and (v[fin] > 0).all():
                q05, q95 = np.percentile(v[fin], [5, 95])
                use_log = q05 > 0 and q95 / q05 > self.log_ratio
            w = np.log10(v) if use_log else v.copy()
            if kind == 'param' and (fin & base).any():
                fill = float(np.median(w[fin & base]))
            else:
                fill = float(np.median(w[fin]))
            w = np.where(fin, w, fill)
            mean, std = float(w.mean()), float(w.std())
            has_missing = not fin.all()
            if std < 1e-12 * max(1., abs(mean)) and not has_missing:
                continue  # constant in training
            self.numeric.append(dict(name=c, kind=kind, log=bool(use_log), fill=fill,
                                     mean=mean, std=std if std > 1e-12 else 1.))
            if has_missing:
                self.indicators.append(c)
        cats = list(SPEC_CATEGORICAL)
        if self.use_names:
            cats += ['feedstock', 'fermentation_product', 'final_product']
        self.categoricals = {}
        for c in cats:
            if c in self.exclude:
                continue
            s = _cat(df[c]) if c in df else pd.Series('none', index=df.index)
            counts = s.value_counts()
            levels = sorted(str(k) for k, v in counts.items() if v >= self.min_level_count)
            self.categoricals[c] = levels + [OTHER]
        r = _to_float(df['r']).to_numpy()
        lr = np.log(r[np.isfinite(r) & (r > 0)])
        ir = np.exp(-lr)
        self.r_stats = dict(log_r_mean=float(lr.mean()), log_r_std=float(max(lr.std(), 1e-3)),
                            inv_r_mean=float(ir.mean()), inv_r_std=float(max(ir.std(), 1e-3)),
                            log_r_min=float(lr.min()), log_r_max=float(lr.max()))
        self.price_stats, self.price_fill = {}, {}
        for c in PRICE_COLUMNS:
            v = _col(df, c).to_numpy()
            fin = np.isfinite(v)
            fill = float(np.median(v[fin & base])) if (fin & base).any() else (
                float(np.median(v[fin])) if fin.any() else 0.)
            vv = np.where(fin, v, fill)
            self.price_fill[c] = fill
            self.price_stats[c] = dict(mean=float(vv.mean()), std=float(max(vv.std(), 1e-12)),
                                       min=float(vv.min()), max=float(vv.max()))
        self.spec_registry = build_spec_registry(samples)
        self.fitted = True
        return self

    # ---- columns ----
    @property
    def columns(self):
        cols = [d['name'] for d in self.numeric]
        cols += [f'missing__{c}' for c in self.indicators]
        for c, levels in self.categoricals.items():
            cols += [f'{c}=={lv}' for lv in levels]
        return cols

    @property
    def n_features(self):
        return len(self.columns)

    def continuous_mask(self):
        """Boolean mask of standardized continuous columns (for noise)."""
        m = np.zeros(self.n_features, bool)
        m[:len(self.numeric)] = True
        return m

    def categorical_blocks(self):
        """[(start, stop, other_offset)] column blocks of each categorical."""
        out, i = [], len(self.numeric) + len(self.indicators)
        for c, levels in self.categoricals.items():
            out.append((i, i + len(levels), len(levels) - 1))
            i += len(levels)
        return out

    # ---- transform ----
    def transform(self, data, descriptors=None, strict=True):
        """Transform samples (DataFrame, dict or list of dicts) into a
        :class:`FeatureBatch`. With ``strict``, missing REQUIRED_FIELDS raise
        :class:`MissingSpecError`; otherwise they become NaN."""
        if not self.fitted:
            raise RuntimeError('FeatureSpec is not fitted')
        df = _as_frame(data)
        df = self._prepare(df, descriptors)
        n = len(df)
        missing_req = [c for c in REQUIRED_FIELDS
                       if c not in df or not np.isfinite(_to_float(df[c]).to_numpy()).all()]
        if missing_req and strict:
            raise MissingSpecError(
                'missing required specification field(s): ' + ', '.join(missing_req)
                + '. Required: ' + ', '.join(REQUIRED_FIELDS)
                + " (feedstock_price may also be given as 'p_feedstock_price' or 'fs_price').")
        imputed = []
        Xn = np.empty((n, len(self.numeric)))
        miss = {}
        for j, d in enumerate(self.numeric):
            c = d['name']
            v = _to_float(df[c]).to_numpy() if c in df else np.full(n, np.nan)
            if d['log']:
                v = _log10_pos(v)
            fin = np.isfinite(v)
            if not fin.all():
                miss[c] = ~fin
                if d['kind'] != 'param' and c not in self.indicators:
                    imputed.append(c)  # never missing in training
            v = np.where(fin, v, d['fill'])
            Xn[:, j] = np.clip((v - d['mean']) / d['std'], -self.clip, self.clip)
        Xi = np.column_stack([miss.get(c, np.zeros(n, bool)).astype(float) for c in self.indicators]) \
            if self.indicators else np.zeros((n, 0))
        blocks = []
        for c, levels in self.categoricals.items():
            s = _cat(df[c]) if c in df else pd.Series('none', index=df.index)
            idx = {lv: i for i, lv in enumerate(levels)}
            k = s.map(lambda v: idx.get(v, len(levels) - 1)).to_numpy(int)
            B = np.zeros((n, len(levels)))
            B[np.arange(n), k] = 1.
            blocks.append(B)
        X = np.hstack([Xn, Xi] + blocks)
        r = _to_float(df['r']).to_numpy() if 'r' in df else np.full(n, np.nan)
        with np.errstate(divide='ignore', invalid='ignore'):
            log_r = np.where(r > 0, np.log(np.where(r > 0, r, 1.)), np.nan)
        P = np.empty((n, len(PRICE_COLUMNS)))
        for j, c in enumerate(PRICE_COLUMNS):
            v = _col(df, c).to_numpy()
            if not np.isfinite(v).all() and c != 'feedstock_price':
                imputed.append(c)
            P[:, j] = np.where(np.isfinite(v), v, self.price_fill[c])
        w = _col(df, W_SUGAR).to_numpy()
        gm = _col(df, G_MAX).to_numpy()
        return FeatureBatch(X, log_r, P, w, gm, df.index.to_numpy(), sorted(set(imputed)))

    def r_features(self, log_r):
        """Standardized (log r, 1/r) features (numpy; the model has a torch twin)."""
        s = self.r_stats
        log_r = np.asarray(log_r, float)
        return np.column_stack([(log_r - s['log_r_mean']) / s['log_r_std'],
                                (np.exp(-log_r) - s['inv_r_mean']) / s['inv_r_std']])

    def full_matrix(self, batch):
        """Encoder features + r-features + standardized prices (when not
        already in ``X``): the design matrix for baselines and the
        applicability-domain check."""
        parts = [batch.X, self.r_features(batch.log_r)]
        if self.price_explicit:
            st = [self.price_stats[c] for c in PRICE_COLUMNS]
            parts.append((batch.prices - np.array([s['mean'] for s in st]))
                         / np.array([s['std'] for s in st]))
        return np.hstack(parts)

    def full_columns(self):
        cols = self.columns + ['log_r', 'inv_r']
        if self.price_explicit:
            cols += list(PRICE_COLUMNS)
        return cols

    # ---- serialization ----
    def to_dict(self):
        return dict(version=self.VERSION, price_explicit=self.price_explicit,
                    use_names=self.use_names, clip=self.clip, exclude=list(self.exclude),
                    log_ratio=self.log_ratio, min_level_count=self.min_level_count,
                    numeric=self.numeric, indicators=self.indicators,
                    categoricals=self.categoricals, r_stats=self.r_stats,
                    price_stats=self.price_stats, price_fill=self.price_fill,
                    registry=self.registry, spec_registry=self.spec_registry,
                    columns=self.columns)

    @classmethod
    def from_dict(cls, d):
        self = cls(price_explicit=d['price_explicit'], use_names=d['use_names'], clip=d['clip'],
                   exclude=d.get('exclude', ()), log_ratio=d.get('log_ratio', 100.),
                   min_level_count=d.get('min_level_count', 1))
        self.numeric = [dict(x) for x in d['numeric']]
        self.indicators = list(d['indicators'])
        self.categoricals = {k: list(v) for k, v in d['categoricals'].items()}
        self.r_stats = dict(d['r_stats'])
        self.price_stats = {k: dict(v) for k, v in d['price_stats'].items()}
        self.price_fill = dict(d['price_fill'])
        self.registry = d.get('registry', {})
        self.spec_registry = d.get('spec_registry', {})
        self.fitted = True
        if d.get('columns') is not None and list(d['columns']) != self.columns:
            raise ValueError('inconsistent FeatureSpec columns')
        return self

    def to_json(self, path=None):
        s = json.dumps(_jsonable(self.to_dict()), indent=1)
        if path is not None:
            with open(path, 'w') as fh:
                fh.write(s)
        return s

    @classmethod
    def from_json(cls, s_or_path):
        s = s_or_path
        if not s.lstrip().startswith('{'):
            with open(s_or_path) as fh:
                s = fh.read()
        return cls.from_dict(json.loads(s))


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, float)):
        o = float(o)
        return o if math.isfinite(o) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def _as_frame(data):
    if isinstance(data, pd.DataFrame):
        return data.reset_index(drop=True) if not data.index.is_unique else data
    if isinstance(data, pd.Series):
        return data.to_frame().T
    if isinstance(data, dict):
        return pd.DataFrame([data])
    return pd.DataFrame(list(data))


# %% Applicability domain

class ApplicabilityDomain:
    """kNN distance in the standardized full feature space. A query is
    flagged when its mean distance to the ``k`` nearest training rows exceeds
    the ``quantile`` of the training rows' leave-one-out kNN distances, or
    when a continuous feature lies outside the training range."""

    def __init__(self, k=5, quantile=0.99, max_rows=5000, seed=0):
        self.k, self.quantile, self.max_rows, self.seed = int(k), float(quantile), int(max_rows), seed
        self.Z = None
        self.threshold = None
        self.lo = self.hi = None
        self.columns = None

    def fit(self, Z, columns=None):
        Z = np.asarray(Z, float)
        rng = np.random.RandomState(self.seed)
        if len(Z) > self.max_rows:
            Z = Z[rng.choice(len(Z), self.max_rows, replace=False)]
        self.Z = Z
        self.lo, self.hi = Z.min(0), Z.max(0)
        self.columns = list(columns) if columns is not None else None
        d = self._knn(Z, exclude_self=True)
        self.threshold = float(np.quantile(d, self.quantile)) if len(d) else float('inf')
        return self

    def _knn(self, Q, exclude_self=False):
        k = min(self.k, len(self.Z) - int(exclude_self))
        if k <= 0:
            return np.zeros(len(Q))
        out = np.empty(len(Q))
        zz = (self.Z ** 2).sum(1)
        for s in range(0, len(Q), 1024):
            q = Q[s:s + 1024]
            D = np.sqrt(np.maximum((q ** 2).sum(1)[:, None] + zz[None] - 2 * q @ self.Z.T, 0.))
            if exclude_self:
                D[np.arange(len(q)), np.arange(s, s + len(q))] = np.inf
            out[s:s + 1024] = np.sort(D, 1)[:, :k].mean(1)
        return out

    def score(self, Q):
        """Return (distance / threshold, list of out-of-range column names per row)."""
        Q = np.asarray(Q, float)
        ratio = self._knn(Q) / max(self.threshold, 1e-12)
        tol = 1e-6 + 1e-3 * (self.hi - self.lo)
        oor = (Q < self.lo - tol) | (Q > self.hi + tol)
        names = self.columns or [str(i) for i in range(Q.shape[1])]
        return ratio, [[names[j] for j in np.flatnonzero(row)] for row in oor]

    def to_dict(self):
        return dict(k=self.k, quantile=self.quantile, max_rows=self.max_rows, seed=self.seed,
                    threshold=self.threshold, lo=self.lo.tolist(), hi=self.hi.tolist(),
                    columns=self.columns)

    @classmethod
    def from_dict(cls, d, Z):
        self = cls(d['k'], d['quantile'], d['max_rows'], d['seed'])
        self.Z = np.asarray(Z, float)
        self.threshold = d['threshold']
        self.lo, self.hi = np.asarray(d['lo']), np.asarray(d['hi'])
        self.columns = d.get('columns')
        return self


def warn_if_outside(ratio, oor, index=None):
    """Issue an ApplicabilityWarning for rows outside the domain."""
    bad = [i for i, (r, o) in enumerate(zip(ratio, oor)) if r > 1 or o]
    if bad:
        parts = []
        for i in bad[:5]:
            lab = index[i] if index is not None else i
            s = f'row {lab}: kNN distance {ratio[i]:.2f}x training threshold'
            if oor[i]:
                s += '; outside training range: ' + ', '.join(oor[i][:6])
            parts.append(s)
        warnings.warn('query outside the applicability domain of the training data '
                      '(extrapolation; predictions may be unreliable): ' + ' | '.join(parts),
                      ApplicabilityWarning, stacklevel=3)
    return bad
