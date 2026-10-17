# -*- coding: utf-8 -*-
"""
Structured final-product and downstream-process specifications
(DESIGN.md section 5, "Downstream specification").

The adapters record the fermentation regime in detail but the downstream
process only as free text (``separation_principle``, which is the product
identity in disguise) and, for a few configurations, ``product_purity_spec`` /
``product_form`` / ``coproducts``. The fields below describe the final
product and the separation route with controlled vocabularies, so that a
model can relate costs to *how* the product is recovered and a user can
describe the downstream process of a new product.

Fields (units / vocabularies in :data:`FIELDS`)
------------------------------------------------
``final_product_mass_fraction``  w_f [kg product / kg final-product stream]
    at the configuration's baseline (purity of the sold stream; the MPSP is
    per kg *pure* product in every adapter). ``None`` for a broth product,
    whose w_f follows the titer (see ``final_product_is_broth``).
``final_product_is_broth``  bool: the final product is the (clarified)
    fermentation broth itself, w_f(t) ~ t/1000 (oxalic acid: 0.005 at 5 g/L,
    0.031 at 30 g/L, 0.117 at 100 g/L in the recorded data).
``final_product_phase``  one of :data:`PHASES`.
``water_removal_route``  one of :data:`ROUTES` -- how the broth water is
    separated from the product: ``evaporate_water`` (non-volatile product,
    water boiled off: multi-effect evaporation / flash), ``crystallize``
    (product crystallized from the evaporated liquor), ``distill_product``
    (product more volatile than water, stripped from the broth), ``none``
    (the broth is sold; no water removal).
``separation_ops``  '+'-joined subset of :data:`SEPARATION_OPS` (separation
    and purification unit operations of the separation and upgrading areas,
    including coproduct trains such as DDGS drying).
``upgrading_ops``  '+'-joined subset of :data:`UPGRADING_OPS` (chemical
    conversions of the fermentation product after fermentation).
``coproducts``  '+'-joined subset of :data:`COPRODUCTS` (products sold besides
    the main product; ``electricity`` = on-site boiler-turbogenerator
    generating electricity from residues/biogas/natural gas, surplus sold).
``neutralization_degree``  fraction of the product's acid groups converted to
    the salt in fermentation [mol base-eq. / mol acid group] (0 = low-pH or
    not an acid; 1 = full neutralization).
``acid_regeneration``  acid used to regenerate the free acid from its salt
    downstream (``'H2SO4'``) or ``'none'``.

Empty multi-valued fields and no acid regeneration are written ``'none'``
(explicit; an omitted field -- None/NaN -- means "not specified" and is
imputed and reported by the feature transformation).

The table :data:`TABLE` holds these fields for every planned configuration
(``experiments.PLANNED_CONFIGS``); ``data_patches`` applies it to config
records (patch ``S1``, versioned by ``downstream_spec_version``) and
:func:`apply_to_frame` to already-processed tables.

Sources
-------
* adapters (``data_generation/adapters``, read-only): ``separation_principle``,
  ``separation_variant``, ``neutralization``/``neutralizing_agent``,
  ``product_purity_spec`` (lactic 0.88), ``product_form`` (oxalic
  'clarified_broth'), ``coproducts`` (isobutanol), ``neutralization_degree``
  (oxalic 0.91), partial pH-control base (succinic 0.42 mol OH-eq / mol);
* the Bioindustrial-Park systems (worktrees ``/home/user/work/bip_<name>``,
  read-only): unit operations of the separation/upgrading areas
  (``HP/process_areas.py``: ``create_HP_separation_improved_process``,
  ``..._hexanol_extraction_process``, ``create_HP_to_acrylic_acid_upgrading_process``;
  ``TAL/process_areas.py``: ``create_TAL_separation_solubility_exploit_process``,
  ``create_TAL_to_sorbic_acid_upgrading_process``; ``succinic/system_*.py``;
  ``lactic/_processes.py``; ``isobutanol/separations.py``; the corn,
  cellulosic and cane ethanol modules) and the boiler-turbogenerator setup;
* the recorded data (baseline samples, per-point purity ``x_*purity*`` of the
  final product stream and ``electricity_production``): w_f is the baseline
  median, ``electricity`` is listed where the CHP generates electricity at
  baseline (not lactic: its boiler covers steam only; not corn ethanol: no
  turbogenerator).

Notes: TAL's w_f varies with titer (baseline samples 0.60-0.93; the value is
the baseline median per feedstock).

succinic_neutral (DESIGN.md section 1, "Known data defects", D15): the table
describes the neutral design (salt formation, acidulation with H2SO4, gypsum
filtration: neutralization_degree 1, acid_regeneration 'H2SO4'). The data of
the v2 models did NOT contain it (succinic adapter v2 left the module default
``base_neutralizes_product=False``: the purchased lime was discarded and
R401 / R401_P / S405 bypassed, so acid, gypsum and acidulation costs were 0).
Adapter v3 simulates the salt formation and the acidulation; the four
configurations were regenerated (v2 rows are rejected, ``unpatchable_D15``),
so the table and the data now agree (product purity at baseline unchanged,
0.8486). Models trained before the regeneration (v2) still carry the mismatch.
"""
from __future__ import annotations

import math
import re

__all__ = ('SPEC_VERSION', 'FIELDS', 'STRING_FIELDS', 'MULTI_FIELDS', 'PHASES', 'ROUTES',
           'SEPARATION_OPS', 'UPGRADING_OPS', 'COPRODUCTS', 'ACIDS', 'BASE_EQ_MASS', 'ACID_EQ_MASS',
           'GYPSUM_EQ_MASS', 'TABLE', 'spec_for', 'parse_ops', 'join_ops', 'validate', 'apply_to_frame',
           'table_frame', 'water_removal_factor', 'broth_mass_fraction')

#: version of :data:`TABLE` (records with an older ``downstream_spec_version``
#: are re-patched)
SPEC_VERSION = 1

PHASES = ('liquid', 'solid', 'aqueous_solution')
ROUTES = ('evaporate_water', 'crystallize', 'distill_product', 'none')
SEPARATION_OPS = ('centrifugation', 'filtration', 'acidulation', 'gypsum_filtration', 'adsorption',
                  'ion_exchange', 'membrane', 'flash', 'evaporation', 'distillation', 'LLE', 'decantation',
                  'molecular_sieves', 'crystallization', 'drying', 'esterification', 'hydrolysis')
UPGRADING_OPS = ('dehydration', 'hydrogenation', 'ring_opening_hydrolysis', 'esterification', 'hydrolysis',
                 'decarboxylation', 'oxidation', 'amination', 'polymerization')
COPRODUCTS = ('electricity', 'DDGS', 'corn_oil', 'ethanol', 'lignin', 'gypsum', 'CO2', 'fertilizer')
ACIDS = ('H2SO4', 'none')

#: g base per mol of base equivalents (OH- or proton acceptor)
BASE_EQ_MASS = {'Ca(OH)2': 74.093 / 2, 'NaOH': 39.997, 'KOH': 56.106, 'NH4OH': 35.046, 'NH3': 17.031,
                'Mg(OH)2': 58.320 / 2, 'CaCO3': 100.087 / 2}
#: g acid per mol H+ (acid regeneration)
ACID_EQ_MASS = {'H2SO4': 98.079 / 2}
#: g gypsum (CaSO4.2H2O) per mol equivalent (calcium base + H2SO4)
GYPSUM_EQ_MASS = 172.17 / 2

#: field -> unit / vocabulary (merged into features.FIELD_UNITS)
FIELDS = {
    'final_product_mass_fraction': 'kg product / kg final-product stream at baseline, in (0, 1]; '
                                   'None for a broth product (see final_product_is_broth)',
    'final_product_is_broth': 'bool: the final product is the fermentation broth (w_f ~ t/1000, '
                              'no water removal)',
    'final_product_phase': 'one of ' + ', '.join(PHASES),
    'water_removal_route': 'one of ' + ', '.join(ROUTES),
    'separation_ops': "'+'-joined subset of: " + ', '.join(SEPARATION_OPS),
    'upgrading_ops': "'+'-joined subset of: " + ', '.join(UPGRADING_OPS) + " or 'none'",
    'coproducts': "'+'-joined subset of: " + ', '.join(COPRODUCTS) + " or 'none'",
    'neutralization_degree': 'mol base-eq. / mol product acid groups converted to salt in fermentation, '
                             'in [0, 1]',
    'acid_regeneration': "acid regenerating the free acid downstream: 'H2SO4' or 'none'",
}
STRING_FIELDS = ('final_product_phase', 'water_removal_route', 'separation_ops', 'upgrading_ops',
                 'coproducts', 'acid_regeneration')
MULTI_FIELDS = {'separation_ops': SEPARATION_OPS, 'upgrading_ops': UPGRADING_OPS, 'coproducts': COPRODUCTS}
#: fields copied into the processed samples (dataset.SPEC_COLUMNS)
SAMPLE_FIELDS = tuple(FIELDS) + ('downstream_spec_version',)


def parse_ops(value, vocab=None):
    """'a+b' / 'a,b' / ['a', 'b'] -> sorted unique list (empty for '', None,
    NaN or 'none'). With ``vocab``, unknown terms raise ValueError."""
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return []
    items = value if isinstance(value, (list, tuple, set)) else re.split(r'[+,;|]', str(value))
    out = sorted({str(v).strip() for v in items if str(v).strip() and str(v).strip().lower() != 'none'})
    if vocab is not None:
        bad = [v for v in out if v not in vocab]
        if bad:
            raise ValueError(f'unknown term(s) {bad}; allowed: {", ".join(vocab)}')
    return out


def join_ops(items):
    """Sorted '+'-joined terms; 'none' for an empty set (an explicit value
    that survives CSV round trips, unlike '')."""
    items = sorted(set(items) - {'none', ''})
    return '+'.join(items) if items else 'none'


def validate(spec):
    """Raise ValueError if a downstream field of ``spec`` (dict) is outside its
    vocabulary / range; returns the normalized dict of the fields present."""
    out = {}
    for k, vocab in MULTI_FIELDS.items():
        v = spec.get(k)
        # None / NaN = not given (filled later); '' or 'none' = explicitly empty
        if v is None or (isinstance(v, float) and not math.isfinite(v)):
            continue
        out[k] = join_ops(parse_ops(v, vocab))
    for k, vocab in (('final_product_phase', PHASES), ('water_removal_route', ROUTES),
                     ('acid_regeneration', ACIDS)):
        v = spec.get(k)
        if v is None or (isinstance(v, float) and not math.isfinite(v)):
            continue
        if k == 'acid_regeneration' and str(v) == '':
            v = 'none'
        if str(v) not in vocab:
            raise ValueError(f'{k}={v!r}; allowed: {vocab}')
        out[k] = str(v)
    for k, lo, hi in (('final_product_mass_fraction', 0., 1.), ('neutralization_degree', 0., 1.)):
        v = spec.get(k)
        if v is None or (isinstance(v, float) and not math.isfinite(v)):
            continue
        v = float(v)
        if not (lo < v <= hi if k == 'final_product_mass_fraction' else lo <= v <= hi):
            raise ValueError(f'{k}={v} outside its range ({FIELDS[k]})')
        out[k] = v
    v = spec.get('final_product_is_broth')
    if v is not None and not (isinstance(v, float) and not math.isfinite(v)):
        if isinstance(v, str):
            if v.strip().lower() not in ('true', 'false', '1', '0', '1.0', '0.0'):
                raise ValueError(f'final_product_is_broth={v!r}; expected a bool')
            v = v.strip().lower() in ('true', '1', '1.0')
        out['final_product_is_broth'] = bool(v)
    return out


def broth_mass_fraction(t):
    """Approximate w_f of a broth product at titer ``t`` [g/L] (t/1000)."""
    return t / 1000.


def water_removal_factor(route):
    """1 if broth water must be removed from the product (any route but
    'none'), 0 for 'none', NaN if unknown."""
    if route is None or (isinstance(route, float) and not math.isfinite(route)) or str(route) == '':
        return float('nan')
    return 0. if str(route) == 'none' else 1.


# %% Table

_FS4 = ('glucose', 'corn', 'sugarcane', 'cornstover')
_SRC_HP = ('adapter HP.py (separation_principle, separation_variant); BIP HP/process_areas.py '
           'create_HP_separation_improved_process (S401 centrifuge, R401 acidulation + S402 gypsum filter '
           'with lime only, adsorption, cation/anion exchange, F403 multi-effect evaporator to ~30 wt% HP), '
           '..._hexanol_extraction_process (F401 evaporator, S404 mixer-settlers, D401 distillation), '
           'create_HP_to_acrylic_acid_upgrading_process (R402 TiO2 dehydration, D408/D409 columns); '
           'w_f = x_AA_purity (0.9967 at baseline)')
_SRC_TAL = ('adapter TAL.py; BIP TAL/process_areas.py create_TAL_separation_solubility_exploit_process '
            '(S401 centrifuge, F401 evaporator, C401 crystallizer, S402 centrifuge, F402 drum dryer, F403 + '
            'D401 mother-liquor distillation), create_TAL_to_sorbic_acid_upgrading_process (R401 '
            'hydrogenation, R402 dehydration, R403 ring-opening hydrolysis with KOH, evaporators, drum '
            'dryers); w_f = baseline median of x_purity (TAL: varies with titer)')
_SRC_SA = ('adapter succinic.py (partial pH control 0.42 mol OH-eq/mol SA in the low-pH case); BIP '
           'succinic/system_*.py (S401 centrifuge, R401 acidulation + S405 gypsum filter, 3 x [F40x '
           'evaporator, C40x crystallizer, S40x pressure filter], F404 drum dryer); w_f = x_product_purity')
_SRC_LA = ('adapter lactic.py (product_purity_spec 0.88, separation_variant); BIP lactic/_processes.py '
           '(S401 cell-mass filter, R401 acidulation + S402 gypsum filter in the neutral case, F401 flash, '
           'D401, R402 esterification with ethanol, D402/D403 distillation, R403 hydrolysis, D404, F402)')
_SRC_ETH = ('adapter ethanol.py (separation_principle distillation+molecular_sieves); BIP corn / '
            'cellulosic / cane ethanol modules (beer column + rectifier, molecular sieves; corn: DDGS '
            'train and corn-oil recovery, no turbogenerator; cornstover: stillage pressure filter, lignin '
            'to the boiler-turbogenerator; cane: bagasse to the boiler-turbogenerator); w_f = '
            'x_product_purity (denatured fuel ethanol)')
_SRC_OX = ('adapter oxalic.py (product_form clarified_broth, neutralization_degree 0.91 mol Ca(OH)2/mol, '
           'centrifugation+acidulation+gypsum_filtration); w_f = x_purity follows the titer (broth)')
_SRC_IBO = ('adapter isobutanol.py (coproducts Ethanol,DDGS,crude_oil; heteroazeotropic distillation + '
            'decantation); BIP isobutanol/separations.py (ethanol molecular sieves, DDGS evaporator and '
            'dryer); x_purity_final_product = 1.0; CHP generates its own electricity demand')


def _hp(variant, feed):
    neutral = 'neutral' in variant
    hexanol = 'hexanol' in variant
    ops = ['centrifugation', 'evaporation', 'distillation'] + (['LLE'] if hexanol else ['adsorption', 'ion_exchange'])
    if neutral:
        ops += ['acidulation', 'gypsum_filtration']
    return dict(final_product_mass_fraction=0.9967, final_product_is_broth=False, final_product_phase='liquid',
                water_removal_route='evaporate_water', separation_ops=join_ops(ops), upgrading_ops='dehydration',
                coproducts='electricity', neutralization_degree=1. if neutral else 0.,
                acid_regeneration='H2SO4' if neutral else 'none', downstream_sources=_SRC_HP)


_TAL_WF = {'glucose': 0.896, 'corn': 0.873, 'sugarcane': 0.896, 'cornstover': 0.818}


def _tal(variant, feed):
    up = variant == 'TAL_SA'
    ops = ['centrifugation', 'evaporation', 'crystallization', 'drying', 'distillation']
    return dict(final_product_mass_fraction=0.949 if up else _TAL_WF[feed], final_product_is_broth=False,
                final_product_phase='solid', water_removal_route='crystallize', separation_ops=join_ops(ops),
                upgrading_ops='dehydration+hydrogenation+ring_opening_hydrolysis' if up else 'none',
                coproducts='electricity', neutralization_degree=0., acid_regeneration='none',
                downstream_sources=_SRC_TAL)


def _succinic(variant, feed):
    neutral = variant == 'succinic_neutral'
    ops = ['centrifugation', 'evaporation', 'crystallization', 'filtration', 'drying']
    if neutral:
        ops += ['acidulation', 'gypsum_filtration']
    return dict(final_product_mass_fraction=0.8487, final_product_is_broth=False, final_product_phase='solid',
                water_removal_route='crystallize', separation_ops=join_ops(ops), upgrading_ops='none',
                coproducts='electricity', neutralization_degree=1. if neutral else 0.21,
                acid_regeneration='H2SO4' if neutral else 'none', downstream_sources=_SRC_SA)


def _lactic(regime, kind):
    neutral = regime == 'neutral'
    ops = ['filtration', 'flash', 'distillation', 'esterification', 'hydrolysis']
    if neutral:
        ops += ['acidulation', 'gypsum_filtration']
    return dict(final_product_mass_fraction=0.88, final_product_is_broth=False,
                final_product_phase='aqueous_solution', water_removal_route='evaporate_water',
                separation_ops=join_ops(ops), upgrading_ops='none', coproducts='none',
                neutralization_degree=1. if neutral else 0., acid_regeneration='H2SO4' if neutral else 'none',
                downstream_sources=_SRC_LA)


_ETH = {
    'corn': dict(final_product_mass_fraction=0.9505,
                 separation_ops='centrifugation+distillation+drying+evaporation+molecular_sieves',
                 coproducts='DDGS+corn_oil'),
    'cornstover': dict(final_product_mass_fraction=0.9693, separation_ops='distillation+filtration+molecular_sieves',
                       coproducts='electricity'),
    'sugarcane': dict(final_product_mass_fraction=0.9704, separation_ops='distillation+molecular_sieves',
                      coproducts='electricity'),
}


def _ethanol(feed):
    return dict(final_product_is_broth=False, final_product_phase='liquid', water_removal_route='distill_product',
                upgrading_ops='none', neutralization_degree=0., acid_regeneration='none', downstream_sources=_SRC_ETH,
                **_ETH[feed])


def _build_table():
    t = {}
    for v in ('HP', 'HP_neutral', 'HP_hexanol', 'HP_neutral_hexanol'):
        for f in _FS4:
            t[f'{v}_{f}'] = _hp(v, f)
    for v in ('TAL', 'TAL_SA'):
        for f in _FS4:
            t[f'{v}_{f}'] = _tal(v, f)
    for v in ('succinic', 'succinic_neutral'):
        for f in _FS4:
            t[f'{v}_{f}'] = _succinic(v, f)
    for r in ('neutral', 'lowpH'):
        for k in ('SSCF', 'SHF'):
            t[f'lactic_{r}_{k}_cornstover'] = _lactic(r, k)
    for f in ('cornstover', 'sugarcane', 'corn'):
        t[f'ethanol_{f}'] = _ethanol(f)
    t['oxalic_neutral_sugarcane'] = dict(
        final_product_mass_fraction=None, final_product_is_broth=True, final_product_phase='aqueous_solution',
        water_removal_route='none', separation_ops='acidulation+centrifugation+gypsum_filtration',
        upgrading_ops='none', coproducts='electricity', neutralization_degree=0.91, acid_regeneration='H2SO4',
        downstream_sources=_SRC_OX)
    t['IBO_EtOH_corn'] = dict(
        final_product_mass_fraction=1.0, final_product_is_broth=False, final_product_phase='liquid',
        water_removal_route='distill_product',
        separation_ops='centrifugation+decantation+distillation+drying+evaporation+molecular_sieves',
        upgrading_ops='none', coproducts='DDGS+corn_oil+electricity+ethanol', neutralization_degree=0.,
        acid_regeneration='none', downstream_sources=_SRC_IBO)
    for cid, e in t.items():
        validate(e)
        e['downstream_spec_version'] = SPEC_VERSION
    return t


#: config_id -> downstream specification (all planned configurations)
TABLE = _build_table()


def spec_for(config_id):
    """Copy of the downstream specification of a known configuration (None if unknown)."""
    e = TABLE.get(str(config_id))
    return dict(e) if e is not None else None


def table_frame():
    """The table as a DataFrame (one row per configuration)."""
    import pandas as pd
    return pd.DataFrame.from_dict(TABLE, orient='index').rename_axis('config_id').reset_index()


def apply_to_frame(df, mark_col='data_patches', patch_id='S1'):
    """Copy of a processed table (``samples`` or ``configs``, with
    ``config_id``) whose downstream fields are set from :data:`TABLE` where the
    row's ``downstream_spec_version`` is missing or older (rows of unknown
    configurations are left unchanged). Idempotent; returns (df, n_rows_set)."""
    import numpy as np
    import pandas as pd
    df = df.copy()
    if 'config_id' not in df or not len(df):
        return df, 0
    ver = pd.to_numeric(df['downstream_spec_version'], errors='coerce') if 'downstream_spec_version' in df \
        else pd.Series(np.nan, index=df.index)
    cid = df['config_id'].astype(str)
    todo = cid.isin(list(TABLE)) & ~(ver >= SPEC_VERSION)
    if not todo.any():
        return df, 0
    fields = list(SAMPLE_FIELDS)
    for k in fields:
        if k not in df:
            df[k] = pd.Series([None] * len(df), index=df.index, dtype=object)
        elif df[k].dtype != object:
            df[k] = df[k].astype(object)
    for c in sorted(set(cid[todo])):
        m = todo & (cid == c)
        for k in fields:
            v = TABLE[c].get(k)
            df.loc[m, k] = pd.Series([v] * int(m.sum()), index=df.index[m], dtype=object)
    if mark_col:
        if mark_col not in df:
            df[mark_col] = ''
        df[mark_col] = [';'.join(sorted((set(str(v).split(';')) - {'', 'nan', 'None'}) | {patch_id}))
                        if t else v for v, t in zip(df[mark_col], todo)]
    for k in ('final_product_mass_fraction', 'neutralization_degree', 'downstream_spec_version'):
        df[k] = pd.to_numeric(df[k], errors='coerce')
    return df, int(todo.sum())
