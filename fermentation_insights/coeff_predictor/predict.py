#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
User-facing prediction of the unifying-equation coefficients.

Example
-------
>>> from fermentation_insights.coeff_predictor.predict import CoefficientPredictor
>>> P = CoefficientPredictor.load('models/release')                 # doctest: +SKIP
>>> spec = dict(feedstock='cornstover', fermentation_product='HP',
...             final_product='AcrylicAcid', r=0.5, feedstock_price=0.08)
>>> P.predict(spec)[['a', 'b', 'c', 'd', 'g']]                      # doctest: +SKIP
>>> P.mpsp(spec, y=[0.3, 0.5], t=[50, 100])                         # doctest: +SKIP
>>> P.mpsp(spec, y=0.5, t=100, recovery=0.7)   # known g            # doctest: +SKIP

See ``README.md`` (quickstart) and ``examples/quickstart.py``.

Input units and basis (``CoefficientPredictor.field_units``,
``features.FIELD_UNITS``)
------------------------------------------------------------------------
* ``r`` -- fermentation productivity [g/L/h], > 0 (required);
* ``feedstock_price`` -- USD per kg of WET (as-received) feedstock, not per
  dry tonne (required; alias ``p_feedstock_price`` / ``fs_price``);
* ``electricity_price`` [USD/kWh], ``natural_gas_price`` [USD/kg] (optional;
  training baseline if omitted -- reported in ``imputed``);
* ``fs_w_*`` -- wet-basis mass fractions [kg/kg wet feedstock];
  ``fs_w_fermentable_sugar_eq`` in (0, 1.12] (glucose-eq.);
* ``final_product_stoich_mass_ratio`` [kg final / kg fermentation product];
* ``theoretical_yield``, ``inflection_yield`` and ``y`` [g product / g
  glucose-eq. sugar]; ``t`` [g/L]; ``T_fermentation`` [K];
  ``p_feedstock_capacity`` [kg/h wet feedstock]; MPSP [USD/kg final product].

Final product and downstream process (``downstream.FIELDS``; used by
checkpoints trained with ``FeatureSpec(downstream=True)``, DESIGN.md section 5)
------------------------------------------------------------------------------
* ``final_product_mass_fraction`` w_f [kg product / kg final-product stream],
  in (0, 1]; for a product sold as the broth give
  ``final_product_is_broth=True`` (w_f ~ t/1000) instead;
* ``final_product_phase`` in {liquid, solid, aqueous_solution};
* ``water_removal_route`` in {evaporate_water, crystallize, distill_product,
  none} ('none' = the broth is sold; the model's water-removal gate then
  removes the titer-dependent downstream costs);
* ``separation_ops``, ``upgrading_ops``, ``coproducts``: '+'-joined terms of
  ``downstream.SEPARATION_OPS`` / ``UPGRADING_OPS`` / ``COPRODUCTS`` (or
  'none'), e.g. ``separation_ops='centrifugation+evaporation+crystallization+drying'``;
* ``neutralization_degree`` [mol base-eq. / mol product acid groups, 0-1] and
  ``acid_regeneration`` ('H2SO4' or 'none').

Values outside the vocabularies raise ValueError. Omitted fields are filled
from the registry when unique for the product pair in training, from
``downstream.TABLE`` when ``config_id`` names a known configuration (e.g.
``config_id='HP_hexanol_corn'``), and otherwise imputed (reported in
``imputed``; ``struct_water_removal`` then assumes water is removed).

Known feedstocks/products can be given by name: omitted fields (feedstock
composition and baseline capacity, theoretical/inflection yields,
stoichiometric mass ratio, chemical IDs, the fermentation regime when it is
unique, baseline contextual ``p_*`` parameters, the baseline electricity
and natural-gas prices of a known configuration, ...) are filled from the
training data's registry (``filled_from_registry``; precedence: a known
``config_id``, then the (product pair, feedstock) combination, then the
pair / feedstock / product entries -- ``features.apply_spec_registry``).
A ``config_id`` with a feedstock composition or type of its own and no
feedstock name describes a NEW feedstock in that configuration's process:
nothing of the configuration's feedstock is filled. A missing
``neutralizing_agent`` without neutralization (or acid group) is 'none', as
in training; boolean fields accept True/False, yes/no and 1/0.
Parameters not applicable to the matched training configuration stay NaN
as in training (``not_applicable``; ``features.NOT_APPLICABLE`` = 'n/a'
marks one explicitly), so a by-name query of a training configuration
encodes exactly as its training rows. Chemical descriptors come from
``registry/chemicals.json``. Unrecognized names are reported (with close
matches). New chemicals are described by passing
``descriptors={chemical_ID: {...}}`` (keys of the registry), by
``features.register_chemical`` or ``features.describe_new_chemical`` (needs
thermosteam), or by ``fp_*``/``pp_*`` fields in the spec. Omitted
categorical regime fields that cannot be filled are encoded as 'unknown'.
Regime fields that vary for the product pair in training should be given
(else they are encoded as unknown and a UserWarning names them); fields the
checkpoint does not use (the downstream
fields without ``FeatureSpec(downstream=True)``, else the free-text
``separation_principle``) are not asked for.

User-supplied recovery (``recovery=`` argument or ``recovery`` spec field)
--------------------------------------------------------------------------
If the recovery of the downstream process is known, pass it as
``recovery=`` (a scalar for every row, or one value per row; None/NaN = use
the prediction) or as a ``recovery`` field of the spec (the argument takes
precedence). It replaces the predicted g in ``g``, in the MPSP and in
:meth:`CoefficientPredictor.unprimed`, while the predicted primed
coefficients a', b', c', d' (and their channels and intervals) are kept.
It must satisfy 0 < g <= g_max = ``final_product_stoich_mass_ratio`` (the
stoichiometric ratio, after the registry fill), else ValueError.

Approximation: g * MPSP = (a' + b'/y + c'/t + d'/(y*t)) is the cost per kg of
*fermentation product*, which, at a fixed process design, does not depend on
how much of that product is recovered downstream; to first order, a
different recovery only changes the mass of final product over which the
same costs are spread, MPSP = N(y, t) / g. Not captured: costs that scale
with the recovered product mass (product purification / drying / upgrading
reagents, storage), the design changes (and their costs) that a different
recovery would need, and changes of coproduct credits from unrecovered
product. The MPSP intervals with a user recovery reflect the coefficient
uncertainty only (their conformal calibration included the predicted g's
error); they are not calibrated for this approximation.

Output columns of :meth:`CoefficientPredictor.predict`
------------------------------------------------------
``a, b, c, d, g`` (ensemble mean; primed basis, a' = g*A; ``g`` = the
user-supplied recovery where given), ``g_predicted`` (the model's g, ensemble
mean), ``g_source`` ('predicted' or 'user'), ``<k>_std``
(ensemble spread), ``<k>_lo``/``<k>_hi`` (split-conformal interval at level
1 - alpha), ``ch__<channel>__<k>`` (per-channel breakdown, sums to the
totals), ``int__<fs|el|ng>__<k>`` (intensities of price-explicit channels),
``ad_ratio`` (kNN distance / training threshold, in the encoder-feature +
productivity space), ``in_domain``, ``outside_range``,
``price_outside_range`` (informational: price channels are exactly linear),
``intervals_calibrated`` (False when the feedstock, fermentation product or
product pair is not in training: the intervals are then not calibrated),
``intervals_undercovered`` (the interval keys whose shipped quantiles cover
the row's fermentation-product family more than 0.1 below nominal in
development CV, with the measured coverage, e.g. ``'c 0.76, g 0.40'`` for a
known TAL configuration; '' otherwise; see "Interval calibration" below),
``imputed`` (per row), ``contextual_defaults`` (contextual parameters not
given by the user; their origin in ``contextual_source``), ``not_applicable``
(contextual parameters kept NaN), ``contextual_outside_range`` (contextual
parameters, incl. ``feedstock_capacity``, outside the training range of the
row's configuration or (pair, feedstock) combination, or never defined for
it; with an ApplicabilityWarning; not part of ``in_domain``, whose range
check is global) and ``filled_from_registry``. A given feedstock capacity
(``p_feedstock_capacity``, alias ``fs_F_mass``) must be a positive number
(ValueError otherwise). Without per-channel heads (``channels=False``) the
breakdown has a ``ch__rest__<k>`` channel. Intervals are omitted when the
calibration set was too small for the nominal level
(``ensemble.conformal['<k>_level']``). With a user recovery, ``g_std`` = 0 and
``g_lo`` = ``g_hi`` = ``g`` on that row. For a DataFrame spec the output keeps
the spec's index (rows in the same order); :meth:`CoefficientPredictor.mpsp`
numbers the spec rows positionally (column ``row`` = 0, 1, ...). Every given
numeric field (``r``, prices and their aliases, contextual ``p_*``
parameters, ``fs_w_*`` fractions, ``theoretical_yield`` etc.) must be a
number, None / NaN / '' (not given) or, for the contextual parameters, 'n/a'
(not applicable); anything else (e.g. ``p_IRR='10%'``) raises ValueError.
An ``fs_w_*`` key that is not a composition column of the checkpoint is
ignored (and reported), so it does not turn a by-name configuration into a
custom feedstock.

Interval calibration
--------------------
``<k>_lo``/``<k>_hi`` = mean -/+ q_k * (spread + eps_k) with the
split-conformal quantile q_k of ``ensemble.conformal`` (MPSP: relative
errors, simultaneous over a sample's titer-yield points). The checkpoints of
the experiment pipeline (``models/final_dev``, ``models/release``; protocol v2
and later) have no calibration hold-out: q_k comes from the out-of-fold scores of the
development leave-one-configuration-out CV of the same training procedure
and ensemble size (``conformal['source']``, ``['n']``,
``['n_calib_groups']``, ``['weighting']``), so the intervals target a new
configuration of known feedstocks / products. ``conformal['coverage_guaranteed']``
is False (the scores come from the CV fold ensembles), and
``conformal['crossfit']`` holds a cross-fitted coverage estimate per key
(mean, pooled, min, q10 over the held-out configurations). Coverage is
marginal over the samples of new configurations (configurations weighted by
their numbers of samples) and varies between configurations and families.
``conformal['<k>_level']`` is then nominal (``conformal['level_kind']``).
``ensemble.meta['conformal_random_cv']`` (when present) holds, for
information only, the quantiles calibrated on the random-CV scores.

The quantiles are pooled over all families, and the coverage differs
strongly between them (REPORT.md section 6). For a KNOWN configuration
(random-CV out-of-fold scores) the intervals are conservative pooled over
samples (MPSP 0.99, g 0.92 for nominal 0.90) but not per family: g covers
about 0.44 of the TAL samples (every TAL configuration; its recovery varies
with a module parameter that is not an input), c' about 0.76 of them and d'
about 0.56 of the lactic ones, while HP, ethanol, isobutanol, oxalic and
succinic acid are covered at 1.0 on every key, i.e. their intervals are too
wide (MPSP half-width typically 15-30 % of the MPSP for a known
configuration, for a random-CV error of about 3.5 %). For a NEW
configuration of a known product (LOCO scores) the MPSP intervals
under-cover ethanol (0.77) and lactic acid (0.69), g under-covers TAL
(0.43) and lactic acid (0.65), and b', c', d' under-cover ethanol and
lactic acid (d' also oxalic acid). ``ensemble.meta['coverage_by_family']``
(written by the experiment pipeline from the development out-of-fold
scores; report-only metadata, the quantiles stay pooled) holds the coverage
of the checkpoint's quantiles per fermentation product and key for a known
configuration (``'known_configuration'``: the lower of the random-CV and
productivity-tail coverage) and for a new one (``'new_configuration'``:
LOCO). ``predict`` / ``mpsp`` list the keys whose coverage for the row's
family is more than :data:`UNDERCOVERAGE_TOL` below nominal in
``intervals_undercovered`` (e.g. ``'c 0.76, g 0.40'`` for a known TAL
configuration; a row named by a training ``config_id``, without its own
feedstock or conflicting names, counts as a known configuration, any other
row of a known product as a new one) and issue an ApplicabilityWarning
quoting the numbers; their intervals are too narrow for that family.
"""
import difflib
import json
import math
import os
import warnings

import numpy as np
import pandas as pd

from . import downstream
from .features import (MissingSpecError, ApplicabilityWarning, REQUIRED_FIELDS, FIELD_UNITS,
                       PRICE_ALIASES, PRICE_COLUMNS, W_SUGAR, W_SUGAR_MAX, G_MAX, DOWNSTREAM_FIELDS,
                       FREE_TEXT_CATEGORICAL, LIGNOCELLULOSE_ONLY, SPEC_REGISTRY_VERSION, SPEC_BOOLEAN,
                       SPEC_CATEGORICAL, SPEC_NUMERIC, CANONICAL_PARAMETERS, EXCLUDED_ALWAYS, FEEDSTOCK_TYPES,
                       FeatureSpec, apply_spec_registry, build_spec_registry, is_not_applicable, parse_bool,
                       warn_if_outside, _as_frame, _feedstock_type, _jsonable, _omitted,
                       CAPACITY_KEY, FREE_SUGARS, STRUCTURAL_SUGARS)
from .model import COEFF_KEYS, mpsp_equation_np
from .train import Ensemble, setup_torch

__all__ = ('CoefficientPredictor', 'MissingSpecError', 'ApplicabilityWarning', 'refresh_spec_registry',
           'FERMENTATION_MODES')

#: documented values of ``fermentation_mode`` (README section 5); like
#: ``feedstock_type`` (``features.FEEDSTOCK_TYPES``) it is checked on user
#: input with case, blanks and '-' / '_' / ' ' ignored (ValueError with close
#: matches otherwise). A documented value not seen in training (e.g.
#: 'continuous') is encoded as the unknown level, with a warning.
FERMENTATION_MODES = ('batch', 'fed-batch', 'continuous')
#: price relative difference above which a price and its alias conflict
_ALIAS_RTOL = 1e-3
#: column prefixes of the training tables and of derived features that a
#: spec may carry without being user fields (never reported as unknown);
#: ``fs_w_*`` keys are checked against the checkpoint's composition columns
#: (:meth:`CoefficientPredictor._known_fields`)
_PASS_PREFIXES = ('m__', 'q__', 'int__', 'ch__', 'x_', 'bp_', 'c__', 'fp_', 'pp_', 'phys_', 'ds_', 'diag_')
#: the two names of the feedstock capacity [kg/h wet feedstock] (the
#: training tables carry both; ``phys_log10_capacity`` reads
#: ``p_feedstock_capacity``, else ``fs_F_mass``)
_CAPACITY_FIELDS = ('p_feedstock_capacity', 'fs_F_mass')
#: relative tolerance of the contextual training-range check
_CTX_RTOL = 1e-6
#: an interval key counts as under-covering a family when its measured
#: coverage (``ensemble.meta['coverage_by_family']``) is more than this below
#: the nominal level 1 - alpha
UNDERCOVERAGE_TOL = 0.1
#: coverage schemes of ``meta['coverage_by_family']`` and their labels
_COVERAGE_SCHEMES = {'known_configuration': 'known configuration: min of random CV and productivity tails',
                     'new_configuration': 'new configuration: leave-one-configuration-out CV'}
#: names, query and bookkeeping fields accepted without being model inputs:
#: registry names, the y / t query, outputs of ``predict``/``mpsp`` and the
#: non-input columns of the processed training tables (so training rows can
#: be passed back as specs)
_PASS_FIELDS = frozenset((
    'config_id', 'feedstock', 'fermentation_product', 'final_product', 'recovery', 'y', 't', 'MPSP',
    'ferm_product_chemical', 'final_product_chemical', 'product_water_solubility',
    'g_predicted', 'g_source', 'ad_ratio', 'in_domain', 'outside_range', 'price_outside_range',
    'intervals_calibrated', 'intervals_undercovered', 'imputed', 'contextual_defaults', 'contextual_source', 'not_applicable',
    'contextual_outside_range', 'filled_from_registry', 'MPSP_lo', 'MPSP_hi', 'MPSP_std',
    'sample_id', 'biorefinery', 'sample_kind', 'seed', 'sample_index', 'source_file', 'separation_variant',
    'y_min', 'y_max', 't_min', 't_max', 'baseline_productivity', 'baseline_yield', 'baseline_titer',
    'downstream_spec_version', 'inflection_yield_config', 'n_points_raw', 'n_ok_raw', 'n_infeasible', 'n_failed',
    'r_mult', 'harness_version', 'adapter_version', 'biosteam_fixes', 'data_patches', 'R2', 'n_points', 'fit_mape', 'fit_rmse',
    'fit_cond', 'phi', 'F_fixed', 'closure_rel', 'sum_check', 'A_mpsp', 'B_mpsp', 'C_mpsp', 'D_mpsp',
    'channels_ok', 'theoretical_yield', 'inflection_yield', 'fs_F_mass', 'fs_price'))


def _norm_token(v):
    return str(v).strip().lower().replace('_', '-').replace(' ', '-')


class CoefficientPredictor:
    """Ensemble predictor of a', b', c', d', g with intervals."""

    def __init__(self, ensemble):
        self.ensemble = ensemble

    @classmethod
    def load(cls, path, threads=None):
        """Load a checkpoint directory written by ``Ensemble.save``.
        ``threads``: if given, ``torch.set_num_threads(threads)`` (a
        process-wide setting; by default the session's setting is kept)."""
        if threads is not None:
            setup_torch(threads)
        self = cls(Ensemble.load(path))
        reg = self.spec.spec_registry or {}
        if reg and reg.get('version', 1) < SPEC_REGISTRY_VERSION:
            warnings.warn(f'{path}: the spec registry of this checkpoint predates the per-configuration '
                          'defaults (version < 2): by-name queries of training configurations use product-pair '
                          'contextual baselines; refresh it with predict.refresh_spec_registry(path, training_samples)',
                          UserWarning, stacklevel=2)
        return self

    def save(self, path):
        return self.ensemble.save(path)

    @property
    def required_fields(self):
        return REQUIRED_FIELDS

    @property
    def field_units(self):
        """Units / basis of the input fields (see the module docstring)."""
        return dict(FIELD_UNITS)

    @property
    def spec(self):
        return self.ensemble.spec

    @property
    def downstream_vocabulary(self):
        """Controlled vocabularies of the downstream specification fields."""
        return dict(final_product_phase=downstream.PHASES, water_removal_route=downstream.ROUTES,
                    separation_ops=downstream.SEPARATION_OPS, upgrading_ops=downstream.UPGRADING_OPS,
                    coproducts=downstream.COPRODUCTS, acid_regeneration=downstream.ACIDS)

    @staticmethod
    def known_downstream(config_id=None):
        """Downstream specification of a known configuration (dict), or the
        whole table (DataFrame) if ``config_id`` is None."""
        return downstream.table_frame() if config_id is None else downstream.spec_for(config_id)

    # ---- inputs ----
    def _frame(self, spec, details=False):
        """Rows after the registry fill (+ contextual defaults). With
        ``details``: also per-row dicts ``filled``, ``unknown``, ``ambiguous``,
        ``contextual``, ``context_source``, ``not_applicable``,
        ``cross_feedstock``, ``conflict``, ``imputed_extra``,
        ``context_outside`` ([(name, value, lo, hi)] outside the training
        range of ``context_ref``, the row's configuration or (pair,
        feedstock) combination)."""
        df = _as_frame(spec).reset_index(drop=True)
        if not len(df):
            raise ValueError('no specification rows: pass a dict, a list of dicts or a DataFrame with at least '
                             'one row (required fields: ' + ', '.join(REQUIRED_FIELDS) + ')')
        reg = self.spec.spec_registry or {}
        params = [d for d in self.spec.numeric if d['kind'] == 'param' and d['name'].startswith('p_')]
        # regime fields this checkpoint does not use are never reported as ambiguous
        unused = set(FREE_TEXT_CATEGORICAL) if self.spec.downstream else set(DOWNSTREAM_FIELDS)
        prices = set(PRICE_COLUMNS) | {a for v in PRICE_ALIASES.values() for a in v}
        rows, det = [], []
        for rec in _records(df):
            info = {}
            rec.update(downstream.validate(rec))     # vocabulary / range check (ValueError)
            unseen = self._check_regime(rec)          # feedstock_type / fermentation_mode / agent (ValueError)
            alias = self._alias_conflicts(rec)
            for c in SPEC_BOOLEAN:                   # True/False, yes/no, 1/0 (else ValueError)
                if c in rec:
                    rec[c] = parse_bool(rec[c], c)
            cap_user = self._check_capacity(rec)     # positive number where given (else ValueError)
            self._check_numeric(rec)                 # numbers where given (else ValueError)
            self._drop_unknown_composition(rec)      # an fs_w_* key the checkpoint does not know: ignored
            rec, f = apply_spec_registry(rec, reg, info)
            # no neutralization (or nothing to neutralize) and no agent given:
            # 'none', as the training rows encode a missing agent
            if _omitted(rec.get('neutralizing_agent')) and (
                    rec.get('neutralization') is False or rec.get('neutralization') == 0
                    or (not _bad(rec.get('n_acid_groups')) and float(rec['n_acid_groups']) == 0)):
                rec['neutralizing_agent'] = 'none'
                f = list(f) + ['neutralizing_agent']
                info.setdefault('source', {})['neutralizing_agent'] = 'none (no neutralization / no acid group)'
            src = dict(info.get('source', {}))
            na = set(info.get('not_applicable', []))
            ctx, extra = [c for c in f if c.startswith('p_') and c != 'p_feedstock_capacity'], []
            pair = f"{rec.get('fermentation_product')}|{rec.get('final_product')}"
            if pair not in reg.get('pair', {}):
                # unknown product pair: global training baseline, as an
                # applicable parameter (missing indicator 0), except where
                # kept NaN (not applicable) and the lignocellulose-only
                # parameters of other feedstocks
                ligno = _feedstock_type(rec.get('feedstock_type')) == 'lignocellulose'
                for d in params:
                    c = d['name']
                    if c in na or is_not_applicable(rec.get(c)) or (c in LIGNOCELLULOSE_ONLY and not ligno):
                        continue
                    if _bad(rec.get(c)):
                        rec[c] = 10 ** d['fill'] if d['log'] else d['fill']
                        ctx.append(c)
                        src[c] = 'global training baseline (unknown product pair)'
            # explicit 'not applicable' (features.NOT_APPLICABLE): NaN as in training
            for c in list(rec):
                if c.startswith('p_') and c not in prices and is_not_applicable(rec[c]):
                    rec[c] = float('nan')
                    na.add(c)
            if cap_user:
                pass                                  # given by the user: not a default
            elif any(c in f for c in _CAPACITY_FIELDS):
                ctx.append('feedstock_capacity')
                src['feedstock_capacity'] = src.get('p_feedstock_capacity', src.get('fs_F_mass'))
            elif _bad(rec.get('p_feedstock_capacity')) and _bad(rec.get('fs_F_mass')):
                ctx.append('feedstock_capacity')
                extra.append('feedstock_capacity')
                src['feedstock_capacity'] = 'training median (imputed)'
            ctx = sorted(set(ctx))
            fs_known = not _bad_name(rec.get('feedstock')) and str(rec['feedstock']) in reg.get('feedstock', {})
            fp_known = str(rec.get('fermentation_product')) in reg.get('fermentation_product', {})
            ref, outside = self._context_outside(rec, info, reg)
            rows.append(rec)
            det.append(dict(filled=list(f), unknown=info.get('unknown', []),
                            ambiguous=[c for c in info.get('ambiguous', []) if c not in unused],
                            contextual=ctx, context_source=sorted({src[c] for c in ctx if src.get(c)}),
                            not_applicable=sorted(c for c in na if c.startswith('p_')),
                            cross_feedstock=list(info.get('cross_feedstock', [])),
                            conflict=list(info.get('conflict', [])), imputed_extra=extra,
                            custom_feedstock=list(info.get('custom_feedstock', [])),
                            unseen_levels=unseen, alias_conflict=alias,
                            config_id=None if _bad_name(rec.get('config_id')) else str(rec['config_id']),
                            context_ref=ref, context_outside=outside,
                            calibrated=bool(fs_known and fp_known and pair in reg.get('pair', {})) if reg else True))
        out = pd.DataFrame(rows)
        return (out, det) if details else (out, [d['filled'] for d in det])

    def _check_regime(self, rec):
        """Validate and normalize the categorical regime fields of a user
        record in place: ``feedstock_type`` against ``features.FEEDSTOCK_TYPES``
        and ``fermentation_mode`` against :data:`FERMENTATION_MODES` (case,
        blanks and '-' / '_' / ' ' ignored; ValueError with close matches
        otherwise), and the case of ``neutralizing_agent`` against the
        training levels. Returns [(field, value)] of the (valid) values not
        seen in training, which the encoder maps to the unknown level."""
        levels = {k: [v for v in lv if v != '__other__'] for k, lv in (self.spec.categoricals or {}).items()}
        unseen = []
        for k, vocab in (('feedstock_type', FEEDSTOCK_TYPES), ('fermentation_mode', FERMENTATION_MODES)):
            v = rec.get(k)
            if k not in rec or _omitted(v) or (isinstance(v, str) and v == ''):
                continue
            canon = {_norm_token(x): x for x in vocab}
            t = canon.get(_norm_token(v))
            if t is None:
                close = difflib.get_close_matches(_norm_token(v), list(canon), n=3, cutoff=0.6)
                raise ValueError(f'{k}={v!r} is not one of {", ".join(vocab)}'
                                 + (f" (did you mean {', '.join(repr(canon[c]) for c in close)}?)" if close else ''))
            rec[k] = t
            if k in levels and t not in levels[k]:
                unseen.append((k, t))
        v = rec.get('neutralizing_agent')
        if 'neutralizing_agent' in rec and isinstance(v, str) and v != '' and 'neutralizing_agent' in levels:
            canon = {x.strip().lower(): x for x in levels['neutralizing_agent']}
            t = canon.get(v.strip().lower())
            if t is not None:
                rec['neutralizing_agent'] = t
            else:
                unseen.append(('neutralizing_agent', v))
        return unseen

    @staticmethod
    def _alias_conflicts(rec):
        """[(price, [(name, value), ...], used)] where a price and its
        aliases (``features.PRICE_ALIASES``) are given with different values;
        the canonical name wins, then the aliases in order (as in
        ``FeatureSpec``)."""
        out = []
        for c, aliases in PRICE_ALIASES.items():
            given = [(n, float(rec[n])) for n in (c,) + tuple(aliases) if n in rec and not _bad(rec[n])]
            if len(given) > 1:
                v0 = given[0][1]
                if any(abs(v - v0) > _ALIAS_RTOL * max(abs(v), abs(v0), 1e-300) for _, v in given[1:]):
                    out.append((c, given, given[0][0]))
        return out

    @staticmethod
    def _check_capacity(rec):
        """Validate a user-given feedstock capacity (``p_feedstock_capacity``
        or its training-table name ``fs_F_mass``, kg/h wet feedstock): a
        positive finite number (ValueError otherwise, as ``r <= 0`` is
        rejected; the capacity is never 'not applicable'). ``fs_F_mass``
        given alone is used as ``p_feedstock_capacity`` (the two are equal in
        the training tables wherever both are defined; otherwise a
        registry-filled ``p_feedstock_capacity`` would take precedence over
        it). Returns the names the user gave."""
        given = {}
        for c in _CAPACITY_FIELDS:
            v = rec.get(c)
            # not given: None, '' or NaN (an infinite value is invalid)
            if v is None or (isinstance(v, str) and not v.strip()) or (isinstance(v, float) and math.isnan(v)):
                continue
            msg = (f'{c}={v!r}: the feedstock capacity must be a positive number [kg/h wet feedstock] '
                   '(omit it for the configuration / feedstock baseline)')
            if is_not_applicable(v):
                raise ValueError(msg + "; it is never 'not applicable'")
            try:
                x = float(v)
            except (TypeError, ValueError):
                raise ValueError(msg) from None
            if not (math.isfinite(x) and x > 0):
                raise ValueError(msg)
            given[c] = x
        if 'fs_F_mass' in given and 'p_feedstock_capacity' not in given:
            rec['p_feedstock_capacity'] = given['fs_F_mass']
        elif len(given) == 2 and abs(given['fs_F_mass'] / given['p_feedstock_capacity'] - 1) > 1e-6:
            warnings.warn(f"p_feedstock_capacity={given['p_feedstock_capacity']:g} and fs_F_mass="
                          f"{given['fs_F_mass']:g} differ -- using p_feedstock_capacity", UserWarning, stacklevel=5)
        return sorted(given)

    def _numeric_fields(self):
        """User fields that must be numbers where given: ``r``, the prices
        and their aliases, the contextual ``p_*`` parameters, the ``fs_w_*``
        fractions and :data:`features.SPEC_NUMERIC` (the feedstock capacity
        is checked by :meth:`_check_capacity`)."""
        nf = getattr(self, '_numeric', None)
        if nf is None:
            user = self._known_fields()[0]
            nf = self._numeric = frozenset(
                ({c for c in user if c.startswith(('p_', 'fs_w_'))} | set(PRICE_COLUMNS) | set(SPEC_NUMERIC)
                 | {a for v in PRICE_ALIASES.values() for a in v} | {'r', W_SUGAR}) - set(_CAPACITY_FIELDS))
        return nf

    def _check_numeric(self, rec):
        """ValueError for a given numeric field (:meth:`_numeric_fields`)
        whose value is not a number: not given are None, NaN and '' (and
        'n/a', not applicable, which the contextual parameters accept and
        the required fields report as missing); a value such as '10%' or
        '300 days' would otherwise be read as missing and silently encoded
        as 'not applicable' / imputed."""
        bad = []
        for c in sorted(self._numeric_fields() & set(rec)):
            v = rec[c]
            if v is None or isinstance(v, (bool, np.bool_)) or (isinstance(v, str) and not v.strip()) \
                    or is_not_applicable(v) or (np.ndim(v) == 0 and pd.isna(v)):
                continue
            try:
                float(v)
            except (TypeError, ValueError):
                bad.append(c)
        if bad:
            unit = lambda c: FIELD_UNITS.get(c, FIELD_UNITS.get('fs_w_*') if c.startswith('fs_w_') else None)
            raise ValueError('not a number: ' + '; '.join(f'{c}={rec[c]!r}' + (f' [{unit(c)}]' if unit(c) else '')
                                                        for c in bad)
                             + ' -- give a number in the field\'s unit (e.g. p_IRR=0.10, p_operating_days=330), '
                             'or None / NaN to omit it' + (", or 'n/a' for a contextual parameter that does not "
                                                           "apply" if any(c.startswith('p_') for c in bad) else ''))

    def _drop_unknown_composition(self, rec):
        """Remove the ``fs_w_*`` keys that are not composition columns of
        the checkpoint (reported as unknown fields, e.g. 'fs_w_watre'): the
        model does not read them, and they must not make a by-name
        configuration a custom feedstock (``features._custom_feedstock_fields``)."""
        _, known, fs_w_checked = self._known_fields()
        if fs_w_checked:
            for c in [c for c in rec if isinstance(c, str) and c.startswith('fs_w_') and c not in known]:
                rec.pop(c)

    @staticmethod
    def _context_outside(rec, info, reg):
        """(reference, [(name, value, lo, hi)]): the contextual parameters
        of the filled row outside the training range of its configuration
        (``config_id``, unless overridden by conflicting names or a custom
        feedstock) or, else, of its (product pair, feedstock) combination
        (``reg['context_range']``; checkpoints without it: no check). The
        feedstock capacity is checked as :data:`features.CAPACITY_KEY`.
        Reported in ``contextual_outside_range``, separately from
        ``in_domain`` (whose range check uses the global training range)."""
        cr = reg.get('context_range') or {}
        ref, rng = None, None
        cid = None if _bad_name(rec.get('config_id')) else str(rec['config_id'])
        if cid is not None and not info.get('conflict') and not info.get('custom_feedstock') \
                and cid in cr.get('config', {}):
            ref, rng = f'configuration {cid}', cr['config'][cid]
        else:
            key = f"{rec.get('fermentation_product')}|{rec.get('final_product')}|{rec.get('feedstock')}"
            if key in cr.get('pair_feedstock', {}):
                fp, pp, fs = key.split('|')
                ref, rng = f'{fp} -> {pp} on {fs}', cr['pair_feedstock'][key]
        if not rng:
            return ref, []
        out = []
        for c, (lo, hi) in rng.items():
            if c == CAPACITY_KEY:
                v = rec.get('p_feedstock_capacity')
                v = rec.get('fs_F_mass') if _bad(v) else v
            else:
                v = rec.get(c)
            if _bad(v):
                continue
            v = float(v)
            tol = _CTX_RTOL * max(abs(lo), abs(hi), 1e-12)
            if v < lo - tol or v > hi + tol:
                out.append((c, v, lo, hi))
        # given, but never defined for the reference in training (NaN there:
        # encoded as missing in its training rows)
        known = set().union(*map(set, cr.get('config', {}).values())) if cr.get('config') else set()
        for c in sorted(known - set(rng) - {CAPACITY_KEY}):
            if not _bad(rec.get(c)):
                out.append((c, float(rec[c]), float('nan'), float('nan')))
        return ref, out

    def _known_fields(self):
        """Fields a spec may carry: the user fields (inputs, aliases,
        documented fields) and the pass-through names (:data:`_PASS_FIELDS`)."""
        kf = getattr(self, '_known', None)
        if kf is None:
            user = (set(REQUIRED_FIELDS) | set(FIELD_UNITS) | set(PRICE_COLUMNS) | set(EXCLUDED_ALWAYS)
                    | {a for v in PRICE_ALIASES.values() for a in v} | set(DOWNSTREAM_FIELDS) | set(SPEC_BOOLEAN)
                    | set(SPEC_CATEGORICAL) | set(SPEC_NUMERIC) | set(CANONICAL_PARAMETERS)
                    | {d['name'] for d in self.spec.numeric} | set(self.spec.categoricals or {}))
            user.discard('fs_w_*')
            # composition columns: the checkpoint's feedstock registry (every
            # fs_w_* of its training tables) and the sugar components read by
            # features.sugar_availability; an fs_w_* key outside them is
            # ignored by the model and reported (a typo, e.g. 'fs_w_watre')
            fs_w = {k for e in ((self.spec.spec_registry or {}).get('feedstock') or {}).values()
                    for k in e if k.startswith('fs_w_')}
            if fs_w:
                user |= fs_w | {f'fs_w_{k}' for k in FREE_SUGARS + STRUCTURAL_SUGARS}
            outs = {f'{k}{s}' for k in COEFF_KEYS + ('g',) for s in ('', '_std', '_lo', '_hi')}
            kf = self._known = (frozenset(user), frozenset(user | outs | _PASS_FIELDS), bool(fs_w))
        return kf

    def _unknown_fields(self, columns):
        """[(name, [close matches])] of the spec keys that are neither user
        fields nor pass-through names: they are ignored by the model (e.g. a
        misspelled price would otherwise silently fall back to the
        configuration's baseline)."""
        user, known, fs_w_checked = self._known_fields()
        out = []
        for c in columns:
            c = str(c)
            if c in known or c.startswith(_PASS_PREFIXES) or (c.startswith('fs_w_') and not fs_w_checked):
                continue
            out.append((c, difflib.get_close_matches(c, sorted(user), n=2, cutoff=0.75)))
        return out

    @staticmethod
    def _format_unknown(unk, n=8):
        return '; '.join(f'{c!r}' + (f" (did you mean {' or '.join(map(repr, m))}?)" if m else '')
                         for c, m in unk[:n]) + (' ...' if len(unk) > n else '')

    def missing_fields(self, spec):
        """Required fields still missing (or physically invalid) after the
        name-registry fill (per row)."""
        df, _ = self._frame(spec)
        return [self._missing_row(rec) for rec in _records(df)]

    @staticmethod
    def _missing_row(rec):
        miss = []
        for c in REQUIRED_FIELDS:
            v = rec.get(c)
            if c == 'feedstock_price' and _bad(v):
                v = next((rec.get(a) for a in PRICE_ALIASES['feedstock_price'] if not _bad(rec.get(a))), None)
            if _bad(v):
                miss.append(c)
            elif (c == 'r' and float(v) <= 0) or (c == G_MAX and float(v) <= 0) \
                    or (c == W_SUGAR and not 0 < float(v) <= W_SUGAR_MAX):
                miss.append(f'{c} (invalid: must be {FIELD_UNITS[c]})')
        return miss

    def _suggest(self, kind, name):
        reg = self.spec.spec_registry or {}
        if kind == 'config_id':
            keys = sorted(set(reg.get('config', {})) | set(downstream.TABLE))
        else:
            keys = list(reg.get(kind, {}))
        close = difflib.get_close_matches(str(name), keys, n=3, cutoff=0.5)
        return f"{kind} {name!r} is not in the registry" + (f" (did you mean {', '.join(map(repr, close))}?)"
                                                             if close else f" (known: {', '.join(keys[:12])})")

    def _batch(self, spec, descriptors, warn=True):
        df, det = self._frame(spec, details=True)
        miss = [self._missing_row(rec) for rec in _records(df)]
        bad = [(i, m) for i, m in enumerate(miss) if m]
        unknown = [(i, k, nm) for i, d in enumerate(det) for k, nm in d['unknown']]
        unk_fields = self._unknown_fields(_as_frame(spec).columns)
        if bad:
            parts = '; '.join(f'row {i}: {", ".join(m)}' for i, m in bad[:10])
            unk = '; '.join(f'row {i}: ' + self._suggest(k, nm) for i, k, nm in unknown[:6])
            raise MissingSpecError(
                f'missing required specification field(s) ({parts}). Required: '
                + ', '.join(f'{c} [{FIELD_UNITS.get(c, "")}]' for c in REQUIRED_FIELDS)
                + (f'. Unrecognized names: {unk}' if unk else '')
                + ('. Unknown field(s), ignored: ' + self._format_unknown(unk_fields) if unk_fields else '')
                + ". Known feedstocks/fermentation products can be given "
                "by name ('feedstock', 'fermentation_product', 'final_product') to fill "
                "composition, yields and stoichiometric ratio from the registry; "
                "feedstock_price may also be given as 'p_feedstock_price'.")
        if warn and unknown:
            warnings.warn('; '.join(f'row {i}: ' + self._suggest(k, nm) for i, k, nm in unknown[:6])
                          + ' -- descriptors/composition not filled from the registry were imputed',
                          UserWarning, stacklevel=3)
        if warn:
            if unk_fields:
                warnings.warn('unknown specification field(s) IGNORED by the model: '
                              + self._format_unknown(unk_fields) + ' -- field reference: README section 5, '
                              'P.field_units', UserWarning, stacklevel=3)
            al = [(i, a) for i, d in enumerate(det) for a in d['alias_conflict']]
            if al:
                warnings.warn('; '.join(f'row {i}: ' + ' vs. '.join(f'{n}={v:g}' for n, v in g) + f' -- using {u}'
                                        for i, (c, g, u) in al[:6])
                              + ' (a price and its alias were both given, with different values)', UserWarning,
                              stacklevel=3)
            un = [(i, k, v) for i, d in enumerate(det) for k, v in d['unseen_levels']]
            if un:
                lv = self.spec.categoricals or {}
                warnings.warn('; '.join(f'row {i}: {k}={v!r} was not seen in training (training levels: '
                                        + ', '.join(x for x in lv.get(k, []) if x != '__other__') + ')'
                                        for i, k, v in un[:6])
                              + ' -- encoded as the unknown level', UserWarning, stacklevel=3)
        conf = [(i, c) for i, d in enumerate(det) for c in d['conflict']]
        if warn and conf:
            warnings.warn('; '.join(f'row {i}: config_id {cid!r} is ' + ', '.join(f'{k}={v!r}' for k, _, v in ks)
                                    + ' in training, the spec says ' + ', '.join(f'{k}={u!r}' for k, u, _ in ks)
                                    for i, (cid, ks) in conf[:6])
                          + ' -- the configuration\'s baseline values were not used (only its downstream '
                          'specification); fix config_id or the names', UserWarning, stacklevel=3)
        amb = [(i, d['ambiguous']) for i, d in enumerate(det) if d['ambiguous']]
        if warn and amb:
            warnings.warn('fermentation regime not unique in training for the product pair (and feedstock); '
                          'please specify (or give config_id): '
                          + '; '.join(f'row {i}: {", ".join(a)}' for i, a in amb[:6])
                          + " (treated as unknown)", UserWarning, stacklevel=3)
        cust = [(i, cid, ks) for i, d in enumerate(det) for cid, ks in d['custom_feedstock']]
        if warn and cust:
            warnings.warn('; '.join(f'row {i}: config_id {cid!r} with its own feedstock composition / type '
                                    f'({", ".join(ks[:4])}{" ..." if len(ks) > 4 else ""}) and no feedstock name'
                                    for i, cid, ks in cust[:6])
                          + ' -- treated as a NEW feedstock: the configuration fills only its process, product, '
                          'regime and downstream fields; feedstock components left out are imputed (training '
                          'median), not taken from the configuration\'s feedstock', UserWarning, stacklevel=3)
        twins = [(i, d['config_id'], self._config_twins().get(d['config_id']))
                 for i, d in enumerate(det) if d.get('config_id') and not d['conflict']]
        twins = [t for t in twins if t[2]]
        if warn and twins:
            warnings.warn('; '.join(f'row {i}: config_id {cid!r} encodes identically to {", ".join(map(repr, tw))}'
                                    for i, cid, tw in twins[:6])
                          + ' -- these training configurations differ only in settings that are not model inputs '
                          '(e.g. the saccharification mode SHF vs. SSCF of the lactic configurations), so the '
                          'model predicts the same coefficients for them', UserWarning, stacklevel=3)
        cross = [(i, d['cross_feedstock']) for i, d in enumerate(det) if d['cross_feedstock']]
        if warn and cross:
            warnings.warn('feedstock-dependent contextual parameter(s) taken from the training baseline of OTHER '
                          'feedstocks of the product pair (e.g. sugarcane mills were simulated with 180 operating '
                          'days, the others with 330): ' + '; '.join(f'row {i}: {", ".join(c)}' for i, c in cross[:6])
                          + ' -- specify them for this feedstock', UserWarning, stacklevel=3)
        batch = self.spec.transform(df, descriptors=descriptors, strict=True, missing_categorical='other')
        rows = batch.imputed_rows or [[] for _ in range(len(batch))]
        batch.imputed_rows = [sorted(set(r) | set(d['imputed_extra'])) for r, d in zip(rows, det)]
        batch.imputed = sorted({c for r in batch.imputed_rows for c in r})
        return df, batch, det

    def _ad(self, batch):
        """(ratio, out-of-range names) of the applicability domain (encoder
        features + productivity; checkpoints whose AD includes prices are
        scored on the full matrix)."""
        ad = self.ensemble.ad
        Z = self.spec.ad_matrix(batch)
        if ad.Z is not None and ad.Z.shape[1] != Z.shape[1]:
            Z = self.spec.full_matrix(batch)
        return ad.score(Z)

    def _warn_imputed(self, batch, det, stacklevel=3):
        names = batch.imputed
        ctx = sorted({c for d in det for c in d['contextual']})
        if names:
            warnings.warn('not specified, assumed training baseline/median values: ' + ', '.join(names),
                          UserWarning, stacklevel=stacklevel)
        if ctx:
            srcs = sorted({s for d in det for s in d['context_source']})
            warnings.warn('contextual parameters not specified, assumed baseline values: ' + ', '.join(ctx)
                          + (' (from: ' + '; '.join(srcs[:6]) + (' ...' if len(srcs) > 6 else '') + ')'
                             if srcs else ''), UserWarning, stacklevel=stacklevel)

    @staticmethod
    def _warn_context_outside(det, stacklevel=4):
        bad = [(i, d) for i, d in enumerate(det) if d.get('context_outside')]
        if not bad:
            return
        parts = []
        for i, d in bad[:6]:
            parts.append(f"row {i}: " + ', '.join(f'{c}={v:.4g} (training {lo:.4g}-{hi:.4g})' if lo == lo
                                                  else f'{c}={v:.4g} (not defined in training)'
                                                  for c, v, lo, hi in d['context_outside'][:4])
                         + f" for {d['context_ref']}")
        cap = any(c == CAPACITY_KEY for _, d in bad for c, *_ in d['context_outside'])
        warnings.warn('contextual parameter(s) outside the training range of the configuration (or product pair '
                      'and feedstock): ' + ' | '.join(parts)
                      + (' -- each configuration\'s capacity was sampled only about +/-10-20 % around its baseline, '
                         'and the model reproduces only about two thirds of the simulated economy of scale (README '
                         'section 7): do not read scale-up effects from it' if cap else
                         ' -- the response to these parameters is extrapolated')
                      + '; in_domain does not include this check (column contextual_outside_range)',
                      ApplicabilityWarning, stacklevel=stacklevel)

    def _warn_uncalibrated(self, det, stacklevel=4):
        bad = [i for i, d in enumerate(det) if not d['calibrated']]
        if bad and (self.ensemble.conformal or {}):
            cov = (self.ensemble.meta or {}).get('coverage_unseen') or {}
            num = '; '.join(f"{lab}: MPSP {cov[k]['mpsp']:.2f}, g {cov[k]['g']:.2f}"
                            for k, lab in (('feedstock', 'new feedstock'), ('product', 'new product'))
                            if isinstance(cov.get(k), dict) and {'mpsp', 'g'} <= set(cov[k]))
            nominal = 1 - (self.ensemble.conformal or {}).get('alpha', 0.1)
            warnings.warn(f'rows {bad[:10]}: feedstock, fermentation product or product pair not in training -- the '
                          'intervals are calibrated for new configurations of KNOWN feedstocks and products; for an '
                          'unseen feedstock or product family their coverage is far below nominal'
                          + (f' (coverage of these quantiles in leave-one-feedstock-out / leave-one-product-out '
                             f'development CV, nominal {nominal:.2f}: {num})' if num else '')
                          + '; REPORT.md section 6', ApplicabilityWarning, stacklevel=stacklevel)

    def _undercovered(self, df, det, keys, user_g=None):
        """Per row: (scheme, fermentation product, {key: coverage}) of the
        interval ``keys`` whose shipped quantiles cover the row's family
        (fermentation product) more than :data:`UNDERCOVERAGE_TOL` below
        nominal in development CV (``ensemble.meta['coverage_by_family']``;
        a known ``config_id`` -> random-CV coverage, else LOCO), or None
        (no intervals, no such metadata, or a feedstock / product not in
        training, which :meth:`_warn_uncalibrated` covers). ``user_g``: rows
        with a user recovery (their g interval is the given value)."""
        conf = self.ensemble.conformal or {}
        cbf = (self.ensemble.meta or {}).get('coverage_by_family') or {}
        out = [None] * len(det)
        if not conf or not cbf:
            return out
        nominal = float(cbf.get('nominal', 1 - conf.get('alpha', 0.1)))
        cfg = (self.spec.spec_registry or {}).get('config', {})
        fps = df['fermentation_product'].tolist() if 'fermentation_product' in df else [None] * len(det)
        for i, (fp, d) in enumerate(zip(fps, det)):
            if not d['calibrated'] or _bad_name(fp):
                continue
            known = (d.get('config_id') in cfg and not d['conflict'] and not d['custom_feedstock'])
            sch = 'known_configuration' if known else 'new_configuration'
            cov = (cbf.get(sch) or {}).get(str(fp)) or {}
            low = {k: float(cov[k]) for k in keys if k in cov and cov[k] is not None
                   and float(cov[k]) < nominal - UNDERCOVERAGE_TOL and conf.get(k) is not None
                   and not (k == 'g' and user_g is not None and user_g[i])}
            out[i] = (sch, str(fp), low)
        return out

    @staticmethod
    def _format_undercovered(u):
        return '' if not u or not u[2] else ', '.join(f'{k} {v:.2f}' for k, v in u[2].items())

    def _warn_undercovered(self, und, stacklevel=4):
        bad = [(i, u) for i, u in enumerate(und) if u and u[2]]
        if not bad:
            return
        nominal = float(((self.ensemble.meta or {}).get('coverage_by_family') or {}).get(
            'nominal', 1 - (self.ensemble.conformal or {}).get('alpha', 0.1)))
        groups = {}
        for i, (sch, fp, low) in bad:
            groups.setdefault((sch, fp, tuple(low.items())), []).append(i)
        parts = [f'rows {rows[:10]}: {fp} ({_COVERAGE_SCHEMES[sch]}) ' + ', '.join(f'{k} {v:.2f}' for k, v in low)
                 for (sch, fp, low), rows in list(groups.items())[:6]]
        warnings.warn('intervals that UNDER-cover the fermentation product\'s family in development CV (coverage of '
                      f'the shipped quantiles, nominal {nominal:.2f}): ' + '; '.join(parts)
                      + ' -- these intervals are too narrow for this family (column intervals_undercovered; '
                      'README section 4.5, REPORT.md section 6)', ApplicabilityWarning, stacklevel=stacklevel)

    def _config_twins(self):
        """{config_id: [other training configurations whose by-name query
        encodes to the same model input]} (computed once per predictor)."""
        tw = getattr(self, '_twins', None)
        if tw is not None:
            return tw
        self._twins = {}
        cids = sorted((self.spec.spec_registry or {}).get('config', {}))
        if len(cids) < 2:
            return self._twins
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                df, _ = self._frame([dict(config_id=c, r=1., feedstock_price=1.) for c in cids])
                X = self.spec.transform(df, strict=False, missing_categorical='other').X
        except Exception:
            return self._twins
        keys = [np.round(x, 9).tobytes() for x in X]
        for i, c in enumerate(cids):
            others = [cids[j] for j in range(len(cids)) if j != i and keys[j] == keys[i]]
            if others:
                self._twins[c] = others
        return self._twins

    @staticmethod
    def _user_recovery(df, batch, recovery=None):
        """Per-row user-supplied recovery g (NaN = use the prediction) from the
        ``recovery=`` argument (scalar or one value per row) or else the
        ``recovery`` spec field. Raises ValueError unless 0 < g <= g_max
        (``final_product_stoich_mass_ratio`` of the row, after the registry fill)."""
        n = len(df)
        if recovery is None:
            if 'recovery' not in df:
                return np.full(n, np.nan)
            vals, src = list(df['recovery']), "spec field 'recovery'"
        else:
            arr = np.asarray(recovery, dtype=object)
            if arr.ndim == 0:
                vals = [arr.item()] * n
            elif arr.ndim == 1 and len(arr) == n:
                vals = list(arr)
            else:
                raise ValueError(f'recovery= must be a scalar or one value per spec row ({n} rows); '
                                 f'got shape {arr.shape}')
            src = 'recovery argument'
        g = np.full(n, np.nan)
        for i, v in enumerate(vals):
            if v is None or (isinstance(v, str) and v.strip() == ''):
                continue
            try:
                x = float(v)
            except (TypeError, ValueError):
                raise ValueError(f'{src}: recovery of row {i} must be a number in (0, g_max] '
                                 f'[{FIELD_UNITS["recovery"]}]; got {v!r}') from None
            if math.isnan(x):
                continue                  # not given: predicted g
            g[i] = x
        given = ~np.isnan(g)
        gmax = np.asarray(batch.g_max, float)
        with np.errstate(invalid='ignore'):
            bad = given & ~(np.isfinite(g) & (g > 0) & (g <= gmax * (1. + 1e-9)))
        if bad.any():
            names = lambda i: f"{df['fermentation_product'].iloc[i]} -> {df['final_product'].iloc[i]}" \
                if {'fermentation_product', 'final_product'} <= set(df) else ''
            rows = '; '.join(f'row {i}: g = {g[i]:.6g}, g_max = {gmax[i]:.6g}'
                             + (f' ({names(i)})' if names(i) else '') for i in np.flatnonzero(bad)[:10])
            raise ValueError(
                f'{src}: the recovery g must satisfy 0 < g <= g_max, where g_max = '
                f'final_product_stoich_mass_ratio is the stoichiometric mass ratio of final to fermentation '
                f'product [kg final product / kg fermentation product]; {rows}')
        return g

    @staticmethod
    def _g_source(g_user):
        return np.where(np.isnan(g_user), 'predicted', 'user')

    # ---- predictions ----
    def predict(self, spec, descriptors=None, warn=True, recovery=None):
        """Predict coefficients for one spec (dict) or many (DataFrame / list
        of dicts). Returns a DataFrame with one row per spec.

        ``recovery``: optional known recovery g [kg final / kg fermentation
        product] (scalar, or one value per row; None/NaN = predicted); else a
        ``recovery`` spec field is used where given. It replaces the predicted
        g (kept as ``g_predicted``) and leaves a', b', c', d' unchanged (a
        first-order approximation, see the module docstring)."""
        ens = self.ensemble
        df, batch, det = self._batch(spec, descriptors, warn)
        g_user = self._user_recovery(df, batch, recovery)
        user = ~np.isnan(g_user)
        outs = ens.member_outputs(batch)
        conf = ens.conformal or {}
        res = {}
        coef, g = outs['coef'], outs['g']
        for k, key in enumerate(COEFF_KEYS):
            mu, sd = coef[:, :, k].mean(0), coef[:, :, k].std(0)
            res[key], res[f'{key}_std'] = mu, sd
            q = conf.get(key)
            if q is not None and np.isfinite(q):
                hw = q * (sd + conf[f'{key}_eps'])
                res[f'{key}_lo'], res[f'{key}_hi'] = mu - hw, mu + hw
        mu, sd = g.mean(0), g.std(0)
        res['g'], res['g_std'] = mu, sd
        q = conf.get('g')
        if q is not None and np.isfinite(q):
            hw = q * (sd + conf['g_eps'])
            res['g_lo'] = np.maximum(mu - hw, 0.)
            hi = mu + hw
            if ens.mcfg.g_bound:
                hi = np.minimum(hi, batch.g_max * (1. + ens.mcfg.g_margin))
            res['g_hi'] = hi
        res['g_predicted'] = mu.copy()
        if user.any():
            res['g'] = np.where(user, g_user, mu)
            res['g_std'] = np.where(user, 0., sd)
            for k in ('g_lo', 'g_hi'):
                if k in res:
                    res[k] = np.where(user, g_user, res[k])
        res['g_source'] = self._g_source(g_user)
        if 'channels' in outs:
            for c, name in enumerate(ens.channel_names()):
                for k, key in enumerate(COEFF_KEYS):
                    res[f'ch__{name}__{key}'] = outs['channels'][:, :, c, k].mean(0)
        if 'intensities' in outs:
            for j, name in enumerate(ens.layout['price_keys']):
                for k, key in enumerate(COEFF_KEYS):
                    res[f'int__{name}__{key}'] = outs['intensities'][:, :, j, k].mean(0)
        out = pd.DataFrame(res)
        if ens.ad is not None:
            ratio, oor = self._ad(batch)
            out['ad_ratio'] = ratio
            out['in_domain'] = np.isfinite(ratio) & (ratio <= 1.) & np.array([not o for o in oor])
            out['outside_range'] = [', '.join(o) for o in oor]
            if warn:
                warn_if_outside(ratio, oor)
        if self.spec.price_explicit:
            out['price_outside_range'] = [', '.join(o) for o in self.spec.price_outside_range(batch)]
        out['intervals_calibrated'] = [d['calibrated'] for d in det]
        und = self._undercovered(df, det, COEFF_KEYS + ('g',), user_g=user)
        out['intervals_undercovered'] = [self._format_undercovered(u) for u in und]
        if warn:
            self._warn_uncalibrated(det)
            self._warn_undercovered(und)
        out['imputed'] = [', '.join(r) for r in batch.imputed_rows]
        out['contextual_defaults'] = [', '.join(d['contextual']) for d in det]
        out['contextual_source'] = ['; '.join(d['context_source']) for d in det]
        out['not_applicable'] = [', '.join(d['not_applicable']) for d in det]
        out['contextual_outside_range'] = [', '.join(c for c, *_ in d['context_outside']) for d in det]
        out['filled_from_registry'] = [', '.join(d['filled']) for d in det]
        if warn:
            self._warn_imputed(batch, det)
            self._warn_context_outside(det)
        if isinstance(spec, pd.DataFrame):
            out.index = spec.index                    # rows in the spec's order
        return out

    def _check_yt(self, df, y, t, warn):
        """Raise for y <= 0 / t <= 0; warn for y above the theoretical (MPSP
        set to NaN there: impossible) or inflection yield, y below the
        training yield range of the product pair (registry ``y_domain``) and t
        outside its training titer range. Returns the (n, P) mask of
        impossible yields."""
        if not (np.all(y > 0) and np.all(t > 0)):
            raise ValueError('y and t must be > 0 (y in g product / g glucose-eq. sugar, t in g/L)')
        num = lambda c: pd.to_numeric(df[c], errors='coerce').to_numpy(float) if c in df else np.full(len(df), np.nan)
        ty, infl = num('theoretical_yield'), num('inflection_yield')
        impossible = y[None, :] > ty[:, None] * (1 + 1e-9)
        msgs = []
        if impossible.any():
            msgs.append(f'y above the theoretical yield (MPSP set to NaN): rows {sorted({int(i) for i in np.nonzero(impossible)[0]})[:10]}')
        above = (y[None, :] > infl[:, None] * (1 + 1e-9)) & ~impossible
        if above.any():
            msgs.append(f'y above the inflection yield (outside the fitted domain): rows '
                        f'{sorted({int(i) for i in np.nonzero(above)[0]})[:10]}')
        reg = self.spec.spec_registry or {}
        dom, ydom = reg.get('domain', {}), reg.get('y_domain', {})
        tout, ylow = [], []
        for i, rec in enumerate(df.to_dict('records')):
            key = f"{rec.get('fermentation_product')}|{rec.get('final_product')}"
            d = dom.get(key) or dom.get('__all__')
            if d and (np.any(t < d['t_min'] * (1 - 1e-9)) or np.any(t > d['t_max'] * (1 + 1e-9))):
                tout.append(i)
            dy = ydom.get(key) or ydom.get('__all__')
            if dy and np.any(y < dy['y_min'] * (1 - 1e-9)):
                ylow.append((i, dy['y_min']))
        if tout:
            msgs.append(f'titer outside the training range: rows {tout[:10]}')
        if ylow:
            msgs.append('y below the training yield range (extrapolated along the 1/y terms): '
                        + ', '.join(f'row {i} (training y >= {v:.3g})' for i, v in ylow[:10]))
        if warn and msgs:
            warnings.warn('; '.join(msgs), ApplicabilityWarning, stacklevel=3)
        return impossible

    def mpsp(self, spec, y, t, descriptors=None, warn=True, recovery=None):
        """MPSP [USD/kg] at yields ``y`` [g/g] and titers ``t`` [g/L]
        (broadcast against each other) for each spec row. Returns a long
        DataFrame with ``row, y, t, MPSP`` (ensemble mean of the members'
        MPSP), ``MPSP_std``, the calibrated ``MPSP_lo``/``MPSP_hi``,
        ``g_source``, ``intervals_calibrated``, ``intervals_undercovered``
        ('mpsp <coverage>' where the MPSP intervals under-cover the row's
        family, see :meth:`predict`) and ``contextual_outside_range``; ``row``
        numbers the spec rows positionally (0, 1, ..., also for a DataFrame
        spec). Raises ValueError for y <= 0 or t <= 0; yields above the row's
        theoretical yield give NaN (with an ApplicabilityWarning); yields
        below the training range of the product pair and titers outside it
        are flagged with an ApplicabilityWarning.
        ``recovery`` as in :meth:`predict`: every member's MPSP is then
        (a' + b'/y + c'/t + d'/(y t)) / recovery (intervals: coefficient
        uncertainty only)."""
        ens = self.ensemble
        y, t = np.broadcast_arrays(np.atleast_1d(np.asarray(y, float)), np.atleast_1d(np.asarray(t, float)))
        y, t = y.ravel(), t.ravel()
        df, batch, det = self._batch(spec, descriptors, warn)
        g_user = self._user_recovery(df, batch, recovery)
        if warn and ens.ad is not None:
            ratio, oor = self._ad(batch)
            warn_if_outside(ratio, oor)
        if warn:
            self._warn_imputed(batch, det)
            self._warn_context_outside(det, stacklevel=3)
        impossible = self._check_yt(df, y, t, warn)
        outs = ens.member_outputs(batch)
        n = len(batch)
        Y, T = np.tile(y, (n, 1)), np.tile(t, (n, 1))
        g = np.where(np.isnan(g_user)[None, :], outs['g'], g_user[None, :])
        mm = mpsp_equation_np(outs['coef'][:, :, None, :], g[:, :, None], Y[None], T[None])
        mm = np.where(impossible[None], np.nan, mm)
        mu, sd = mm.mean(0), mm.std(0)
        rows = dict(row=np.repeat(np.arange(n), len(y)), y=Y.ravel(), t=T.ravel(),
                    MPSP=mu.ravel(), MPSP_std=sd.ravel())
        q = (ens.conformal or {}).get('mpsp')
        if q is not None and np.isfinite(q):
            hw = q * (sd + ens.conformal.get('mpsp_eps_rel', 0.01) * np.abs(mu))
            rows['MPSP_lo'], rows['MPSP_hi'] = (mu - hw).ravel(), (mu + hw).ravel()
        rows['g_source'] = np.repeat(self._g_source(g_user), len(y))
        rows['intervals_calibrated'] = np.repeat([d['calibrated'] for d in det], len(y))
        und = self._undercovered(df, det, ('mpsp',))
        rows['intervals_undercovered'] = np.repeat([self._format_undercovered(u) for u in und], len(y))
        rows['contextual_outside_range'] = np.repeat([', '.join(c for c, *_ in d['context_outside'])
                                                      for d in det], len(y))
        if warn:
            self._warn_uncalibrated(det, stacklevel=3)
            self._warn_undercovered(und, stacklevel=3)
        return pd.DataFrame(rows)

    def unprimed(self, spec, **kw):
        """MPSP-level coefficients A, B, C, D = a'/g, ... (ensemble means; g =
        the user-supplied ``recovery`` where given, see :meth:`predict`)."""
        p = self.predict(spec, **kw)
        return pd.DataFrame({k.upper(): p[k] / p['g'] for k in COEFF_KEYS})


def refresh_spec_registry(path, samples, strict=True):
    """Rebuild the spec registry (the by-name defaults used by
    :class:`CoefficientPredictor`, ``features.build_spec_registry``) of the
    checkpoint ``path`` in the current layout from its training samples
    ``samples`` and rewrite ``feature_spec.json``; the networks, feature
    encoding, scaling, intervals and applicability domain are unchanged.
    The registry parts the checkpoint already has must be reproduced exactly
    by ``samples`` (i.e. they are its training samples): else ValueError
    (``strict``) or nothing is written. Returns True if the file changed."""
    fpath = os.path.join(path, 'feature_spec.json')
    spec = FeatureSpec.from_json(fpath)
    old = spec.spec_registry or {}
    new = json.loads(json.dumps(_jsonable(build_spec_registry(samples))))
    bad = [k for k in old if k != 'version' and old.get(k) != new.get(k)]
    if bad:
        if strict:
            raise ValueError(f'{path}: the samples do not reproduce the stored spec registry ({", ".join(bad)}); '
                             'pass the training samples of this checkpoint')
        return False
    if old == new:
        return False
    spec.spec_registry = new
    spec.to_json(fpath)
    return True


def _bad_name(v):
    return v is None or (isinstance(v, float) and not math.isfinite(v)) or str(v).strip() == ''


def _records(df):
    """Row dicts (also for a frame without columns, e.g. from ``{}``)."""
    return df.to_dict('records') if df.shape[1] else [{} for _ in range(len(df))]


def _bad(v):
    if v is None:
        return True
    try:
        return not np.isfinite(float(v))
    except (TypeError, ValueError):
        return True
