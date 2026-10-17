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
>>> P = CoefficientPredictor.load('checkpoints/ensemble')            # doctest: +SKIP
>>> spec = dict(feedstock='cornstover', fermentation_product='HP',
...             final_product='AcrylicAcid', r=0.5, feedstock_price=0.08)
>>> P.predict(spec)[['a', 'b', 'c', 'd', 'g']]                      # doctest: +SKIP
>>> P.mpsp(spec, y=[0.3, 0.5], t=[50, 100])                         # doctest: +SKIP

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
unique for the product pair in training, baseline contextual ``p_*``
parameters, ...) are filled from the training data's registry
(``filled_from_registry``), and chemical descriptors come from
``registry/chemicals.json``. Unrecognized names are reported (with close
matches). New chemicals are described by passing
``descriptors={chemical_ID: {...}}`` (keys of the registry), by
``features.register_chemical`` or ``features.describe_new_chemical`` (needs
thermosteam), or by ``fp_*``/``pp_*`` fields in the spec. Omitted
categorical regime fields that cannot be filled are encoded as 'unknown'.

Output columns of :meth:`CoefficientPredictor.predict`
------------------------------------------------------
``a, b, c, d, g`` (ensemble mean; primed basis, a' = g*A), ``<k>_std``
(ensemble spread), ``<k>_lo``/``<k>_hi`` (split-conformal interval at level
1 - alpha), ``ch__<channel>__<k>`` (per-channel breakdown, sums to the
totals), ``int__<fs|el|ng>__<k>`` (intensities of price-explicit channels),
``ad_ratio`` (kNN distance / training threshold, in the encoder-feature +
productivity space), ``in_domain``, ``outside_range``,
``price_outside_range`` (informational: price channels are exactly linear),
``imputed`` (per row), ``contextual_defaults`` (contextual parameters not
given by the user: product baseline or training baseline) and
``filled_from_registry``. Without per-channel heads (``channels=False``) the
breakdown has a ``ch__rest__<k>`` channel. Intervals are omitted when the
calibration set was too small for the nominal level
(``ensemble.conformal['<k>_level']``).
"""
import difflib
import warnings

import numpy as np
import pandas as pd

from . import downstream
from .features import (MissingSpecError, ApplicabilityWarning, REQUIRED_FIELDS, FIELD_UNITS,
                       PRICE_ALIASES, W_SUGAR, W_SUGAR_MAX, G_MAX, apply_spec_registry,
                       warn_if_outside, _as_frame)
from .model import COEFF_KEYS, mpsp_equation_np
from .train import Ensemble, setup_torch

__all__ = ('CoefficientPredictor', 'MissingSpecError', 'ApplicabilityWarning')


class CoefficientPredictor:
    """Ensemble predictor of a', b', c', d', g with intervals."""

    def __init__(self, ensemble):
        self.ensemble = ensemble

    @classmethod
    def load(cls, path):
        """Load a checkpoint directory written by ``Ensemble.save``."""
        setup_torch()
        return cls(Ensemble.load(path))

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
        ``contextual``, ``imputed_extra``."""
        df = _as_frame(spec).reset_index(drop=True)
        reg = self.spec.spec_registry or {}
        params = [d for d in self.spec.numeric if d['kind'] == 'param' and d['name'].startswith('p_')]
        rows, det = [], []
        for rec in df.to_dict('records'):
            info = {}
            rec.update(downstream.validate(rec))     # vocabulary / range check (ValueError)
            rec, f = apply_spec_registry(rec, reg, info)
            ctx, extra = [c for c in f if c.startswith('p_') and c != 'p_feedstock_capacity'], []
            pair = f"{rec.get('fermentation_product')}|{rec.get('final_product')}"
            if pair not in reg.get('pair', {}):
                # unknown product pair: global training baseline, as an
                # applicable parameter (missing indicator 0)
                for d in params:
                    if _bad(rec.get(d['name'])):
                        rec[d['name']] = 10 ** d['fill'] if d['log'] else d['fill']
                        ctx.append(d['name'])
            if any(c in f for c in ('p_feedstock_capacity', 'fs_F_mass')):
                ctx.append('feedstock_capacity')
            elif _bad(rec.get('p_feedstock_capacity')) and _bad(rec.get('fs_F_mass')):
                ctx.append('feedstock_capacity')
                extra.append('feedstock_capacity')
            rows.append(rec)
            det.append(dict(filled=list(f), unknown=info.get('unknown', []), ambiguous=info.get('ambiguous', []),
                            contextual=sorted(set(ctx)), imputed_extra=extra))
        out = pd.DataFrame(rows)
        return (out, det) if details else (out, [d['filled'] for d in det])

    def missing_fields(self, spec):
        """Required fields still missing (or physically invalid) after the
        name-registry fill (per row)."""
        df, _ = self._frame(spec)
        return [self._missing_row(rec) for rec in df.to_dict('records')]

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
        keys = list((self.spec.spec_registry or {}).get(kind, {}))
        close = difflib.get_close_matches(str(name), keys, n=3, cutoff=0.5)
        return f"{kind} {name!r} is not in the registry" + (f" (did you mean {', '.join(map(repr, close))}?)"
                                                             if close else f" (known: {', '.join(keys[:12])})")

    def _batch(self, spec, descriptors, warn=True):
        df, det = self._frame(spec, details=True)
        miss = [self._missing_row(rec) for rec in df.to_dict('records')]
        bad = [(i, m) for i, m in enumerate(miss) if m]
        unknown = [(i, k, nm) for i, d in enumerate(det) for k, nm in d['unknown']]
        if bad:
            parts = '; '.join(f'row {i}: {", ".join(m)}' for i, m in bad[:10])
            unk = '; '.join(f'row {i}: ' + self._suggest(k, nm) for i, k, nm in unknown[:6])
            raise MissingSpecError(
                f'missing required specification field(s) ({parts}). Required: '
                + ', '.join(f'{c} [{FIELD_UNITS.get(c, "")}]' for c in REQUIRED_FIELDS)
                + (f'. Unrecognized names: {unk}' if unk else '')
                + ". Known feedstocks/fermentation products can be given "
                "by name ('feedstock', 'fermentation_product', 'final_product') to fill "
                "composition, yields and stoichiometric ratio from the registry; "
                "feedstock_price may also be given as 'p_feedstock_price'.")
        if warn and unknown:
            warnings.warn('; '.join(f'row {i}: ' + self._suggest(k, nm) for i, k, nm in unknown[:6])
                          + ' -- descriptors/composition not filled from the registry were imputed',
                          UserWarning, stacklevel=3)
        amb = [(i, d['ambiguous']) for i, d in enumerate(det) if d['ambiguous']]
        if warn and amb:
            warnings.warn('fermentation regime not unique for the product pair in training; please '
                          'specify: ' + '; '.join(f'row {i}: {", ".join(a)}' for i, a in amb[:6])
                          + " (treated as unknown)", UserWarning, stacklevel=3)
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
            warnings.warn('contextual parameters not specified, assumed baseline values: ' + ', '.join(ctx),
                          UserWarning, stacklevel=stacklevel)

    # ---- predictions ----
    def predict(self, spec, descriptors=None, warn=True):
        """Predict coefficients for one spec (dict) or many (DataFrame / list
        of dicts). Returns a DataFrame with one row per spec."""
        ens = self.ensemble
        df, batch, det = self._batch(spec, descriptors, warn)
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
        out['imputed'] = [', '.join(r) for r in batch.imputed_rows]
        out['contextual_defaults'] = [', '.join(d['contextual']) for d in det]
        out['filled_from_registry'] = [', '.join(d['filled']) for d in det]
        if warn:
            self._warn_imputed(batch, det)
        return out

    def _check_yt(self, df, y, t, warn):
        """Raise for y <= 0 / t <= 0; warn for y above the theoretical (MPSP
        set to NaN there: impossible) or inflection yield and t outside the
        training titer range. Returns the (n, P) mask of impossible yields."""
        if not (np.all(y > 0) and np.all(t > 0)):
            raise ValueError('y and t must be > 0 (y in g product / g glucose-eq. sugar, t in g/L)')
        num = lambda c: pd.to_numeric(df[c], errors='coerce').to_numpy(float) if c in df else np.full(len(df), np.nan)
        ty, infl = num('theoretical_yield'), num('inflection_yield')
        impossible = y[None, :] > ty[:, None] * (1 + 1e-9)
        msgs = []
        if impossible.any():
            msgs.append(f'y above the theoretical yield (MPSP set to NaN): rows {sorted(set(np.nonzero(impossible)[0]))[:10]}')
        above = (y[None, :] > infl[:, None] * (1 + 1e-9)) & ~impossible
        if above.any():
            msgs.append(f'y above the inflection yield (outside the fitted domain): rows '
                        f'{sorted(set(np.nonzero(above)[0]))[:10]}')
        dom = (self.spec.spec_registry or {}).get('domain', {})
        tout = []
        for i, rec in enumerate(df.to_dict('records')):
            d = dom.get(f"{rec.get('fermentation_product')}|{rec.get('final_product')}") or dom.get('__all__')
            if d and (np.any(t < d['t_min'] * (1 - 1e-9)) or np.any(t > d['t_max'] * (1 + 1e-9))):
                tout.append(i)
        if tout:
            msgs.append(f'titer outside the training range: rows {tout[:10]}')
        if warn and msgs:
            warnings.warn('; '.join(msgs), ApplicabilityWarning, stacklevel=3)
        return impossible

    def mpsp(self, spec, y, t, descriptors=None, warn=True):
        """MPSP [USD/kg] at yields ``y`` [g/g] and titers ``t`` [g/L]
        (broadcast against each other) for each spec row. Returns a long
        DataFrame with ``row, y, t, MPSP`` (ensemble mean of the members'
        MPSP), ``MPSP_std`` and the calibrated ``MPSP_lo``/``MPSP_hi``.
        Raises ValueError for y <= 0 or t <= 0; yields above the row's
        theoretical yield give NaN (with an ApplicabilityWarning)."""
        ens = self.ensemble
        y, t = np.broadcast_arrays(np.atleast_1d(np.asarray(y, float)), np.atleast_1d(np.asarray(t, float)))
        y, t = y.ravel(), t.ravel()
        df, batch, det = self._batch(spec, descriptors, warn)
        if warn and ens.ad is not None:
            ratio, oor = self._ad(batch)
            warn_if_outside(ratio, oor)
        if warn:
            self._warn_imputed(batch, det)
        impossible = self._check_yt(df, y, t, warn)
        outs = ens.member_outputs(batch)
        n = len(batch)
        Y, T = np.tile(y, (n, 1)), np.tile(t, (n, 1))
        mm = mpsp_equation_np(outs['coef'][:, :, None, :], outs['g'][:, :, None], Y[None], T[None])
        mm = np.where(impossible[None], np.nan, mm)
        mu, sd = mm.mean(0), mm.std(0)
        rows = dict(row=np.repeat(np.arange(n), len(y)), y=Y.ravel(), t=T.ravel(),
                    MPSP=mu.ravel(), MPSP_std=sd.ravel())
        q = (ens.conformal or {}).get('mpsp')
        if q is not None and np.isfinite(q):
            hw = q * (sd + ens.conformal.get('mpsp_eps_rel', 0.01) * np.abs(mu))
            rows['MPSP_lo'], rows['MPSP_hi'] = (mu - hw).ravel(), (mu + hw).ravel()
        return pd.DataFrame(rows)

    def unprimed(self, spec, **kw):
        """MPSP-level coefficients A, B, C, D = a'/g, ... (ensemble means)."""
        p = self.predict(spec, **kw)
        return pd.DataFrame({k.upper(): p[k] / p['g'] for k in COEFF_KEYS})


def _bad(v):
    if v is None:
        return True
    try:
        return not np.isfinite(float(v))
    except (TypeError, ValueError):
        return True
