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

Known feedstocks/products can be given by name: omitted fields (feedstock
composition, theoretical/inflection yields, stoichiometric mass ratio,
chemical IDs, ...) are filled from the training data's registry, and chemical
descriptors come from ``registry/chemicals.json``. New chemicals are described
by passing ``descriptors={chemical_ID: {...}}`` (keys of the registry), by
``features.register_chemical`` or ``features.describe_new_chemical`` (needs
thermosteam), or by ``fp_*``/``pp_*`` fields in the spec.

Output columns of :meth:`CoefficientPredictor.predict`
------------------------------------------------------
``a, b, c, d, g`` (ensemble mean; primed basis, a' = g*A), ``<k>_std``
(ensemble spread), ``<k>_lo``/``<k>_hi`` (split-conformal interval at level
1 - alpha), ``ch__<channel>__<k>`` (per-channel breakdown, sums to the
totals), ``int__<fs|el|ng>__<k>`` (intensities of price-explicit channels),
``ad_ratio`` (kNN distance / training threshold), ``in_domain``,
``imputed`` and ``filled_from_registry``.
"""
import warnings

import numpy as np
import pandas as pd

from .features import (MissingSpecError, ApplicabilityWarning, REQUIRED_FIELDS,
                       apply_spec_registry, warn_if_outside, _as_frame)
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
    def spec(self):
        return self.ensemble.spec

    # ---- inputs ----
    def _frame(self, spec):
        df = _as_frame(spec).reset_index(drop=True)
        reg = self.spec.spec_registry or {}
        rows, filled = [], []
        for rec in df.to_dict('records'):
            rec, f = apply_spec_registry(rec, reg)
            rows.append(rec)
            filled.append(f)
        return pd.DataFrame(rows), filled

    def missing_fields(self, spec):
        """Required fields still missing after the name-registry fill (per row)."""
        df, _ = self._frame(spec)
        out = []
        for i, rec in enumerate(df.to_dict('records')):
            miss = []
            for c in REQUIRED_FIELDS:
                v = rec.get(c)
                if c == 'feedstock_price' and _bad(v):
                    v = rec.get('p_feedstock_price', rec.get('fs_price'))
                if _bad(v):
                    miss.append(c)
            out.append(miss)
        return out

    def _batch(self, spec, descriptors):
        df, filled = self._frame(spec)
        miss = self.missing_fields(df)
        bad = [(i, m) for i, m in enumerate(miss) if m]
        if bad:
            parts = '; '.join(f'row {i}: {", ".join(m)}' for i, m in bad[:10])
            raise MissingSpecError(
                f'missing required specification field(s) ({parts}). Required: '
                + ', '.join(REQUIRED_FIELDS) + ". Known feedstocks/fermentation products can be given "
                "by name ('feedstock', 'fermentation_product', 'final_product') to fill "
                "composition, yields and stoichiometric ratio from the registry; "
                "feedstock_price may also be given as 'p_feedstock_price'.")
        batch = self.spec.transform(df, descriptors=descriptors, strict=True)
        return df, batch, filled

    # ---- predictions ----
    def predict(self, spec, descriptors=None, warn=True):
        """Predict coefficients for one spec (dict) or many (DataFrame / list
        of dicts). Returns a DataFrame with one row per spec."""
        ens = self.ensemble
        df, batch, filled = self._batch(spec, descriptors)
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
            ratio, oor = ens.ad.score(self.spec.full_matrix(batch))
            out['ad_ratio'] = ratio
            out['in_domain'] = (ratio <= 1.) & np.array([not o for o in oor])
            out['outside_range'] = [', '.join(o) for o in oor]
            if warn:
                warn_if_outside(ratio, oor)
        out['imputed'] = ', '.join(batch.imputed)
        out['filled_from_registry'] = [', '.join(f) for f in filled]
        if warn and batch.imputed:
            warnings.warn('not specified, assumed training baseline/median values: ' + ', '.join(batch.imputed),
                          UserWarning, stacklevel=2)
        return out

    def mpsp(self, spec, y, t, descriptors=None, warn=True):
        """MPSP [USD/kg] at yields ``y`` [g/g] and titers ``t`` [g/L]
        (broadcast against each other) for each spec row. Returns a long
        DataFrame with ``row, y, t, MPSP`` (ensemble mean of the members'
        MPSP), ``MPSP_std`` and the calibrated ``MPSP_lo``/``MPSP_hi``."""
        ens = self.ensemble
        df, batch, _ = self._batch(spec, descriptors)
        if warn and ens.ad is not None:
            ratio, oor = ens.ad.score(self.spec.full_matrix(batch))
            warn_if_outside(ratio, oor)
        y, t = np.broadcast_arrays(np.atleast_1d(np.asarray(y, float)), np.atleast_1d(np.asarray(t, float)))
        y, t = y.ravel(), t.ravel()
        outs = ens.member_outputs(batch)
        n = len(batch)
        Y, T = np.tile(y, (n, 1)), np.tile(t, (n, 1))
        mm = mpsp_equation_np(outs['coef'][:, :, None, :], outs['g'][:, :, None], Y[None], T[None])
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
