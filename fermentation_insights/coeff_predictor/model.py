#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Physics-informed network for the coefficients of the unifying equation
(DESIGN.md section 3)::

    MPSP(y, t) = (a' + b'/y + c'/t + d'/(y*t)) / g

Architecture
------------
* encoder: MLP on [standardized features, standardized log r, standardized
  1/r] with ``depth`` blocks of Linear -> LayerNorm -> SiLU -> Dropout;
* price-explicit heads: intensity coefficients ``q_j`` (4 each) for the
  feedstock (``fs``), electricity (``el``) and natural gas (``ng``) channels;
  channel coefficients are ``p_j * q_j`` with the price ``p_j`` *not* seen by
  the encoder, so predictions are exactly linear in prices;
  the feedstock-intensity ``b'`` is structured as ``exp(z)/w_sugar`` with
  ``z = softplus(.) >= 0`` (``exp(z) = 1/eta``, eta = fraction of potential
  sugar reaching fermentation);
* channel heads: coefficients (4 each) of every other cost channel;
* totals ``a', b', c', d'`` = sum over channels (hard identity);
* recovery ``g = g_max * sigmoid(z_g)`` with ``g_max`` the stoichiometric mass
  ratio of final to fermentation product (hard upper bound; an optional small
  relative ``g_margin`` absorbs rounding-level excesses found in the data);
* equation layer: :func:`mpsp_equation`;
* optional input clipping (``clip_inputs >= 0``; off by default): standardized
  continuous encoder features are clipped to their training range, which bounds the
  network's extrapolation to descriptors outside the training products
  (productivity is never clipped, nor are the prices of the price-explicit
  channels, which do not enter the encoder; the plain MLP's price features are);
* optional water-removal gate (``water_gate``; needs a FeatureSpec with
  ``downstream=True``): the c' and d' coefficients of the channels that scale
  with the broth water processed downstream (``GATED_CHANNELS``: downstream,
  utilities and wastewater capital, non-electric utilities, and the
  natural-gas intensity) are multiplied by the raw structural input
  ``struct_water_removal`` (0 when the broth itself is sold, route 'none').
  A hard physical constraint: without evaporation, distillation or
  wastewater treatment of the broth water, these titer-proportional costs
  vanish. Structural inputs are never fed to the encoder.

All raw head outputs are affine-scaled by per-target training means and
standard deviations (buffers), so that the network works on O(1) numbers.
``ModelConfig(kind='plain')`` is the plain-MLP baseline that predicts
``a', b', c', d', g`` directly (prices as ordinary features); the flags
``channels``, ``price_explicit``, ``structured_fs`` and ``g_bound`` ablate the
individual structural components.
"""
import math
from dataclasses import dataclass, asdict

import numpy as np
import torch
from torch import nn
import torch.nn.functional as F

try:
    from .physics import CHANNELS, PRICE_CHANNELS, INTENSITY_KEYS
except Exception:  # pragma: no cover
    CHANNELS = ('feedstock', 'electricity', 'natural_gas', 'neutralization',
                'materials_other', 'utilities_other', 'capex_upstream',
                'capex_fermentation', 'capex_downstream', 'capex_wastewater',
                'capex_utilities', 'fixed', 'coproducts', 'closure')
    PRICE_CHANNELS = {'feedstock': 'fs', 'electricity': 'el', 'natural_gas': 'ng'}
    INTENSITY_KEYS = ('fs', 'el', 'ng')

__all__ = ('ModelConfig', 'CoefficientNet', 'mpsp_equation', 'mpsp_equation_np',
           'make_layout', 'COEFF_KEYS', 'CHANNELS', 'PRICE_CHANNELS', 'INTENSITY_KEYS')

COEFF_KEYS = ('a', 'b', 'c', 'd')
#: channels / intensities whose c', d' are gated by the water-removal input
GATED_CHANNELS = ('capex_downstream', 'capex_utilities', 'capex_wastewater', 'utilities_other')
GATED_INTENSITIES = ('ng',)
GATED_COEFFS = (2, 3)
#: intensity key -> channel name
INTENSITY_CHANNEL = {v: k for k, v in PRICE_CHANNELS.items()}


def mpsp_equation(coef, g, y, t):
    """Equation layer (torch): ``coef`` (..., 4) broadcast against ``y``, ``t``
    and ``g``."""
    return (coef[..., 0] + coef[..., 1] / y + coef[..., 2] / t + coef[..., 3] / (y * t)) / g


def mpsp_equation_np(coef, g, y, t):
    """Numpy twin of :func:`mpsp_equation`."""
    coef = np.asarray(coef, float)
    return (coef[..., 0] + coef[..., 1] / y + coef[..., 2] / t + coef[..., 3] / (y * t)) / g


@dataclass
class ModelConfig:
    """Network and structural options."""
    kind: str = 'physics'        # 'physics' | 'plain'
    width: int = 128
    depth: int = 3
    dropout: float = 0.05
    channels: bool = True        # predict per-channel coefficients, totals = sum
    price_explicit: bool = True  # price x intensity channels (prices not in encoder)
    structured_fs: bool = True   # b'_fs = exp(z)/w_sugar
    fs_eta_bound: bool = True    # z >= 0  (eta <= 1)
    g_bound: bool = True         # g = g_max * (1 + g_margin) * sigmoid(.)
    g_margin: float = 0.0        # relative slack of the bound (set from data checks)
    head_init: float = 0.1       # scale of the initial head weights
    #: clip the standardized continuous encoder features to the training range
    #: widened by this margin (in training standard deviations); < 0 = off.
    #: Bounds the network's extrapolation in descriptor space (e.g. a new
    #: product whose descriptors lie outside the training products'); the
    #: applicability-domain report still flags such queries. Productivity and
    #: the price-explicit prices (not encoder inputs) are never clipped.
    clip_inputs: float = -1.0
    #: multiply c', d' of GATED_CHANNELS / GATED_INTENSITIES by the structural
    #: input ``struct_water_removal`` (needs FeatureSpec(downstream=True);
    #: without that column it has no effect)
    water_gate: bool = False

    def __post_init__(self):
        if self.kind not in ('physics', 'plain'):
            raise ValueError(f'unknown model kind {self.kind!r}')
        if self.kind == 'plain':
            self.channels = self.price_explicit = self.structured_fs = False
            self.fs_eta_bound = self.g_bound = False
        if not self.price_explicit:
            self.structured_fs = False

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        d = dict(d)
        # checkpoints saved before input clipping existed were trained without it
        d.setdefault('clip_inputs', -1.0)
        d.setdefault('water_gate', False)
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


def make_layout(cfg, columns):
    """Which channels / intensities the model predicts given the available
    target columns: ``{'channels': [...], 'price_keys': [...]}``. Price-explicit
    channels become intensities when ``cfg.price_explicit`` and their
    ``int__<key>__*`` targets exist; all other channels with ``ch__`` targets
    are predicted directly (only if ``cfg.channels``)."""
    columns = set(columns)
    price_keys = []
    if cfg.price_explicit:
        price_keys = [k for k in INTENSITY_KEYS
                      if all(f'int__{k}__{q}' in columns for q in COEFF_KEYS)]
    channels = []
    if cfg.channels:
        skip = {INTENSITY_CHANNEL[k] for k in price_keys}
        extra = sorted({c.split('__')[1] for c in columns
                        if c.startswith('ch__') and c.count('__') == 2} - set(CHANNELS))
        channels = [c for c in list(CHANNELS) + extra if c not in skip
                    and all(f'ch__{c}__{q}' in columns for q in COEFF_KEYS)]
    return {'channels': channels, 'price_keys': price_keys}


def _inv_softplus(x):
    x = max(float(x), 1e-4)
    return x + math.log(-math.expm1(-x))


class CoefficientNet(nn.Module):
    """Physics-informed coefficient network (see module docstring).

    Parameters
    ----------
    n_features : int
        Number of encoder features (FeatureSpec.n_features).
    cfg : ModelConfig
    layout : dict
        Output of :func:`make_layout`.
    scales : dict
        Target scales (``train.fit_target_scales``): ``coef_mean/std`` (4),
        ``g_mean/std``, ``ch_mean/std`` (C x 4), ``int_mean/std`` (J x 4),
        ``rest_mean/std`` (4), ``fs_logb_center/scale``, ``g_logit0``.
    r_stats : dict
        Standardization of log r and 1/r (FeatureSpec.r_stats).
    price_index : list of int
        Columns of the ``prices`` input for each of ``layout['price_keys']``.
    enc_index : list of int
        Columns of ``X`` fed to the encoder (default: all;
        ``FeatureSpec.encoder_index`` leaves out the raw structural inputs).
    gate_index : int
        Column of ``X`` holding ``struct_water_removal`` (used if
        ``cfg.water_gate``).
    """

    def __init__(self, n_features, cfg, layout, scales, r_stats, price_index=None, enc_index=None,
                 gate_index=None):
        super().__init__()
        self.cfg = cfg
        self.layout = {k: list(v) for k, v in layout.items()}
        self.n_features = int(n_features)
        C, J = len(layout['channels']), len(layout['price_keys'])
        self.n_channels, self.n_price = C, J
        self.price_index = list(price_index if price_index is not None
                                else [INTENSITY_KEYS.index(k) for k in layout['price_keys']])
        enc = list(range(self.n_features)) if enc_index is None else [int(i) for i in enc_index]
        self.enc_all = enc == list(range(self.n_features))
        # derived from the FeatureSpec: not part of the state_dict (old checkpoints load)
        self.register_buffer('enc_index', torch.as_tensor(enc, dtype=torch.long), persistent=False)
        self.gate_index = None if (gate_index is None or not cfg.water_gate or cfg.kind == 'plain') \
            else int(gate_index)
        layers, d_in = [], len(enc) + 2
        for _ in range(cfg.depth):
            layers += [nn.Linear(d_in, cfg.width), nn.LayerNorm(cfg.width), nn.SiLU(),
                       nn.Dropout(cfg.dropout)]
            d_in = cfg.width
        self.encoder = nn.Sequential(*layers)
        w = cfg.width if cfg.depth > 0 else d_in

        def buf(name, value, shape=None):
            t = torch.as_tensor(np.asarray(value, float), dtype=torch.float32)
            if shape is not None:
                t = t.reshape(shape)
            self.register_buffer(name, t.clone())

        if self.gate_index is not None:
            gc = np.ones((C, 4))
            for i, c in enumerate(layout['channels']):
                if c in GATED_CHANNELS:
                    gc[i, list(GATED_COEFFS)] = 0.
            gi = np.ones((J, 4))
            for i, k in enumerate(layout['price_keys']):
                if k in GATED_INTENSITIES:
                    gi[i, list(GATED_COEFFS)] = 0.
            # 1 = ungated entry; 0 = entry multiplied by the gate (derived, not saved)
            self.register_buffer('gate_ch', torch.as_tensor(gc, dtype=torch.float32), persistent=False)
            self.register_buffer('gate_int', torch.as_tensor(gi, dtype=torch.float32), persistent=False)
        rs = r_stats
        if cfg.clip_inputs >= 0:
            lo = np.asarray(scales.get('x_lo', np.full(self.n_features, -np.inf)), float)
            hi = np.asarray(scales.get('x_hi', np.full(self.n_features, np.inf)), float)
            buf('x_lo', np.where(np.isfinite(lo), lo - cfg.clip_inputs, -1e30))
            buf('x_hi', np.where(np.isfinite(hi), hi + cfg.clip_inputs, 1e30))
        buf('r_mu', [rs['log_r_mean'], rs['inv_r_mean']])
        buf('r_sd', [rs['log_r_std'], rs['inv_r_std']])
        if cfg.kind == 'plain':
            self.head = nn.Linear(w, 5)
            buf('out_mean', list(scales['coef_mean']) + [scales['g_mean']])
            buf('out_std', list(scales['coef_std']) + [scales['g_std']])
        else:
            if J:
                self.int_head = nn.Linear(w, 4 * J)
                buf('int_mean', scales['int_mean'], (J, 4))
                buf('int_std', scales['int_std'], (J, 4))
                self.j_fs = layout['price_keys'].index('fs') if 'fs' in layout['price_keys'] else -1
                c, s = float(scales.get('fs_logb_center', 0.)), float(scales.get('fs_logb_scale', 0.1))
                if cfg.fs_eta_bound:
                    # u -> softplus(u) = z; linearize scale at the center
                    c = _inv_softplus(c)
                    s = s / max(1. - math.exp(-max(float(scales.get('fs_logb_center', 0.)), 1e-4)), 1e-2)
                buf('fs_u_center', c)
                buf('fs_u_scale', s)
            if C:
                self.ch_head = nn.Linear(w, 4 * C)
                buf('ch_mean', scales['ch_mean'], (C, 4))
                buf('ch_std', scales['ch_std'], (C, 4))
            else:
                self.rest_head = nn.Linear(w, 4)
                buf('rest_mean', scales['rest_mean'])
                buf('rest_std', scales['rest_std'])
            self.g_head = nn.Linear(w, 1)
            if cfg.g_bound:
                buf('g_logit0', scales['g_logit0'])
            else:
                buf('g_mean', scales['g_mean'])
                buf('g_std', scales['g_std'])
        for m in self.modules():
            if isinstance(m, nn.Linear) and m is not None and not any(m is l for l in self.encoder):
                with torch.no_grad():
                    m.weight.mul_(cfg.head_init)
                    m.bias.zero_()

    # ---- inputs ----
    def r_features(self, log_r):
        return (torch.stack([log_r, torch.exp(-log_r)], -1) - self.r_mu) / self.r_sd

    def encode(self, X, log_r):
        if self.cfg.clip_inputs >= 0:
            X = torch.maximum(torch.minimum(X, self.x_hi.to(X.dtype)), self.x_lo.to(X.dtype))
        if not self.enc_all:
            X = X[:, self.enc_index]
        return self.encoder(torch.cat([X, self.r_features(log_r)], -1))

    def _gate(self, X, ungated):
        """Multiplier (n, K, 4): 1 where ``ungated``, else the gate input."""
        gv = X[:, self.gate_index].clamp(0., 1.)[:, None, None]
        u = ungated.to(X.dtype)[None]
        return u + (1. - u) * gv

    # ---- forward ----
    def forward(self, X, log_r, prices, w_sugar, g_max):
        """Return a dict with ``coef`` (n, 4) totals a', b', c', d' (primed
        basis), ``g`` (n,), and, depending on the configuration,
        ``channels`` (n, C, 4), ``intensities`` (n, J, 4),
        ``price_channels`` (n, J, 4) = price x intensity and, without channel
        heads, ``rest`` (n, 4) = all non-price costs."""
        h = self.encode(X, log_r)
        n = X.shape[0]
        if self.cfg.kind == 'plain':
            o = self.out_mean + self.out_std * self.head(h)
            return {'coef': o[:, :4], 'g': o[:, 4]}
        out = {}
        total = 0.
        if self.n_price:
            raw = self.int_head(h).view(n, self.n_price, 4)
            q = self.int_mean + self.int_std * raw
            if self.cfg.structured_fs and self.j_fs >= 0:
                u = self.fs_u_center + self.fs_u_scale * raw[:, self.j_fs, 1]
                z = F.softplus(u) if self.cfg.fs_eta_bound else u
                qb = torch.exp(z) / w_sugar
                qfs = torch.stack([q[:, self.j_fs, 0], qb, q[:, self.j_fs, 2], q[:, self.j_fs, 3]], -1)
                q = torch.cat([q[:, :self.j_fs], qfs[:, None], q[:, self.j_fs + 1:]], 1)
            if self.gate_index is not None:
                q = q * self._gate(X, self.gate_int)
            p = prices[:, self.price_index]
            pc = p[:, :, None] * q
            out['intensities'] = q
            out['price_channels'] = pc
            total = total + pc.sum(1)
        if self.n_channels:
            ch = self.ch_mean + self.ch_std * self.ch_head(h).view(n, self.n_channels, 4)
            if self.gate_index is not None:
                ch = ch * self._gate(X, self.gate_ch)
            out['channels'] = ch
            total = total + ch.sum(1)
        else:
            # all non-price costs in one head; returned so that the reported
            # channel breakdown still sums to the totals
            rest = self.rest_mean + self.rest_std * self.rest_head(h)
            out['rest'] = rest
            total = total + rest
        zg = self.g_head(h)[:, 0]
        if self.cfg.g_bound:
            g = g_max * (1. + self.cfg.g_margin) * torch.sigmoid(self.g_logit0 + zg)
        else:
            g = self.g_mean + self.g_std * zg
        out['coef'] = total
        out['g'] = g
        return out

    @staticmethod
    def mpsp(out, y, t):
        """MPSP at points ``y``, ``t`` of shape (n, P) for model output ``out``."""
        return mpsp_equation(out['coef'][:, None, :], out['g'][:, None], y, t)

    def channel_names(self):
        """Names of all channels in output order (price-explicit first; 'rest'
        when the other channels are not resolved, ``channels=False``)."""
        rest = ['rest'] if self.cfg.kind != 'plain' and not self.n_channels else []
        return [INTENSITY_CHANNEL[k] for k in self.layout['price_keys']] + list(self.layout['channels']) + rest

    def describe(self):
        n = sum(p.numel() for p in self.parameters())
        return (f'CoefficientNet(kind={self.cfg.kind}, width={self.cfg.width}, depth={self.cfg.depth}, '
                f'channels={self.layout["channels"]}, price_keys={self.layout["price_keys"]}, '
                f'params={n})')
