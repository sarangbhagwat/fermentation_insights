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
  ``clip_groups`` clips only (or also) the columns of the given feature groups,
  e.g. the feedstock composition, whatever ``clip_inputs``;
* optional extrapolation safeguards (off by default; need a FeatureSpec with
  ``downstream=True``): ``fs_prior`` writes the feedstock efficiency as
  ``z = log(1/s_avail) + dz`` with the composition's sugar availability
  ``s_avail`` (structural input ``struct_sugar_availability``), so the network
  predicts only the residual ``dz``; ``neutralization_gate`` multiplies the
  neutralization channel by ``struct_neutralization`` (0 when there is
  nothing to neutralize: no acid group, no neutralizing agent);
* optional power gate (``power_gate``; needs ``FeatureSpec(downstream=True,
  power_export=True)``): where electricity is not among the coproducts
  (structural input ``struct_power_export`` = 0: no turbogenerator, so no
  power export) every coefficient of the net electricity intensity is
  mapped to ``s softplus(q / s) >= 0`` (``s`` = ``POWER_GATE_SCALE`` x its
  training std), so the purchased electricity and dMPSP/dp_el are >= 0 at
  every y, t > 0; it acts after the productivity structure (a monotone map
  of a decreasing function of r keeps dMPSP/dr <= 0) and before the price
  (exact price linearity);
* optional water-removal gate (``water_gate``; needs a FeatureSpec with
  ``downstream=True``): the c' and d' coefficients of the channels that scale
  with the broth water processed downstream (``GATED_CHANNELS``: downstream,
  utilities and wastewater capital, non-electric utilities, and the
  natural-gas intensity) are multiplied by the raw structural input
  ``struct_water_removal`` (0 when the broth itself is sold, route 'none').
  A hard physical constraint: without evaporation, distillation or
  wastewater treatment of the broth water, these titer-proportional costs
  vanish. Structural inputs are never fed to the encoder;
* optional channel-wise input masking (``channel_inputs``; off by default):
  instead of one shared encoder, every distinct set of feature groups
  (``features.feature_group``: feedstock, scale, financial, fermentation,
  neutralization, product, final product, upgrading, separation operations,
  water-removal route, coproducts, utilities, and the pseudo-group
  'productivity' = log r, 1/r) gets its own sub-network, and each output
  (intensity, channel, recovery) is computed only from the sub-network of its
  groups (``CHANNEL_INPUT_PRESETS['physics']``: e.g. the feedstock intensity
  and upstream capital see no product information, the recovery no feedstock
  information), so an excluded input has exactly zero influence on it. The
  hybrid (``residual_width > 0``) adds a shared residual network on all
  inputs (zero-initialized head, strong weight decay
  ``residual_weight_decay``).
* optional hard productivity structure (``r_structure``; off by default):
  the four coefficients of every listed output (``R_STRUCTURE_PRESETS
  ['physics']``: fermentation capital and the electricity intensity) are
  ``theta0(x) + theta1(x) * (r_ref / r)**beta`` with ``theta1 >= 0``
  (softplus), ``beta`` fixed (``r_beta``) or learned within
  ``r_beta_bounds``, and ``theta0``, ``theta1`` computed without the
  productivity inputs (masked model: 'productivity' removed from those
  outputs' groups; shared encoder: no productivity inputs at all), so that
  ``dMPSP/dr <= 0`` holds by construction (fermenter volume ~ product rate /
  r; agitation / aeration power ~ volume). ``r_ref`` = exp(training mean of
  log r) (``FeatureSpec.r_stats``).

All raw head outputs are affine-scaled by per-target training means and
standard deviations (buffers), so that the network works on O(1) numbers.
``ModelConfig(kind='plain')`` is the plain-MLP baseline that predicts
``a', b', c', d', g`` directly (prices as ordinary features); the flags
``channels``, ``price_explicit``, ``structured_fs`` and ``g_bound`` ablate the
individual structural components.
"""
import math
import warnings
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
           'make_layout', 'COEFF_KEYS', 'CHANNELS', 'PRICE_CHANNELS', 'INTENSITY_KEYS',
           'CHANNEL_INPUT_PRESETS', 'resolve_channel_inputs', 'R_STRUCTURE_PRESETS', 'resolve_r_structure')

COEFF_KEYS = ('a', 'b', 'c', 'd')
#: channels / intensities whose c', d' are gated by the water-removal input
GATED_CHANNELS = ('capex_downstream', 'capex_utilities', 'capex_wastewater', 'utilities_other')
GATED_INTENSITIES = ('ng',)
#: smoothing scale of the power gate (fraction of the training std of each
#: electricity-intensity coefficient): s softplus(q / s) ~ q for q >> s
POWER_GATE_SCALE = 0.1
GATED_COEFFS = (2, 3)
#: intensity key -> channel name
INTENSITY_CHANNEL = {v: k for k, v in PRICE_CHANNELS.items()}

#: pseudo feature group of the channel-wise input masking: the productivity
#: features (log r, 1/r)
R_GROUP = 'productivity'
#: output key of the recovery in ``channel_inputs``
G_KEY = 'g'
#: Channel-wise input masking presets (``ModelConfig.channel_inputs``):
#: output -> feature groups (``features.FEATURE_GROUPS`` + ``R_GROUP``) it may
#: depend on. Outputs are channel names (a price-explicit channel by its
#: channel name: its intensity head), and ``G_KEY`` for the recovery; an
#: output not listed sees every group. DESIGN.md section 3 gives the physics
#: rationale and the data checks. The product-identifying groups (product and
#: final-product chemistry, separation operations, neutralization chemistry,
#: fermentation regime, coproducts) reach only the channels whose process
#: area handles the product.
CHANNEL_INPUT_PRESETS = {
    'physics': {
        # purchase of sugar: composition and sugar release only
        'feedstock': ('feedstock', 'pretreatment'),
        # feedstock handling, pretreatment, saccharification: no product
        'capex_upstream': ('feedstock', 'pretreatment', 'scale', 'financial'),
        # labor: plant scale
        'fixed': ('feedstock', 'scale', 'financial'),
        'closure': ('feedstock', 'scale', 'financial'),
        # fermenter volume ~ product rate / r, broth volume ~ 1/t; aeration, pH control
        'capex_fermentation': ('fermentation', R_GROUP, 'neutralization', 'reagent_prices', 'feedstock', 'scale',
                               'financial'),
        # base + acid stoichiometry and their prices
        'neutralization': ('neutralization', 'reagent_prices'),
        # product recovery and upgrading
        'capex_downstream': ('product', 'final_product', 'upgrading', 'separation', 'water_route',
                             'neutralization', 'reagent_prices', 'scale', 'financial'),
        'utilities_other': ('product', 'final_product', 'upgrading', 'separation', 'water_route', 'fermentation',
                            'utilities', 'scale'),
        'natural_gas': ('product', 'final_product', 'upgrading', 'separation', 'water_route', 'feedstock',
                        'utilities', 'scale'),
        'materials_other': ('feedstock', 'pretreatment', 'fermentation', 'product', 'final_product', 'upgrading',
                            'separation', 'water_route', 'scale'),
        # broth water, residual solids and organics; boiler / turbogenerator
        'capex_wastewater': ('feedstock', 'fermentation', 'water_route', 'coproducts', 'scale', 'financial'),
        'capex_utilities': ('feedstock', 'water_route', 'upgrading', 'coproducts', 'utilities', 'scale',
                            'financial'),
        'electricity': ('feedstock', 'fermentation', R_GROUP, 'water_route', 'upgrading', 'coproducts',
                        'utilities', 'scale'),
        'coproducts': ('feedstock', 'fermentation', 'coproducts', 'scale'),
        # recovery: separation of the product, never the feedstock
        G_KEY: ('product', 'final_product', 'upgrading', 'separation', 'water_route', 'neutralization',
                'reagent_prices'),
    },
}
#: 'physics' with the recovery blind to the reagent prices (they set what the
#: base and acid cost, not how much product is recovered; DESIGN.md section 6,
#: round 4). A new name, so that checkpoints saved with 'physics' load unchanged.
CHANNEL_INPUT_PRESETS['physics_v2'] = {
    **CHANNEL_INPUT_PRESETS['physics'],
    G_KEY: tuple(g for g in CHANNEL_INPUT_PRESETS['physics'][G_KEY] if g != 'reagent_prices'),
}


#: Hard productivity structure presets (``ModelConfig.r_structure``): outputs
#: (channel names; a price-explicit channel by its channel name = its
#: intensity) whose coefficients are theta0(x) + theta1(x) (r_ref / r)^beta.
#: Development data (DESIGN.md section 3): over the productivity samples of a
#: configuration only the fermentation capital changes materially (13 % of
#: the total numerator, up to 41 %; exactly ~ 1/r: beta = 1.0, R^2 = 0.9999,
#: theta1 >= 0 in every coefficient costs R^2 0.0004); the electricity
#: intensity changes by <= 1.4 % (agitation / aeration power ~ volume ~ 1/r);
#: every other channel by <= 1.2 % without a consistent sign (noise of the
#: per-sample fits).
R_STRUCTURE_PRESETS = {'physics': ('capex_fermentation', 'electricity')}


def resolve_r_structure(r_structure):
    """``ModelConfig.r_structure`` -> tuple of output names (() = off): a
    preset name of ``R_STRUCTURE_PRESETS`` or a list of output names."""
    if r_structure in (None, '', False) or (not isinstance(r_structure, str) and not len(r_structure)):
        return ()
    if isinstance(r_structure, str):
        if r_structure not in R_STRUCTURE_PRESETS:
            raise ValueError(f'unknown r_structure preset {r_structure!r} (known: {sorted(R_STRUCTURE_PRESETS)}); '
                             'pass a list of output names for a single output')
        return tuple(R_STRUCTURE_PRESETS[r_structure])
    return tuple(str(x) for x in r_structure)


def _known_group(g):
    try:
        from .features import FEATURE_GROUPS
    except Exception:  # pragma: no cover
        return True
    return g in FEATURE_GROUPS


def resolve_channel_inputs(channel_inputs):
    """``ModelConfig.channel_inputs`` -> {output: tuple of groups} (None = off).
    A preset name, a dict (output -> groups) or a dict with key ``'preset'``
    (a preset name) whose other entries override the preset."""
    if channel_inputs in (None, '', False):
        return None
    if isinstance(channel_inputs, str):
        if channel_inputs not in CHANNEL_INPUT_PRESETS:
            raise ValueError(f'unknown channel_inputs preset {channel_inputs!r} '
                             f'(known: {sorted(CHANNEL_INPUT_PRESETS)})')
        return {k: tuple(v) for k, v in CHANNEL_INPUT_PRESETS[channel_inputs].items()}
    d = dict(channel_inputs)
    base = resolve_channel_inputs(d.pop('preset')) if 'preset' in d else {}
    base.update({k: tuple(v) for k, v in d.items()})
    return base


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
    #: with structured_fs: the feedstock intensity is the b' term alone
    #: (a' = c' = d' = 0): kg feedstock per kg fermentation product =
    #: 1 / (eta w_sugar y), the sugar-balance identity at constant recovery.
    #: Its departures in the data come from the per-point recovery varying
    #: around the sample's median g (a pure b'/y fit is within 5 % RMS of
    #: every development sample's intensity but lactic's, 6-9 %) and from
    #: ill-conditioned OLS splits (TAL c' = -17 ... -90 with a pure-b' error
    #: of 4 %), which a product-blind feedstock head cannot learn anyway
    fs_b_only: bool = False
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
    #: feature groups (``features.FEATURE_GROUPS``) whose continuous columns
    #: are always clipped to their training range (margin
    #: ``max(clip_inputs, 0)``), also when ``clip_inputs < 0`` leaves the other
    #: columns unclipped; None = no group-specific clipping. E.g.
    #: ('feedstock', 'pretreatment'): a new feedstock's composition is held at
    #: the edge of the training feedstocks' range (minor mass fractions such
    #: as 'other' or ash have a tiny training spread, so a new feedstock
    #: lies many standard deviations outside it and saturates the feedstock
    #: sub-networks), while the product descriptors keep their linear
    #: extrapolation.
    clip_groups: object = None
    #: encoder columns never clipped by ``clip_inputs`` / ``clip_groups``
    #: (names of FeatureSpec.columns; None = none). For exact structural
    #: quantities, e.g. the final-product stoichiometric mass ratio: clipping
    #: a held-out product's ratio to the training range maps it onto the
    #: value of other products (round 5: HP's 0.80 to 1.00 when HP is held out)
    clip_exempt: object = None
    #: multiply c', d' of GATED_CHANNELS / GATED_INTENSITIES by the structural
    #: input ``struct_water_removal`` (needs FeatureSpec(downstream=True);
    #: without that column it has no effect)
    water_gate: bool = False
    #: multiply the whole neutralization channel by the structural input
    #: ``struct_neutralization`` (0 when the fermentation product has no
    #: acid group and no neutralizing / pH-control agent is used: nothing to
    #: neutralize or acidulate; needs FeatureSpec(downstream=True))
    neutralization_gate: bool = False
    #: non-negative net electricity intensity (each coefficient s softplus(q/s))
    #: where the structural input ``struct_power_export`` is 0 (electricity
    #: not among the coproducts: no on-site generation, nothing to export);
    #: needs FeatureSpec(downstream=True, power_export=True), which
    #: ``train.train_ensemble`` requests itself when this is on
    power_gate: bool = False
    #: structured feedstock intensity with a physics prior of the sugar
    #: availability (needs FeatureSpec(downstream=True)):
    #: z = log(1/eta) = log(1/s_avail) + dz, s_avail = structural input
    #: ``struct_sugar_availability`` (free sugars and starch fully, structural
    #: polysaccharides by the saccharification yield, none without
    #: pretreatment); the network predicts only the residual dz (centered on
    #: its training median; ``fs_eta_bound`` does not apply), so a new
    #: feedstock's efficiency follows its composition instead of the
    #: network's extrapolation
    fs_prior: bool = False
    #: channel-wise input masking (physics model only; None = one shared
    #: encoder): a preset name of CHANNEL_INPUT_PRESETS or a dict output ->
    #: feature groups. Each distinct group set gets its own sub-network fed
    #: only the encoder columns of those groups (FeatureSpec.column_groups)
    #: and, with the group 'productivity', log r and 1/r.
    channel_inputs: object = None
    masked_width: int = 64
    masked_depth: int = 2
    #: hybrid (with channel_inputs): a shared residual network of this width
    #: on all inputs whose standardized outputs are added to those of the
    #: masked sub-networks (0 = none); its weights get their own (strong)
    #: weight decay ``residual_weight_decay`` (train.train_model)
    residual_width: int = 0
    residual_depth: int = 2
    residual_weight_decay: float = 1.0
    #: hard productivity structure (physics model only; None = off): a preset
    #: name of R_STRUCTURE_PRESETS or a list of output names whose
    #: coefficients are theta0(x) + theta1(x) * (r_ref / r)^beta, theta1 >= 0,
    #: theta0 and theta1 without productivity inputs (dMPSP/dr <= 0 by
    #: construction when no other output sees the productivity)
    r_structure: object = None
    #: exponent beta (fixed), or its initial value when ``r_beta_learn``
    r_beta: float = 1.0
    #: learn one beta per structured output within ``r_beta_bounds``
    r_beta_learn: bool = False
    r_beta_bounds: tuple = (0.5, 1.5)
    #: initial theta1 at r_ref: the pooled within-configuration slope of the
    #: training coefficients on (r_ref / r) (target scales ``*_r_slope``), at
    #: least ``r_slope_floor`` x the target's training std; without slopes in
    #: the scales ``r_slope_init`` x std
    r_slope_init: float = 0.3
    r_slope_floor: float = 0.01

    def __post_init__(self):
        if self.kind not in ('physics', 'plain'):
            raise ValueError(f'unknown model kind {self.kind!r}')
        if self.kind == 'plain':
            self.channels = self.price_explicit = self.structured_fs = False
            self.fs_eta_bound = self.g_bound = False
            self.channel_inputs = None
            self.r_structure = None
            # physics-model structures (the pipeline merges them into every NN method)
            self.fs_prior = self.neutralization_gate = self.power_gate = False
            self.clip_groups = None
        if not self.price_explicit:
            self.structured_fs = False
        if isinstance(self.channel_inputs, dict):
            self.channel_inputs = {k: (v if k == 'preset' else list(v)) for k, v in self.channel_inputs.items()}
        resolve_channel_inputs(self.channel_inputs)   # validates
        if self.r_structure is not None and not isinstance(self.r_structure, str):
            self.r_structure = list(self.r_structure) or None
        resolve_r_structure(self.r_structure)          # validates
        self.r_beta_bounds = tuple(float(b) for b in self.r_beta_bounds)
        lo, hi = self.r_beta_bounds
        if not 0. < lo <= hi:
            raise ValueError(f'r_beta_bounds must satisfy 0 < lo <= hi, got {self.r_beta_bounds}')
        if not self.r_beta > 0.:
            raise ValueError(f'r_beta must be > 0 (dMPSP/dr <= 0), got {self.r_beta}')
        if self.r_beta_learn and not lo <= self.r_beta <= hi:
            raise ValueError(f'r_beta {self.r_beta} outside r_beta_bounds {self.r_beta_bounds}')
        if self.clip_groups in ((), [], ''):
            self.clip_groups = None
        if self.clip_groups is not None:
            if isinstance(self.clip_groups, str):
                self.clip_groups = [self.clip_groups]
            self.clip_groups = [str(g) for g in self.clip_groups]
            unknown = [g for g in self.clip_groups if not _known_group(g)]
            if unknown:
                raise ValueError(f'clip_groups: unknown feature group(s) {unknown}')
        if self.clip_exempt in ((), [], ''):
            self.clip_exempt = None
        if self.clip_exempt is not None:
            self.clip_exempt = [str(c) for c in ([self.clip_exempt] if isinstance(self.clip_exempt, str)
                                                 else self.clip_exempt)]

    @property
    def masked(self):
        return self.kind != 'plain' and resolve_channel_inputs(self.channel_inputs) is not None

    @property
    def r_structured(self):
        """Names of the r-structured outputs (() = off)."""
        return () if self.kind == 'plain' else resolve_r_structure(self.r_structure)

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        d = dict(d)
        # checkpoints saved before input clipping existed were trained without it
        d.setdefault('clip_inputs', -1.0)
        d.setdefault('water_gate', False)
        # ... and before channel-wise input masking: one shared encoder
        d.setdefault('channel_inputs', None)
        # ... and before the hard productivity structure: off
        d.setdefault('r_structure', None)
        # ... and before group-restricted clipping / the neutralization gate
        d.setdefault('clip_groups', None)
        d.setdefault('clip_exempt', None)
        d.setdefault('neutralization_gate', False)
        d.setdefault('fs_prior', False)
        d.setdefault('power_gate', False)
        d.setdefault('fs_b_only', False)
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
    input_groups : dict
        {feature group: [columns of X]} (``FeatureSpec.column_groups``);
        needed by the channel-wise input masking (``cfg.channel_inputs``)
        and by ``cfg.clip_groups``.
    neut_index : int
        Column of ``X`` holding ``struct_neutralization`` (used if
        ``cfg.neutralization_gate``).
    avail_index : int
        Column of ``X`` holding ``struct_sugar_availability`` (used if
        ``cfg.fs_prior``).
    power_index : int
        Column of ``X`` holding ``struct_power_export`` (used if
        ``cfg.power_gate``).
    """

    def __init__(self, n_features, cfg, layout, scales, r_stats, price_index=None, enc_index=None,
                 gate_index=None, input_groups=None, neut_index=None, avail_index=None, power_index=None,
                 clip_exempt_index=None):
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
        self.neut_index = None if (neut_index is None or not cfg.neutralization_gate or cfg.kind == 'plain'
                                   or 'neutralization' not in layout['channels']) else int(neut_index)
        self.j_neut = layout['channels'].index('neutralization') if self.neut_index is not None else -1
        self.avail_index = None if (avail_index is None or not cfg.fs_prior or not cfg.structured_fs
                                    or cfg.kind == 'plain') else int(avail_index)
        self.power_index = None if (power_index is None or not getattr(cfg, 'power_gate', False)
                                    or cfg.kind == 'plain' or 'el' not in layout['price_keys']) else int(power_index)
        self.j_el = layout['price_keys'].index('el') if self.power_index is not None else -1
        self.masked = cfg.masked
        self._setup_r_structure(layout)
        if self.masked:
            w = None
            body = self._build_masked(enc, input_groups)
        else:
            # hard productivity structure: the shared encoder sees no productivity
            layers, d_in = [], len(enc) + (0 if self.r_free else 2)
            for _ in range(cfg.depth):
                layers += [nn.Linear(d_in, cfg.width), nn.LayerNorm(cfg.width), nn.SiLU(),
                           nn.Dropout(cfg.dropout)]
                d_in = cfg.width
            self.encoder = nn.Sequential(*layers)
            body = set(map(id, self.encoder))
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
        self.clip_on = cfg.clip_inputs >= 0 or cfg.clip_groups is not None
        if self.clip_on:
            lo = np.asarray(scales.get('x_lo', np.full(self.n_features, -np.inf)), float)
            hi = np.asarray(scales.get('x_hi', np.full(self.n_features, np.inf)), float)
            margin = max(float(cfg.clip_inputs), 0.)
            if cfg.clip_inputs < 0:
                # only the columns of clip_groups
                if input_groups is None:
                    raise ValueError('clip_groups needs input_groups (FeatureSpec.column_groups)')
                sel = np.zeros(self.n_features, bool)
                for g in cfg.clip_groups:
                    sel[[int(i) for i in input_groups.get(g, [])]] = True
                lo, hi = np.where(sel, lo, -np.inf), np.where(sel, hi, np.inf)
            if clip_exempt_index:
                lo, hi = lo.copy(), hi.copy()
                lo[[int(i) for i in clip_exempt_index]] = -np.inf
                hi[[int(i) for i in clip_exempt_index]] = np.inf
            buf('x_lo', np.where(np.isfinite(lo), lo - margin, -1e30))
            buf('x_hi', np.where(np.isfinite(hi), hi + margin, 1e30))
        buf('r_mu', [rs['log_r_mean'], rs['inv_r_mean']])
        buf('r_sd', [rs['log_r_std'], rs['inv_r_std']])
        if cfg.kind == 'plain':
            self.head = nn.Linear(w, 5)
            buf('out_mean', list(scales['coef_mean']) + [scales['g_mean']])
            buf('out_std', list(scales['coef_std']) + [scales['g_std']])
        else:
            if J:
                if not self.masked:
                    self.int_head = nn.Linear(w, 4 * J)
                buf('int_mean', scales['int_mean'], (J, 4))
                buf('int_std', scales['int_std'], (J, 4))
                self.j_fs = layout['price_keys'].index('fs') if 'fs' in layout['price_keys'] else -1
                c, s = float(scales.get('fs_logb_center', 0.)), float(scales.get('fs_logb_scale', 0.1))
                if self.avail_index is not None:
                    # residual dz = z - log(1/s_avail), unbounded
                    c = float(scales.get('fs_prior_center', 0.))
                    s = float(scales.get('fs_prior_scale', 0.1))
                elif cfg.fs_eta_bound:
                    # u -> softplus(u) = z; linearize scale at the center
                    c = _inv_softplus(c)
                    s = s / max(1. - math.exp(-max(float(scales.get('fs_logb_center', 0.)), 1e-4)), 1e-2)
                buf('fs_u_center', c)
                buf('fs_u_scale', s)
                if self.power_index is not None:
                    sd = np.asarray(scales['int_std'], float).reshape(J, 4)[self.j_el]
                    # derived from the target scales: not saved
                    self.register_buffer('power_scale', torch.as_tensor(
                        POWER_GATE_SCALE * np.maximum(sd, 1e-8), dtype=torch.float32), persistent=False)
            if C:
                if not self.masked:
                    self.ch_head = nn.Linear(w, 4 * C)
                buf('ch_mean', scales['ch_mean'], (C, 4))
                buf('ch_std', scales['ch_std'], (C, 4))
            else:
                if not self.masked:
                    self.rest_head = nn.Linear(w, 4)
                buf('rest_mean', scales['rest_mean'])
                buf('rest_std', scales['rest_std'])
            if not self.masked:
                self.g_head = nn.Linear(w, 1)
                if self.n_rs:
                    self.r1_head = nn.Linear(w, 4 * self.n_rs)
            if self.n_rs:
                self._build_r_buffers(scales, C, J)
            if cfg.g_bound:
                buf('g_logit0', scales['g_logit0'])
            else:
                buf('g_mean', scales['g_mean'])
                buf('g_std', scales['g_std'])
        for m in self.modules():
            if isinstance(m, nn.Linear) and m is not None and id(m) not in body:
                with torch.no_grad():
                    m.weight.mul_(cfg.head_init)
                    m.bias.zero_()
        if self.masked and cfg.residual_width > 0:
            # the hybrid starts as the purely masked model
            with torch.no_grad():
                self.resid_head.weight.zero_()

    # ---- hard productivity structure ----
    def _setup_r_structure(self, layout):
        """Structured output slots (``rs_slots``, indices into
        ``_output_names()[:-1]``) of ``cfg.r_structure``. Without resolved
        channels (``channels=False``) a listed non-price channel structures
        the 'rest' head."""
        cfg = self.cfg
        rs = set(cfg.r_structured)
        names = self._output_names()[:-1] if cfg.kind != 'plain' else []
        if rs and 'rest' in names and any(n not in PRICE_CHANNELS for n in rs):
            rs.add('rest')
        self.rs_slots = [s for s, nm in enumerate(names) if nm in rs]
        self.rs_names = [names[s] for s in self.rs_slots]
        self.n_rs = len(self.rs_slots)
        # no productivity input anywhere but through the structure (shared encoder / residual network)
        self.r_free = self.n_rs > 0
        if rs and not self.n_rs:
            warnings.warn(f'r_structure: none of {sorted(rs)} is an output of this model; the structure is off')
        if self.n_rs and cfg.structured_fs and 'fs' in layout['price_keys'] and 'feedstock' in self.rs_names:
            raise ValueError('r_structure: the feedstock intensity is structured (structured_fs); '
                             'it cannot also be r-structured')

    def _build_r_buffers(self, scales, C, J):
        cfg = self.cfg
        S = self.n_rs
        scale = []
        E_int, E_ch = np.zeros((J, S)), np.zeros((max(C, 1), S))

        def pick(std_key, slope_key, shape, row):
            sd = np.asarray(scales[std_key], float).reshape(shape)[row]
            sl = scales.get(slope_key)
            if sl is None or not len(sl):
                return cfg.r_slope_init * sd
            sl = np.asarray(sl, float).reshape(shape)[row]
            return np.maximum(np.nan_to_num(sl), cfg.r_slope_floor * sd)

        for i, s in enumerate(self.rs_slots):
            if s < J:
                scale.append(pick('int_std', 'int_r_slope', (J, 4), s))
                E_int[s, i] = 1.
            elif C:
                scale.append(pick('ch_std', 'ch_r_slope', (C, 4), s - J))
                E_ch[s - J, i] = 1.
            else:
                scale.append(pick('rest_std', 'rest_r_slope', (1, 4), 0))
                E_ch[0, i] = 1.
        # derived from the config and the target scales: not saved
        # (theta1 = rs_std * softplus(raw1 + u0) = rs_std at initialization)
        for name, v in (('rs_std', np.array(scale)), ('rs_E_int', E_int), ('rs_E_ch', E_ch)):
            self.register_buffer(name, torch.as_tensor(v, dtype=torch.float32), persistent=False)
        self.rs_u0 = _inv_softplus(1.)
        lo, hi = cfg.r_beta_bounds
        if cfg.r_beta_learn:
            f = min(max((cfg.r_beta - lo) / (hi - lo), 1e-4), 1 - 1e-4) if hi > lo else 0.5
            self.r_beta_raw = nn.Parameter(torch.full((S,), math.log(f / (1 - f))))
        else:
            self.register_buffer('r_beta_fixed', torch.full((S,), float(cfg.r_beta)), persistent=False)

    def r_beta_values(self):
        """Exponent beta of every structured output, (n_rs,) tensor."""
        if self.cfg.r_beta_learn:
            lo, hi = self.cfg.r_beta_bounds
            return lo + (hi - lo) * torch.sigmoid(self.r_beta_raw)
        return self.r_beta_fixed

    @property
    def r_ref(self):
        """Reference productivity of theta1 (g/L/h): exp(training mean of log r)."""
        return float(torch.exp(self.r_mu[0]))

    def _r_terms(self, raw1, log_r):
        """(theta1 (n, S, 4) >= 0 at r_ref, additive term (n, S, 4)): the
        structured part theta1 * ((r_ref / r)^beta - 1); the '- 1' keeps
        the raw head of the output equal to the coefficient at r_ref."""
        th1 = self.rs_std * F.softplus(raw1 + self.rs_u0)
        beta = self.r_beta_values()
        rho = torch.exp(beta[None, :] * (self.r_mu[0] - log_r[:, None])) - 1.
        return th1, th1 * rho[:, :, None]

    # ---- channel-wise input masking ----
    def _output_names(self):
        """Names of the standardized output slots, in the order of ``_raw``:
        intensities (by channel name), channels (or 'rest'), and G_KEY."""
        names = [INTENSITY_CHANNEL[k] for k in self.layout['price_keys']]
        names += list(self.layout['channels']) if self.n_channels else ['rest']
        return names + [G_KEY]

    def _block(self, d_in, width, depth):
        layers = []
        for _ in range(depth):
            layers += [nn.Linear(d_in, width), nn.LayerNorm(width), nn.SiLU(), nn.Dropout(self.cfg.dropout)]
            d_in = width
        return nn.Sequential(*layers), (width if depth > 0 else d_in)

    def _build_masked(self, enc, input_groups):
        """Sub-networks of the channel-wise input masking (one per distinct
        group set) and, for the hybrid, the shared residual network. Returns
        the ids of the hidden (non-head) layers."""
        cfg = self.cfg
        if self.cfg.kind == 'plain':
            raise ValueError('channel_inputs needs the physics model')
        if input_groups is None:
            raise ValueError('channel_inputs needs input_groups (FeatureSpec.column_groups)')
        enc_set = set(enc)
        groups = {str(g): sorted(int(i) for i in idx if int(i) in enc_set) for g, idx in input_groups.items()}
        every = tuple(sorted(groups)) + (R_GROUP,)
        spec = resolve_channel_inputs(cfg.channel_inputs)
        names = self._output_names()
        sets, slots = [], []
        n_int = len(self.layout['price_keys'])
        for s, name in enumerate(names):
            gs = set(spec.get(name, every))
            if s in self.rs_slots:
                # hard productivity structure: theta0, theta1 without productivity inputs
                gs.discard(R_GROUP)
            if s >= n_int and name in PRICE_CHANNELS:
                # a price channel resolved as an ordinary channel (price_explicit=False):
                # its price is an encoder feature
                gs.add('prices')
            gs = tuple(sorted(gs))
            unknown = [g for g in gs if g != R_GROUP and g not in groups and not _known_group(g)]
            if unknown:
                raise ValueError(f'channel_inputs[{name!r}]: unknown feature group(s) {unknown}')
            if gs not in sets:
                sets.append(gs)
                slots.append([])
            slots[sets.index(gs)].append(s)
        self.n_slots = len(names) - 1
        self.sub_groups, self.sub_slots, self.sub_r = sets, slots, [R_GROUP in gs for gs in sets]
        nets, heads = [], []
        self.output_inputs = {}
        for k, gs in enumerate(sets):
            idx = sorted({i for g in gs for i in groups.get(g, [])})
            # derived from the FeatureSpec: not saved
            self.register_buffer(f'sub_index_{k}', torch.as_tensor(idx, dtype=torch.long), persistent=False)
            d_in = len(idx) + 2 * self.sub_r[k]
            net, w = self._block(max(d_in, 1), cfg.masked_width, cfg.masked_depth)
            nets.append(net)
            n_out = sum(1 if s == self.n_slots else 8 if s in self.rs_slots else 4 for s in slots[k])
            heads.append(nn.Linear(w, n_out))
            for s in slots[k]:
                self.output_inputs[names[s]] = dict(columns=idx, productivity=self.sub_r[k], groups=list(gs))
        self.sub_nets, self.sub_heads = nn.ModuleList(nets), nn.ModuleList(heads)
        body = {id(m) for net in nets for m in net}
        if cfg.residual_width > 0:
            self.resid_encoder, w = self._block(len(enc) + (0 if self.r_free else 2), cfg.residual_width,
                                                cfg.residual_depth)
            self.resid_head = nn.Linear(w, 4 * self.n_slots + 1)
            body |= {id(m) for m in self.resid_encoder}
        used = {i for v in self.output_inputs.values() for i in v['columns']}
        self.unused_inputs = sorted(enc_set - used)
        if self.unused_inputs and cfg.residual_width <= 0:
            warnings.warn(f'channel_inputs: {len(self.unused_inputs)} encoder input(s) belong to no group of any '
                          f'output and are ignored (e.g. columns {self.unused_inputs[:5]}; '
                          f'groups: {sorted(g for g in groups if not any(g in gs for gs in sets))})')
        return body

    def residual_parameters(self):
        """Parameters of the hybrid's shared residual network (train.py gives
        them ``cfg.residual_weight_decay``)."""
        if not (self.masked and self.cfg.residual_width > 0):
            return []
        return list(self.resid_encoder.parameters()) + list(self.resid_head.parameters())

    def _raw_masked(self, X, log_r):
        """Standardized outputs of the masked model: (n, S, 4) for the slots of
        ``_output_names`` but the last, (n,) for the recovery and, with the
        hard productivity structure, (n, n_rs, 4) raw theta1 heads (else None)."""
        n = X.shape[0]
        if self.clip_on:
            X = torch.maximum(torch.minimum(X, self.x_hi.to(X.dtype)), self.x_lo.to(X.dtype))
        rf = self.r_features(log_r)
        out, zg = [None] * self.n_slots, None
        out1 = [None] * self.n_rs
        for k, (net, head) in enumerate(zip(self.sub_nets, self.sub_heads)):
            idx = getattr(self, f'sub_index_{k}')
            parts = ([X[:, idx]] if len(idx) else []) + ([rf] if self.sub_r[k] else [])
            o = head(net(torch.cat(parts, -1) if parts else X.new_ones(n, 1)))
            j = 0
            for s in self.sub_slots[k]:
                if s == self.n_slots:
                    zg, j = o[:, j], j + 1
                else:
                    out[s], j = o[:, j:j + 4], j + 4
                    if s in self.rs_slots:
                        out1[self.rs_slots.index(s)], j = o[:, j:j + 4], j + 4
        raw = torch.stack(out, 1)
        if self.cfg.residual_width > 0:
            Xe = X if self.enc_all else X[:, self.enc_index]
            r = self.resid_head(self.resid_encoder(torch.cat([Xe] if self.r_free else [Xe, rf], -1)))
            raw = raw + r[:, :4 * self.n_slots].view(n, self.n_slots, 4)
            zg = zg + r[:, 4 * self.n_slots]
        return raw, zg, (torch.stack(out1, 1) if self.n_rs else None)

    # ---- inputs ----
    def r_features(self, log_r):
        return (torch.stack([log_r, torch.exp(-log_r)], -1) - self.r_mu) / self.r_sd

    def encode(self, X, log_r):
        if self.clip_on:
            X = torch.maximum(torch.minimum(X, self.x_hi.to(X.dtype)), self.x_lo.to(X.dtype))
        if not self.enc_all:
            X = X[:, self.enc_index]
        if self.r_free:
            return self.encoder(X)
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
        n = X.shape[0]
        raw1 = None
        if self.masked:
            h = None
            raw_all, zg, raw1 = self._raw_masked(X, log_r)
            J = self.n_price
            raw_int, raw_ch = raw_all[:, :J], raw_all[:, J:]
        else:
            h = self.encode(X, log_r)
        if self.cfg.kind == 'plain':
            o = self.out_mean + self.out_std * self.head(h)
            return {'coef': o[:, :4], 'g': o[:, 4]}
        out = {}
        total = 0.
        if self.n_rs:
            # hard productivity structure: theta(x, r) = base(x) + theta1(x) ((r_ref/r)^beta - 1)
            if raw1 is None:
                raw1 = self.r1_head(h).view(n, self.n_rs, 4)
            th1, radd = self._r_terms(raw1, log_r)
            base_rs = 0.
        if self.n_price:
            raw = raw_int if self.masked else self.int_head(h).view(n, self.n_price, 4)
            q = self.int_mean + self.int_std * raw
            if self.n_rs:
                base_rs = base_rs + torch.einsum('js,njk->nsk', self.rs_E_int, q)
                q = q + torch.einsum('js,nsk->njk', self.rs_E_int, radd)
            if self.cfg.structured_fs and self.j_fs >= 0:
                u = self.fs_u_center + self.fs_u_scale * raw[:, self.j_fs, 1]
                if self.avail_index is not None:
                    z = u - torch.log(X[:, self.avail_index].clamp(0.05, 1.))
                else:
                    z = F.softplus(u) if self.cfg.fs_eta_bound else u
                qb = torch.exp(z) / w_sugar
                if getattr(self.cfg, 'fs_b_only', False):
                    z0 = torch.zeros_like(qb)
                    qfs = torch.stack([z0, qb, z0, z0], -1)
                else:
                    qfs = torch.stack([q[:, self.j_fs, 0], qb, q[:, self.j_fs, 2], q[:, self.j_fs, 3]], -1)
                q = torch.cat([q[:, :self.j_fs], qfs[:, None], q[:, self.j_fs + 1:]], 1)
            if self.power_index is not None:
                # no power export: net electricity intensity >= 0 coefficient-wise
                e = X[:, self.power_index].clamp(0., 1.)[:, None]
                qe = q[:, self.j_el]
                sc_ = self.power_scale.to(qe.dtype)
                qe = e * qe + (1. - e) * sc_ * F.softplus(qe / sc_)
                q = torch.cat([q[:, :self.j_el], qe[:, None], q[:, self.j_el + 1:]], 1)
            if self.gate_index is not None:
                q = q * self._gate(X, self.gate_int)
            p = prices[:, self.price_index]
            pc = p[:, :, None] * q
            out['intensities'] = q
            out['price_channels'] = pc
            total = total + pc.sum(1)
        if self.n_channels:
            rc = raw_ch if self.masked else self.ch_head(h).view(n, self.n_channels, 4)
            ch = self.ch_mean + self.ch_std * rc
            if self.n_rs:
                base_rs = base_rs + torch.einsum('cs,nck->nsk', self.rs_E_ch, ch)
                ch = ch + torch.einsum('cs,nsk->nck', self.rs_E_ch, radd)
            if self.gate_index is not None:
                ch = ch * self._gate(X, self.gate_ch)
            if self.neut_index is not None:
                nv = X[:, self.neut_index].clamp(0., 1.)
                ch = torch.cat([ch[:, :self.j_neut], (ch[:, self.j_neut] * nv[:, None])[:, None],
                                ch[:, self.j_neut + 1:]], 1)
            out['channels'] = ch
            total = total + ch.sum(1)
        else:
            # all non-price costs in one head; returned so that the reported
            # channel breakdown still sums to the totals
            rest = self.rest_mean + self.rest_std * (raw_ch[:, 0] if self.masked else self.rest_head(h))
            if self.n_rs:
                base_rs = base_rs + self.rs_E_ch[0][None, :, None] * rest[:, None, :]
                rest = rest + torch.einsum('s,nsk->nk', self.rs_E_ch[0], radd)
            out['rest'] = rest
            total = total + rest
        if not self.masked:
            zg = self.g_head(h)[:, 0]
        if self.cfg.g_bound:
            g = g_max * (1. + self.cfg.g_margin) * torch.sigmoid(self.g_logit0 + zg)
        else:
            g = self.g_mean + self.g_std * zg
        out['coef'] = total
        out['g'] = g
        if self.n_rs:
            # theta(x, r) = r_theta0 + r_theta1 (r_ref / r)^beta for the outputs ``rs_names``
            # (before the water gate); both are independent of r by construction
            out['r_theta1'] = th1
            out['r_theta0'] = base_rs - th1
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
        arch = (f'masked: {len(self.sub_nets)} sub-networks of width {self.cfg.masked_width}, '
                f'residual width {self.cfg.residual_width}, unused inputs {len(self.unused_inputs)}'
                if self.masked else f'width={self.cfg.width}, depth={self.cfg.depth}')
        if self.n_rs:
            b = [round(float(v), 3) for v in self.r_beta_values().detach()]
            arch += (f', r-structured {self.rs_names} (beta={b}{" learned" if self.cfg.r_beta_learn else ""}, '
                     f'r_ref={self.r_ref:.3g} g/L/h)')
        return (f'CoefficientNet(kind={self.cfg.kind}, {arch}, '
                f'channels={self.layout["channels"]}, price_keys={self.layout["price_keys"]}, '
                f'params={n})')
