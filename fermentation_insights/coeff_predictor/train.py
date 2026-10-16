#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Training of the physics-informed coefficient network (DESIGN.md sections 4
and 6): loss, physics residuals, collocation, grouped splits, deep ensembles,
split-conformal calibration, constraint checks and checkpoints.

Loss (per mini-batch of samples; all terms are means so that the weights are
comparable)::

    L = w_mpsp * mean_s mean_i [(MPSP_hat_si - MPSP_si) / MPSP_si]^2
      + w_tot  * mean_k [(theta_hat_k - theta_k) / s_k]^2
      + w_ch   * mean_{c,k} [(theta_hat_ck - theta_ck) / s_ck]^2   (channels;
                 intensities for price-explicit channels)
      + w_g    * [(g_hat - g) / s_g]^2
      + w_phys * mean_r mean_i relu(-residual_r)^2  (labelled + collocation)

with AdamW weight decay (on weight matrices only). Physics residuals, each
normalized by a training scale and evaluated at the simulated TY points of
the (source) sample:

* ``mono_y``  ``b' + d'/t >= 0``  (dMPSP/dy <= 0)
* ``mono_t``  ``c' + d'/y >= 0``  (dMPSP/dt <= 0)
* ``positive`` ``a' + b'/y + c'/t + d'/(y t) > 0``
* ``mono_r``  ``dMPSP/dlog r <= 0`` by autograd w.r.t. the productivity input
* ``nonneg_fs``, ``nonneg_ng`` feedstock / natural-gas intensities >= 0

Collocation inputs copy training samples with log-uniform productivity over
a range ``colloc_r_factor`` wider than the data, log-uniform prices over a
widened training range and Gaussian noise on continuous features; only the
physics residuals are evaluated there.

Run ``python -m fermentation_insights.coeff_predictor.train --processed DIR
--out CKPT`` to train and save an ensemble.
"""
import copy
import json
import math
import os
import time
from dataclasses import dataclass, field, asdict

import numpy as np
import pandas as pd
import torch

from .features import FeatureSpec, ApplicabilityDomain, PRICE_COLUMNS, W_SUGAR, G_MAX
from .model import (CoefficientNet, ModelConfig, make_layout, mpsp_equation_np,
                    COEFF_KEYS, INTENSITY_CHANNEL, INTENSITY_KEYS)

__all__ = (
    'LossConfig', 'TrainConfig', 'Ensemble', 'fit_target_scales', 'make_tensors',
    'pad_points', 'compute_loss', 'physics_residuals', 'mpsp_dlogr', 'train_model', 'train_ensemble',
    'calibrate_conformal', 'check_constraints', 'apply_constraint_report',
    'split_random', 'split_group', 'split_productivity', 'make_splits', 'inner_split',
    'SCHEMES', 'load_tables', 'setup_torch',
)

RESIDUALS = ('mono_y', 'mono_t', 'positive', 'mono_r', 'nonneg_fs', 'nonneg_ng')


def setup_torch(threads=1):
    """CPU-only, few threads (simulation agents share the machine). One
    thread is the default: for these small networks it was ~10x faster than
    two threads on the shared, oversubscribed machine."""
    try:
        torch.set_num_threads(int(threads))
    except Exception:  # pragma: no cover
        pass


# %% Configurations

@dataclass
class LossConfig:
    """Loss weights and physics/collocation switches."""
    w_mpsp: float = 1.0
    w_tot: float = 1.0
    w_ch: float = 0.5
    w_g: float = 1.0
    w_phys: float = 0.1
    physics: bool = True
    residuals: dict = field(default_factory=lambda: {r: True for r in RESIDUALS})
    collocation: bool = True
    colloc_ratio: float = 1.0
    colloc_r_factor: float = 3.0
    colloc_price_range: tuple = (0.5, 1.5)
    colloc_noise: float = 0.1

    def active(self):
        """Residuals actually evaluated."""
        return [r for r in RESIDUALS if self.physics and self.residuals.get(r, False)]

    def to_dict(self):
        d = asdict(self)
        d['colloc_price_range'] = list(self.colloc_price_range)
        return d

    @classmethod
    def from_dict(cls, d):
        d = {k: v for k, v in d.items() if k in cls.__dataclass_fields__}
        if 'colloc_price_range' in d:
            d['colloc_price_range'] = tuple(d['colloc_price_range'])
        return cls(**d)


@dataclass
class TrainConfig:
    """Optimization settings."""
    epochs: int = 300
    batch_size: int = 128
    lr: float = 2e-3
    weight_decay: float = 1e-4
    patience: int = 40
    min_epochs: int = 20
    warmup_frac: float = 0.05
    lr_min_frac: float = 0.05
    val_frac: float = 0.15
    cat_dropout: float = 0.05
    grad_clip: float = 5.0
    max_points: int = 48
    seed: int = 0
    threads: int = 1
    dtype: str = 'float32'
    verbose: bool = False

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


def _dtype(name):
    return {'float32': torch.float32, 'float64': torch.float64}[name]


# %% Data handling

def load_tables(path):
    """Load ``(samples, points, configs)`` from a processed-data directory."""
    try:
        from . import dataset
        d = dataset.load(path)
        return d.samples, d.points, d.configs
    except Exception:
        pass
    out = []
    for n in ('samples', 'points', 'configs'):
        for ext in ('.csv.gz', '.csv', '.parquet'):
            p = os.path.join(path, n + ext)
            if os.path.exists(p):
                out.append(pd.read_parquet(p) if ext == '.parquet' else pd.read_csv(p, low_memory=False))
                break
        else:
            if n == 'configs':
                out.append(pd.DataFrame())
            else:
                raise FileNotFoundError(f'{n} table not found in {path}')
    return tuple(out)


def pad_points(samples, points, max_points=48):
    """Padded per-sample point arrays aligned with ``samples`` rows:
    ``Y, T, M, R, mask`` of shape (n, P) (``R`` = recovery, NaN if absent).
    Samples with more than ``max_points`` points keep an evenly spaced
    deterministic subset."""
    sid = samples['sample_id'].astype(str).to_numpy()
    pos = {s: i for i, s in enumerate(sid)}
    n = len(sid)
    if points is None or not len(points):
        z = np.ones((n, 1))
        return z, z.copy(), z.copy(), np.full((n, 1), np.nan), np.zeros((n, 1), bool)
    p = points[points['sample_id'].astype(str).isin(pos)]
    y, t, M = (pd.to_numeric(p[c], errors='coerce').to_numpy(float) for c in ('y', 't', 'MPSP'))
    rec = pd.to_numeric(p['recovery'], errors='coerce').to_numpy(float) if 'recovery' in p else np.full(len(p), np.nan)
    ok = np.isfinite(y) & np.isfinite(t) & np.isfinite(M) & (y > 0) & (t > 0)
    rows = p['sample_id'].astype(str).map(pos).to_numpy()[ok]
    y, t, M, rec = y[ok], t[ok], M[ok], rec[ok]
    order = np.argsort(rows, kind='stable')
    rows, y, t, M, rec = rows[order], y[order], t[order], M[order], rec[order]
    counts = np.bincount(rows, minlength=n)
    P = int(max(1, min(max_points, counts.max() if n else 1)))
    Y, T, MM = np.ones((n, P)), np.ones((n, P)), np.ones((n, P))
    RR = np.full((n, P), np.nan)
    mask = np.zeros((n, P), bool)
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
    for i in np.flatnonzero(counts):
        c = counts[i]
        sel = starts[i] + (np.arange(c) if c <= P else np.unique(np.linspace(0, c - 1, P).round().astype(int)))
        k = len(sel)
        Y[i, :k], T[i, :k], MM[i, :k], RR[i, :k] = y[sel], t[sel], M[sel], rec[sel]
        mask[i, :k] = True
    return Y, T, MM, RR, mask


def _targets(samples, layout):
    n = len(samples)

    def col(c):
        return pd.to_numeric(samples[c], errors='coerce').to_numpy(float) if c in samples else np.full(n, np.nan)

    coef = np.column_stack([col(k) for k in COEFF_KEYS])
    g = col('g')
    ch = np.stack([np.column_stack([col(f'ch__{c}__{k}') for k in COEFF_KEYS]) for c in layout['channels']], 1) \
        if layout['channels'] else np.zeros((n, 0, 4))
    q = np.stack([np.column_stack([col(f'int__{j}__{k}') for k in COEFF_KEYS]) for j in layout['price_keys']], 1) \
        if layout['price_keys'] else np.zeros((n, 0, 4))
    return coef, g, ch, q


def _robust_std(x, floor):
    x = x[np.isfinite(x)]
    s = float(x.std()) if len(x) > 1 else 0.
    return max(s, floor)


def fit_target_scales(samples, points, layout, batch):
    """Per-target means/standard deviations and physics-residual scales
    computed on the training rows only. ``batch`` is the FeatureBatch of
    ``samples`` (for prices, w_sugar and g_max)."""
    coef, g, ch, q = _targets(samples, layout)
    sc = {}
    cm = np.nanmean(coef, 0)
    cs = np.array([_robust_std(coef[:, k], 1e-6 * max(1., abs(cm[k]))) for k in range(4)])
    sc['coef_mean'], sc['coef_std'] = cm.tolist(), cs.tolist()
    sc['g_mean'] = float(np.nanmean(g))
    # floor: recovery errors below ~1 % are immaterial (MPSP ~ 1/g) and a nearly
    # constant g must not dominate the loss / early-stopping criterion
    sc['g_std'] = _robust_std(g, 0.01 * max(abs(sc['g_mean']), 0.1))
    if ch.shape[1]:
        chm = np.nan_to_num(np.nanmean(ch, 0))
        chs = np.array([[_robust_std(ch[:, c, k], 0.01 * cs[k]) for k in range(4)] for c in range(ch.shape[1])])
        sc['ch_mean'], sc['ch_std'] = chm.tolist(), chs.tolist()
    else:
        sc['ch_mean'], sc['ch_std'] = [], []
    if q.shape[1]:
        qm = np.nan_to_num(np.nanmean(q, 0))
        qs = []
        for j in range(q.shape[1]):
            s = np.array([_robust_std(q[:, j, k], 0.) for k in range(4)])
            floor = max(0.01 * s.max(), 1e-8)
            qs.append(np.maximum(s, floor))
        sc['int_mean'], sc['int_std'] = qm.tolist(), np.array(qs).tolist()
        pidx = [INTENSITY_KEYS.index(k) for k in layout['price_keys']]
        rest = coef - np.einsum('nj,njk->nk', batch.prices[:, pidx], np.nan_to_num(q))
    else:
        sc['int_mean'], sc['int_std'] = [], []
        rest = coef
    rm = np.nanmean(rest, 0)
    sc['rest_mean'] = rm.tolist()
    sc['rest_std'] = [_robust_std(rest[:, k], 1e-6 * max(1., abs(rm[k]))) for k in range(4)]
    # structured feedstock intensity: z = log(b'_fs * w_sugar) = log(1/eta)
    if 'fs' in layout['price_keys']:
        z = np.log(np.clip(q[:, layout['price_keys'].index('fs'), 1] * batch.w_sugar, 1e-6, None))
        z = z[np.isfinite(z)]
        sc['fs_logb_center'] = float(np.median(z)) if len(z) else 0.1
        sc['fs_logb_scale'] = float(max(z.std() if len(z) > 1 else 0., 0.02))
    ratio = g / batch.g_max
    ratio = ratio[np.isfinite(ratio)]
    rmed = float(np.clip(np.median(ratio) if len(ratio) else 0.9, 0.05, 0.995))
    sc['g_logit0'] = float(math.log(rmed / (1 - rmed)))
    Y, T, M, _, mask = pad_points(samples, points)
    Mv = M[mask]
    sc['mpsp_scale'] = float(np.median(np.abs(Mv))) if len(Mv) else 1.
    sc['num_scale'] = float(np.median(np.abs(M * g[:, None])[mask])) if len(Mv) else 1.
    ips = []
    for j in range(q.shape[1]):
        I = mpsp_equation_np(np.nan_to_num(q[:, j])[:, None, :], 1., Y, T)
        v = np.abs(I[mask])
        ips.append(float(max(v.mean() if len(v) else 0., 1e-6)))
    sc['int_point_scale'] = ips
    return sc


def make_tensors(spec, samples, points, layout, batch=None, max_points=48, dtype=torch.float32):
    """Tensors for training/evaluation (inputs, targets, masks, points)."""
    if batch is None:
        batch = spec.transform(samples)
    coef, g, ch, q = _targets(samples, layout)
    Y, T, M, R, mask = pad_points(samples, points, max_points)

    def tt(a):
        return torch.tensor(np.asarray(a, float), dtype=dtype)

    d = dict(X=tt(batch.X), log_r=tt(batch.log_r), prices=tt(batch.prices),
             w_sugar=tt(batch.w_sugar), g_max=tt(batch.g_max),
             coef=tt(np.nan_to_num(coef)), coef_m=tt(np.isfinite(coef)),
             g=tt(np.nan_to_num(g, nan=1.)), g_m=tt(np.isfinite(g)),
             ch=tt(np.nan_to_num(ch)), ch_m=tt(np.isfinite(ch)),
             q=tt(np.nan_to_num(q)), q_m=tt(np.isfinite(q)),
             Y=tt(Y), T=tt(T), M=tt(M), mask=tt(mask))
    return d


def _index(tb, idx):
    return {k: v[idx] for k, v in tb.items()}


# %% Loss

def _mmean(x, m):
    return (torch.where(m > 0, x, torch.zeros_like(x)) * m).sum() / m.sum().clamp(min=1.)


def data_loss(model, out, b, lcfg):
    """Data terms of the loss (dict of scalar tensors)."""
    sc = model.scales_np
    dt = out['coef'].dtype
    L = {}
    mp = model.mpsp(out, b['Y'], b['T'])
    den = b['M'].abs().clamp(min=0.05 * sc['mpsp_scale'])
    rel2 = ((mp - b['M']) / den) ** 2
    m = b['mask']
    per = (rel2 * m).sum(1) / m.sum(1).clamp(min=1.)
    has = (m.sum(1) > 0).to(dt)
    L['mpsp'] = (per * has).sum() / has.sum().clamp(min=1.)
    cs = torch.as_tensor(sc['coef_std'], dtype=dt)
    L['tot'] = _mmean(((out['coef'] - b['coef']) / cs) ** 2, b['coef_m'])
    L['g'] = _mmean(((out['g'] - b['g']) / sc['g_std']) ** 2, b['g_m'])
    parts, masks = [], []
    if 'channels' in out and out['channels'].shape[1]:
        s = torch.as_tensor(sc['ch_std'], dtype=dt)
        parts.append((((out['channels'] - b['ch']) / s) ** 2).flatten(1))
        masks.append(b['ch_m'].flatten(1))
    if 'intensities' in out and out['intensities'].shape[1]:
        s = torch.as_tensor(sc['int_std'], dtype=dt)
        parts.append((((out['intensities'] - b['q']) / s) ** 2).flatten(1))
        masks.append(b['q_m'].flatten(1))
    L['ch'] = _mmean(torch.cat(parts, 1), torch.cat(masks, 1)) if parts else torch.zeros((), dtype=dt)
    L['data'] = lcfg.w_mpsp * L['mpsp'] + lcfg.w_tot * L['tot'] + lcfg.w_g * L['g'] + lcfg.w_ch * L['ch']
    return L


def mpsp_dlogr(out, log_r, Y, T):
    """dMPSP/dlog r at points (n, P) by autograd (differentiable,
    ``create_graph=True``). ``log_r`` is the leaf tensor of the forward pass;
    each sample's outputs depend only on its own ``log_r``, so gradients of
    batch sums give per-sample derivatives of a', b', c', d', g, combined by
    the chain rule through the equation layer."""
    coef, g = out['coef'], out['g']
    outs = [coef[:, k] for k in range(4)] + [g]
    grads = [torch.autograd.grad(o.sum(), log_r, create_graph=True, allow_unused=True)[0] for o in outs]
    grads = [torch.zeros_like(log_r) if gr is None else gr for gr in grads]
    da, db, dc, dd, dg = (gr[:, None] for gr in grads)
    a, b, c, d = (coef[:, k:k + 1] for k in range(4))
    N = a + b / Y + c / T + d / (Y * T)
    return (da + db / Y + dc / T + dd / (Y * T)) / g[:, None] - N * dg / g[:, None] ** 2


def physics_residuals(model, out, Y, T, mask, active, log_r=None):
    """Mean squared normalized constraint violations at points (n, P).
    ``log_r`` must be the leaf tensor (requires_grad) used in the forward pass
    for ``mono_r``. Returns {name: scalar tensor}."""
    sc = model.scales_np
    coef, g = out['coef'], out['g']
    a, b, c, d = (coef[:, k:k + 1] for k in range(4))
    res = {}
    if 'mono_y' in active:
        res['mono_y'] = _mmean(torch.relu(-(b + d / T)) ** 2, mask) / sc['coef_std'][1] ** 2
    if 'mono_t' in active:
        res['mono_t'] = _mmean(torch.relu(-(c + d / Y)) ** 2, mask) / sc['coef_std'][2] ** 2
    if 'positive' in active:
        N = a + b / Y + c / T + d / (Y * T)
        res['positive'] = _mmean(torch.relu(-N) ** 2, mask) / sc['num_scale'] ** 2
    for res_name, key, name in (('nonneg_fs', 'fs', 'feedstock'), ('nonneg_ng', 'ng', 'natural_gas')):
        if res_name not in active:
            continue
        pk, chn = model.layout['price_keys'], model.layout['channels']
        if key in pk and 'intensities' in out:
            j = pk.index(key)
            q, scale = out['intensities'][:, j], sc['int_point_scale'][j]
        elif name in chn and 'channels' in out:
            q, scale = out['channels'][:, chn.index(name)], sc['num_scale']
        else:
            continue
        I = q[:, 0:1] + q[:, 1:2] / Y + q[:, 2:3] / T + q[:, 3:4] / (Y * T)
        res[res_name] = _mmean(torch.relu(-I) ** 2, mask) / scale ** 2
    if 'mono_r' in active and log_r is not None and log_r.requires_grad:
        dM = mpsp_dlogr(out, log_r, Y, T)
        res['mono_r'] = _mmean(torch.relu(dM) ** 2, mask) / sc['mpsp_scale'] ** 2
    return res


def _collocation(tb, idx, model, spec_info, lcfg, gen):
    """Perturbed copies of the training rows ``idx`` (inputs + points only)."""
    n = max(1, int(round(lcfg.colloc_ratio * len(idx))))
    src = idx[torch.randint(len(idx), (n,), generator=gen)]
    b = {k: tb[k][src] for k in ('X', 'log_r', 'prices', 'w_sugar', 'g_max', 'Y', 'T', 'mask')}
    dt = b['X'].dtype
    lo = spec_info['log_r_min'] - math.log(lcfg.colloc_r_factor)
    hi = spec_info['log_r_max'] + math.log(lcfg.colloc_r_factor)
    b['log_r'] = lo + (hi - lo) * torch.rand(n, generator=gen, dtype=dt)
    P = b['prices'].clone()
    f0, f1 = lcfg.colloc_price_range
    for j, (pmin, pmax) in enumerate(spec_info['price_range']):
        if pmin > 0 and pmax > 0:
            l0, l1 = math.log(pmin * f0), math.log(pmax * f1)
            P[:, j] = torch.exp(l0 + (l1 - l0) * torch.rand(n, generator=gen, dtype=dt))
        else:
            span = pmax - pmin
            P[:, j] = (pmin - (1 - f0) * span) + (span * (1 + (f1 - f0))) * torch.rand(n, generator=gen, dtype=dt)
    b['prices'] = P
    if lcfg.colloc_noise > 0:
        cm = spec_info['cont_mask']
        noise = torch.randn(b['X'].shape, generator=gen, dtype=dt) * lcfg.colloc_noise
        b['X'] = b['X'] + noise * cm
    return b


def compute_loss(model, tb, idx, lcfg, spec_info=None, gen=None, train=True, cat_dropout=0.):
    """Total loss for the rows ``idx`` of tensors ``tb``; returns
    ``(loss, parts)``. ``cat_dropout`` replaces categorical blocks by their
    ``__other__`` level with this probability (training only)."""
    b = _index(tb, idx)
    if train and cat_dropout > 0 and spec_info is not None:
        b['X'] = _cat_dropout(b['X'], spec_info['blocks'], cat_dropout, gen)
    active = lcfg.active()
    need_r = 'mono_r' in active
    use_col = train and lcfg.collocation and bool(active) and spec_info is not None
    X, lr, P, w, gm = b['X'], b['log_r'], b['prices'], b['w_sugar'], b['g_max']
    if use_col:
        cb = _collocation(tb, idx, model, spec_info, lcfg, gen)
        X = torch.cat([X, cb['X']])
        lr = torch.cat([lr, cb['log_r']])
        P, w, gm = torch.cat([P, cb['prices']]), torch.cat([w, cb['w_sugar']]), torch.cat([gm, cb['g_max']])
    if need_r and train:
        lr = lr.detach().clone().requires_grad_(True)
    out = model(X, lr, P, w, gm)
    n = len(idx)
    out_l = {k: v[:n] for k, v in out.items()}
    parts = data_loss(model, out_l, b, lcfg)
    loss = parts['data']
    if active and train:
        if use_col:
            Y = torch.cat([b['Y'], cb['Y']])
            T = torch.cat([b['T'], cb['T']])
            M = torch.cat([b['mask'], cb['mask']])
        else:
            Y, T, M = b['Y'], b['T'], b['mask']
        res = physics_residuals(model, out, Y, T, M, active, log_r=lr)
        if res:
            phys = sum(res.values()) / len(res)
            parts.update({f'phys_{k}': v for k, v in res.items()})
            parts['phys'] = phys
            loss = loss + lcfg.w_phys * phys
    parts['loss'] = loss
    return loss, parts


# %% Training

def _spec_info(spec, tb, dtype):
    P = tb['prices'].double().numpy()
    lr = tb['log_r'].double().numpy()
    pr = [(float(P[:, j].min()), float(P[:, j].max())) for j in range(P.shape[1])]
    return dict(log_r_min=float(np.nanmin(lr)), log_r_max=float(np.nanmax(lr)),
                price_range=pr, cont_mask=torch.as_tensor(spec.continuous_mask(), dtype=dtype),
                blocks=spec.categorical_blocks())


def _cat_dropout(X, blocks, p, gen):
    if p <= 0 or not blocks:
        return X
    X = X.clone()
    for (s, e, o) in blocks:
        drop = torch.rand(X.shape[0], generator=gen) < p
        if drop.any():
            X[drop, s:e] = 0.
            X[drop, s + o] = 1.
    return X


def build_model(spec, layout, scales, mcfg, dtype=torch.float32):
    m = CoefficientNet(spec.n_features, mcfg, layout, scales, spec.r_stats).to(dtype)
    m.scales_np = scales
    return m


def train_model(spec, samples, points, mcfg=None, lcfg=None, tcfg=None, val=None, layout=None,
                scales=None, batch=None, seed=None, log=None, tensors=None, val_tensors=None):
    """Train one network.

    Parameters
    ----------
    spec : FeatureSpec fitted on (a superset of) the training samples.
    samples, points : training tables.
    val : (samples_val, points_val) for early stopping (optional).
    layout, scales : precomputed (default: from the training data).
    tensors, val_tensors : precomputed :func:`make_tensors` dicts (then
        ``samples``/``points``/``val`` are only used for defaults).

    Returns ``(model, history)`` with the best-validation state restored.
    """
    mcfg = mcfg or ModelConfig()
    lcfg = lcfg or LossConfig()
    tcfg = tcfg or TrainConfig()
    setup_torch(tcfg.threads)
    seed = tcfg.seed if seed is None else seed
    torch.manual_seed(seed)
    np.random.seed(seed % (2 ** 32))
    gen = torch.Generator().manual_seed(seed)
    dt = _dtype(tcfg.dtype)
    layout = layout or make_layout(mcfg, samples.columns)
    if scales is None or tensors is None:
        batch = spec.transform(samples) if batch is None else batch
    scales = scales or fit_target_scales(samples, points, layout, batch)
    tb = tensors if tensors is not None else make_tensors(spec, samples, points, layout, batch,
                                                          tcfg.max_points, dt)
    vb = val_tensors
    if vb is None and val is not None and len(val[0]):
        vb = make_tensors(spec, val[0], val[1], layout, None, tcfg.max_points, dt)
    if vb is not None and len(vb['X']) == 0:
        vb = None
    model = build_model(spec, layout, scales, mcfg, dt)
    info = _spec_info(spec, tb, dt)
    decay, no_decay = [], []
    for name, p in model.named_parameters():
        (decay if p.ndim >= 2 else no_decay).append(p)
    opt = torch.optim.AdamW([{'params': decay, 'weight_decay': tcfg.weight_decay},
                             {'params': no_decay, 'weight_decay': 0.}], lr=tcfg.lr)
    n = len(tb['X'])
    steps_per_epoch = max(1, math.ceil(n / tcfg.batch_size))
    total = steps_per_epoch * tcfg.epochs
    warm = max(1, int(tcfg.warmup_frac * total))

    def lr_lambda(step):
        if step < warm:
            return (step + 1) / warm
        f = min(1., (step - warm) / max(1, total - warm))
        return tcfg.lr_min_frac + (1 - tcfg.lr_min_frac) * 0.5 * (1 + math.cos(math.pi * f))

    sched = torch.optim.lr_scheduler.LambdaLR(opt, lr_lambda)
    best, best_state, bad, hist = float('inf'), None, 0, []
    t0 = time.time()
    for epoch in range(tcfg.epochs):
        model.train()
        perm = torch.randperm(n, generator=gen)
        agg = {}
        for s in range(0, n, tcfg.batch_size):
            idx = perm[s:s + tcfg.batch_size]
            loss, parts = compute_loss(model, tb, idx, lcfg, info, gen, train=True,
                                       cat_dropout=tcfg.cat_dropout)
            opt.zero_grad(set_to_none=True)
            loss.backward()
            if tcfg.grad_clip:
                torch.nn.utils.clip_grad_norm_(model.parameters(), tcfg.grad_clip)
            opt.step()
            sched.step()
            for k, v in parts.items():
                agg[k] = agg.get(k, 0.) + float(v.detach()) * len(idx) / n
        rec = {'epoch': epoch, **{f'train_{k}': v for k, v in agg.items()}}
        crit = agg['data']
        if vb is not None:
            model.eval()
            with torch.no_grad():
                _, vp = compute_loss(model, vb, torch.arange(len(vb['X'])), lcfg, train=False)
            rec.update({f'val_{k}': float(v) for k, v in vp.items()})
            crit = float(vp['data'])
        hist.append(rec)
        if not np.isfinite(crit):
            break
        if not np.isfinite(best) or crit < best - 1e-6 * abs(best):
            best, bad = crit, 0
            best_state = copy.deepcopy(model.state_dict())
        else:
            bad += 1
        if tcfg.verbose and (epoch % 20 == 0 or epoch == tcfg.epochs - 1):
            msg = f'[seed {seed}] epoch {epoch:4d} train {agg["loss"]:.4g} crit {crit:.4g} ({time.time() - t0:.0f}s)'
            (log or print)(msg)
        if epoch >= tcfg.min_epochs and bad >= tcfg.patience:
            break
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    return model, hist


# %% Splits

def split_random(samples, n_splits=5, seed=0):
    """Random sample-level K-fold (TY points always follow their sample)."""
    n = len(samples)
    perm = np.random.RandomState(seed).permutation(n)
    folds = np.array_split(perm, n_splits)
    return [(np.sort(np.concatenate([f for j, f in enumerate(folds) if j != i])), np.sort(folds[i]), f'fold{i}')
            for i in range(n_splits)]


def split_group(samples, col, n_splits=None, seed=0):
    """Grouped splits on ``col``: leave-one-group-out when ``n_splits`` is
    None, else K folds of shuffled groups."""
    g = samples[col].astype(str).to_numpy()
    groups = np.array(sorted(set(g)))
    if n_splits is None or n_splits >= len(groups):
        parts = [[x] for x in groups]
    else:
        perm = np.random.RandomState(seed).permutation(len(groups))
        parts = [list(groups[p]) for p in np.array_split(perm, n_splits)]
    out = []
    for p in parts:
        te = np.isin(g, p)
        out.append((np.flatnonzero(~te), np.flatnonzero(te), '+'.join(p)))
    return out


def _rel_productivity(samples):
    r = pd.to_numeric(samples['r'], errors='coerce').to_numpy(float)
    if 'baseline_productivity' in samples:
        r0 = pd.to_numeric(samples['baseline_productivity'], errors='coerce').to_numpy(float)
        rel = np.where(np.isfinite(r0) & (r0 > 0), r / np.where(r0 > 0, r0, 1.), np.nan)
        if np.isfinite(rel).mean() > 0.9:
            return np.log(rel)
    return np.log(r)


def split_productivity(samples, central=(0.2, 0.8)):
    """Productivity extrapolation: train on the central quantile range of
    log(r / baseline productivity) (log r if unavailable), test outside."""
    lr = _rel_productivity(samples)
    lo, hi = np.nanquantile(lr, central)
    tr = (lr >= lo) & (lr <= hi)
    return [(np.flatnonzero(tr), np.flatnonzero(~tr), f'central{central[0]:.2f}-{central[1]:.2f}')]


SCHEMES = {
    'random': lambda s, k, seed: split_random(s, k, seed),
    'config': lambda s, k, seed: split_group(s, 'config_id', None, seed),
    'feedstock': lambda s, k, seed: split_group(s, 'feedstock', None, seed),
    'product': lambda s, k, seed: split_group(s, 'fermentation_product', None, seed),
    'productivity': lambda s, k, seed: split_productivity(s),
}
#: grouping column used for inner validation splits of each scheme
SCHEME_GROUP = {'random': None, 'config': 'config_id', 'feedstock': 'config_id',
                'product': 'config_id', 'productivity': None}


def make_splits(samples, scheme, n_splits=5, seed=0):
    """List of ``(train_idx, test_idx, label)`` (positional indices)."""
    if scheme not in SCHEMES:
        raise ValueError(f'unknown scheme {scheme!r}; choose from {sorted(SCHEMES)}')
    return SCHEMES[scheme](samples, n_splits, seed)


def inner_split(samples, idx, frac=0.15, seed=0, group_col=None):
    """Split positional indices ``idx`` into (fit, val), grouped by
    ``group_col`` when given and when there are >= 4 groups."""
    idx = np.asarray(idx)
    rng = np.random.RandomState(seed)
    if frac <= 0 or len(idx) < 4:
        return idx, idx[:0]
    if group_col is not None and group_col in samples:
        g = samples[group_col].astype(str).to_numpy()[idx]
        u = np.array(sorted(set(g)))
        if len(u) >= 4:
            k = max(1, int(round(frac * len(u))))
            vg = rng.choice(u, k, replace=False)
            m = np.isin(g, vg)
            return idx[~m], idx[m]
    perm = rng.permutation(len(idx))
    k = max(1, int(round(frac * len(idx))))
    return np.sort(idx[perm[k:]]), np.sort(idx[perm[:k]])


# %% Constraint checks on data

def _np_targets(samples, points, max_points=48):
    coef = np.column_stack([pd.to_numeric(samples[k], errors='coerce').to_numpy(float) for k in COEFF_KEYS])
    g = pd.to_numeric(samples['g'], errors='coerce').to_numpy(float)
    Y, T, M, _, mask = pad_points(samples, points, max_points)
    return coef, g, Y, T, mask


def _all_points(cond, mask):
    """Per-sample: condition holds at all valid points (NaN if no points)."""
    ok = np.where(mask, cond, True).all(1)
    has = mask.any(1)
    return np.where(has, ok.astype(float), np.nan)


def constraint_case_values(coef, g, Y, T, mask, scales=None, q_fs=None, q_ng=None, tol=1e-3):
    """Per-sample satisfaction (1/0/NaN) of the residual constraints for
    coefficient arrays (targets or predictions)."""
    s_b = scales['coef_std'][1] if scales else np.nanstd(coef[:, 1]) or 1.
    s_c = scales['coef_std'][2] if scales else np.nanstd(coef[:, 2]) or 1.
    a, b, c, d = (coef[:, k:k + 1] for k in range(4))
    N = a + b / Y + c / T + d / (Y * T)
    out = {'mono_y': _all_points(b + d / T >= -tol * s_b, mask),
           'mono_t': _all_points(c + d / Y >= -tol * s_c, mask),
           'positive': _all_points(N > 0, mask)}
    for name, q in (('nonneg_fs', q_fs), ('nonneg_ng', q_ng)):
        if q is not None:
            I = mpsp_equation_np(q[:, None, :], 1., Y, T)
            sc = np.nanmean(np.abs(I[mask])) if mask.any() else 1.
            v = _all_points(I >= -tol * max(sc, 1e-12), mask)
            v[~np.isfinite(q).all(1)] = np.nan
            out[name] = v
    return out


def _mono_r_pairs(samples, coef, g, Y, T, mask, tol=1e-3):
    """Fraction of consecutive productivity pairs (same configuration and
    identical parameters otherwise) with non-increasing fitted MPSP."""
    keycols = [c for c in samples.columns if (c.startswith(('p_', 'bp_')) or c in PRICE_COLUMNS)]
    key = samples['config_id'].astype(str)
    if keycols:
        V = np.column_stack([pd.to_numeric(samples[c], errors='coerce').to_numpy(float)
                             for c in dict.fromkeys(keycols)])
        with np.errstate(invalid='ignore'):
            scale = np.nanmax(np.abs(V), 0) if len(V) else np.ones(V.shape[1])
        scale = np.where(np.isfinite(scale) & (scale > 0), scale, 1.)
        Vr = np.round(V / scale, 9)
        Vr = np.where(np.isfinite(Vr), Vr, -9e9) + 0.  # +0. maps -0.0 to 0.0
        key = key + '|' + pd.Series([row.tobytes().hex() for row in Vr], index=samples.index)
    r = pd.to_numeric(samples['r'], errors='coerce').to_numpy(float)
    res = []
    for _, ix in pd.Series(np.arange(len(samples))).groupby(key.to_numpy()):
        ix = ix.to_numpy()
        if len(ix) < 2:
            continue
        ix = ix[np.argsort(r[ix])]
        for i, j in zip(ix[:-1], ix[1:]):
            if not (r[j] > r[i] * (1 + 1e-9)) or not mask[i].any():
                continue
            m = mask[i]
            M1 = mpsp_equation_np(coef[i], g[i], Y[i, m], T[i, m])
            M2 = mpsp_equation_np(coef[j], g[j], Y[i, m], T[i, m])
            if np.isfinite(M1).all() and np.isfinite(M2).all():
                res.append(bool(np.all(M2 <= M1 + tol * np.abs(M1))))
    return np.array(res, float)


def check_constraints(samples, points, threshold=0.99, tol=1e-3):
    """Check whether each physics constraint holds in the training targets.

    Returns a DataFrame with columns ``constraint, type, n_cases,
    frac_satisfied, enabled, note``. A residual is enabled only if it holds in
    >= ``threshold`` of the checkable cases (``mono_r`` is enabled on the
    physical prior when no productivity series exist). Structural constraints
    (``g_le_gmax``, ``fs_eta_le_1``) decide whether the hard bounds of the
    model may be used; identities (``channel_sum``, ``price_linear``) only
    validate the data processing."""
    coef, g, Y, T, mask = _np_targets(samples, points)
    n = len(samples)

    def col(c):
        return pd.to_numeric(samples[c], errors='coerce').to_numpy(float) if c in samples else np.full(n, np.nan)

    def qarr(key):
        cols = [f'int__{key}__{k}' for k in COEFF_KEYS]
        return np.column_stack([col(c) for c in cols]) if all(c in samples for c in cols) else None

    q_fs, q_ng = qarr('fs'), qarr('ng')
    vals = constraint_case_values(coef, g, Y, T, mask, None, q_fs, q_ng, tol)
    rows = []

    def add(name, kind, v, note='', default=False, stat=float('nan')):
        v = np.asarray(v, float)
        v = v[np.isfinite(v)]
        frac = float(v.mean()) if len(v) else float('nan')
        en = bool(frac >= threshold) if len(v) else default
        rows.append(dict(constraint=name, type=kind, n_cases=int(len(v)), frac_satisfied=frac,
                         enabled=en, stat=stat, note=note if len(v) else (note or 'not checkable')))

    s_tot = np.nanstd(coef, 0) + 1e-12  # identities: errors relative to the total's spread

    for k in ('mono_y', 'mono_t', 'positive'):
        add(k, 'residual', vals[k])
    pairs = _mono_r_pairs(samples, coef, g, Y, T, mask, tol)
    add('mono_r', 'residual', pairs, 'consecutive productivity pairs', default=True)
    for k, lab in (('nonneg_fs', 'feedstock intensity'), ('nonneg_ng', 'natural-gas intensity')):
        add(k, 'residual', vals.get(k, []), lab if k in vals else 'no intensity targets')
    gm = col(G_MAX)
    ok = np.isfinite(g) & np.isfinite(gm)
    gr = float(np.max(g[ok] / gm[ok])) if ok.any() else float('nan')
    add('g_le_gmax', 'structural', np.where(ok, (g <= gm * (1 + tol)).astype(float), np.nan),
        f'max g/g_max = {gr:.4f}', stat=gr)
    if q_fs is not None:
        w = col(W_SUGAR)
        v = q_fs[:, 1] * w
        add('fs_eta_le_1', 'structural', np.where(np.isfinite(v), (v >= 1 - tol).astype(float), np.nan),
            "b'_fs * w_sugar >= 1")
    chs = sorted({c.split('__')[1] for c in samples.columns if c.startswith('ch__') and c.count('__') == 2})
    if chs:
        S = np.column_stack([np.nansum(np.column_stack([col(f'ch__{c}__{k}') for c in chs]), 1) for k in COEFF_KEYS])
        err = np.abs(S - coef) / s_tot
        add('channel_sum', 'identity', np.where(np.isfinite(err).all(1), (err.max(1) <= 1e-4).astype(float), np.nan))
    lin = []
    for key, ch in INTENSITY_CHANNEL.items():
        q = qarr(key)
        pc = PRICE_COLUMNS[INTENSITY_KEYS.index(key)]
        if q is None or pc not in samples or f'ch__{ch}__a' not in samples:
            continue
        C = np.column_stack([col(f'ch__{ch}__{k}') for k in COEFF_KEYS])
        e = np.abs(C - col(pc)[:, None] * q) / s_tot
        lin.append(np.where(np.isfinite(e).all(1), (e.max(1) <= 1e-4).astype(float), np.nan))
    if lin:
        add('price_linear', 'identity', np.nanmin(np.column_stack(lin), 1))
    return pd.DataFrame(rows)


def apply_constraint_report(report, mcfg, lcfg, max_g_margin=0.02):
    """Return copies of (mcfg, lcfg) with residuals / hard bounds disabled
    where the training data do not satisfy them. If ``g > g_max`` occurs but
    by at most ``max_g_margin`` (relative), the bound is relaxed to
    ``g_max * (1 + g_margin)`` instead of being dropped."""
    mcfg, lcfg = copy.deepcopy(mcfg), copy.deepcopy(lcfg)
    en = dict(zip(report['constraint'], report['enabled']))
    res = dict(lcfg.residuals)
    for r in RESIDUALS:
        if r in en and not en[r]:
            res[r] = False
    lcfg.residuals = res
    if 'g_le_gmax' in en and not en['g_le_gmax']:
        # small excesses (rounding, impurities) -> widen the hard bound instead
        gr = report.set_index('constraint').loc['g_le_gmax', 'stat'] if 'stat' in report else np.nan
        if np.isfinite(gr) and gr <= 1 + max_g_margin:
            mcfg.g_margin = max(mcfg.g_margin, float(gr - 1) + 1e-3)
        else:
            mcfg.g_bound = False
    if 'fs_eta_le_1' in en and not en['fs_eta_le_1']:
        mcfg.fs_eta_bound = False
    return mcfg, lcfg


# %% Ensemble

class Ensemble:
    """Deep ensemble with shared FeatureSpec, layout and target scales."""

    def __init__(self, members, spec, layout, scales, mcfg, lcfg, tcfg, conformal=None,
                 constraint_report=None, ad=None, meta=None):
        self.members, self.spec, self.layout, self.scales = members, spec, layout, scales
        self.mcfg, self.lcfg, self.tcfg = mcfg, lcfg, tcfg
        self.conformal = conformal or {}
        self.constraint_report = constraint_report
        self.ad = ad
        self.meta = meta or {}

    @property
    def dtype(self):
        return next(self.members[0].parameters()).dtype

    def channel_names(self):
        return self.members[0].channel_names()

    def _inputs(self, batch):
        dt = self.dtype
        return [torch.tensor(np.asarray(a, float), dtype=dt)
                for a in (batch.X, batch.log_r, batch.prices, batch.w_sugar, batch.g_max)]

    def member_outputs(self, batch):
        """Stacked numpy outputs over members: ``coef`` (M, n, 4), ``g`` (M, n),
        ``channels`` (M, n, C_all, 4) incl. price channels first, and
        ``intensities`` (M, n, J, 4)."""
        inp = self._inputs(batch)
        res = {}
        with torch.no_grad():
            for m in self.members:
                m.eval()
                o = m(*inp)
                chans = []
                if 'price_channels' in o:
                    chans.append(o['price_channels'])
                if 'channels' in o:
                    chans.append(o['channels'])
                d = {'coef': o['coef'], 'g': o['g']}
                if chans:
                    d['channels'] = torch.cat(chans, 1)
                if 'intensities' in o:
                    d['intensities'] = o['intensities']
                for k, v in d.items():
                    res.setdefault(k, []).append(v.double().numpy())
        return {k: np.stack(v) for k, v in res.items()}

    def mpsp_members(self, outs, Y, T):
        """(M, n, P) MPSP of every member at points (n, P)."""
        return mpsp_equation_np(outs['coef'][:, :, None, :], outs['g'][:, :, None], Y[None], T[None])

    def predict_samples(self, samples, strict=True, descriptors=None):
        """Mean/std of a', b', c', d', g for a samples DataFrame."""
        batch = self.spec.transform(samples, descriptors=descriptors, strict=strict)
        return batch, self.member_outputs(batch)

    # ---- persistence ----
    def save(self, path):
        """Save to a directory: ``meta.json``, ``feature_spec.json``,
        ``target_scales.json``, ``member_<i>.pt`` (state_dicts) and
        ``ad_train.npy``."""
        os.makedirs(path, exist_ok=True)
        self.spec.to_json(os.path.join(path, 'feature_spec.json'))
        with open(os.path.join(path, 'target_scales.json'), 'w') as fh:
            json.dump(self.scales, fh, indent=1)
        for i, m in enumerate(self.members):
            torch.save(m.state_dict(), os.path.join(path, f'member_{i}.pt'))
        meta = dict(format=1, n_members=len(self.members), model_config=self.mcfg.to_dict(),
                    loss_config=self.lcfg.to_dict(), train_config=self.tcfg.to_dict(),
                    layout=self.layout, conformal=self.conformal,
                    constraint_report=(self.constraint_report.to_dict('records')
                                       if isinstance(self.constraint_report, pd.DataFrame) else self.constraint_report),
                    ad=self.ad.to_dict() if self.ad is not None else None,
                    dtype=str(self.dtype).replace('torch.', ''), meta=self.meta)
        with open(os.path.join(path, 'meta.json'), 'w') as fh:
            json.dump(_jsonable(meta), fh, indent=1)
        if self.ad is not None:
            np.save(os.path.join(path, 'ad_train.npy'), self.ad.Z.astype(np.float32))
        return path

    @classmethod
    def load(cls, path):
        with open(os.path.join(path, 'meta.json')) as fh:
            meta = json.load(fh)
        spec = FeatureSpec.from_json(os.path.join(path, 'feature_spec.json'))
        with open(os.path.join(path, 'target_scales.json')) as fh:
            scales = json.load(fh)
        mcfg = ModelConfig.from_dict(meta['model_config'])
        lcfg = LossConfig.from_dict(meta['loss_config'])
        tcfg = TrainConfig.from_dict(meta['train_config'])
        dt = _dtype(meta.get('dtype', 'float32'))
        members = []
        for i in range(meta['n_members']):
            m = build_model(spec, meta['layout'], scales, mcfg, dt)
            m.load_state_dict(torch.load(os.path.join(path, f'member_{i}.pt'), weights_only=True))
            m.eval()
            members.append(m)
        ad = None
        if meta.get('ad') is not None and os.path.exists(os.path.join(path, 'ad_train.npy')):
            ad = ApplicabilityDomain.from_dict(meta['ad'], np.load(os.path.join(path, 'ad_train.npy')))
        rep = meta.get('constraint_report')
        rep = pd.DataFrame(rep) if rep is not None else None
        return cls(members, spec, meta['layout'], scales, mcfg, lcfg, tcfg, meta.get('conformal'),
                   rep, ad, meta.get('meta'))


def _jsonable(o):
    from .features import _jsonable as j
    return j(o)


def train_ensemble(samples, points, n_members=5, mcfg=None, lcfg=None, tcfg=None, calib_frac=0.15,
                   group_col=None, check=True, bootstrap=False, alpha=0.1, seed=0, log=None,
                   spec_kwargs=None, check_kwargs=None):
    """Train a deep ensemble.

    * the FeatureSpec and target scales are fit on the training rows only;
    * ``calib_frac`` of the samples (grouped by ``group_col``) are held out
      of training for split-conformal calibration;
    * each member uses its own inner validation split (``tcfg.val_frac``) for
      early stopping and its own seed (optionally a bootstrap resample);
    * with ``check``, physics residuals / hard bounds that the training
      targets violate in > 1 % of cases are disabled (see
      :func:`check_constraints`; ``check_kwargs`` sets its ``threshold``/``tol``).
    """
    mcfg = copy.deepcopy(mcfg or ModelConfig())
    lcfg = copy.deepcopy(lcfg or LossConfig())
    tcfg = copy.deepcopy(tcfg or TrainConfig())
    setup_torch(tcfg.threads)
    samples = samples.reset_index(drop=True)
    all_idx = np.arange(len(samples))
    fit_idx, cal_idx = inner_split(samples, all_idx, calib_frac, seed + 1000, group_col)
    s_fit = samples.iloc[fit_idx].reset_index(drop=True)
    report = None
    if check:
        report = check_constraints(s_fit, points, **(check_kwargs or {}))
        if mcfg.kind == 'physics':
            mcfg, lcfg = apply_constraint_report(report, mcfg, lcfg)
    spec = FeatureSpec(price_explicit=mcfg.price_explicit, **(spec_kwargs or {})).fit(s_fit)
    batch = spec.transform(s_fit)
    layout = make_layout(mcfg, samples.columns)
    scales = fit_target_scales(s_fit, points, layout, batch)
    dt = _dtype(tcfg.dtype)
    tb_fit = make_tensors(spec, s_fit, points, layout, batch, tcfg.max_points, dt)
    members, hists = [], []
    for mi in range(n_members):
        ms = seed + 7919 * mi
        tr, va = inner_split(s_fit, np.arange(len(s_fit)), tcfg.val_frac, ms, group_col)
        if bootstrap:
            tr = np.random.RandomState(ms).choice(tr, len(tr), replace=True)
        m, h = train_model(spec, s_fit, points, mcfg, lcfg, tcfg, layout=layout, scales=scales,
                           seed=ms, log=log, tensors=_index(tb_fit, torch.as_tensor(tr)),
                           val_tensors=_index(tb_fit, torch.as_tensor(va)) if len(va) else None)
        members.append(m)
        hists.append(h)
    ad = ApplicabilityDomain().fit(spec.full_matrix(batch), spec.full_columns())
    ens = Ensemble(members, spec, layout, scales, mcfg, lcfg, tcfg, constraint_report=report, ad=ad,
                   meta=dict(n_train=int(len(s_fit)), n_calib=int(len(cal_idx)),
                             epochs=[len(h) for h in hists], seed=seed, group_col=group_col))
    ens.histories = hists
    if len(cal_idx):
        calibrate_conformal(ens, samples.iloc[cal_idx].reset_index(drop=True), points, alpha)
    return ens


# %% Conformal calibration

def _conf_q(scores, alpha):
    """Split-conformal quantile: the ceil((n+1)(1-alpha))-th smallest score
    (the largest score if that rank exceeds n)."""
    s = np.sort(scores[np.isfinite(scores)])
    n = len(s)
    if n == 0:
        return float('nan')
    k = int(math.ceil((n + 1) * (1 - alpha)))
    return float(s[min(k, n) - 1])


def calibrate_conformal(ens, samples, points, alpha=0.1):
    """Split-conformal scaling of ensemble spreads on held-out samples:
    score = |target - mean| / (std + eps) with eps = 0.1 * training target
    std (so a single member still gets valid, constant-width intervals);
    MPSP uses relative errors with eps = 1 %. Stores quantiles in
    ``ens.conformal``."""
    batch, outs = ens.predict_samples(samples, strict=False)
    coef, g = outs['coef'], outs['g']
    sc = ens.scales
    conf = {'alpha': float(alpha), 'n': int(len(samples))}
    tgt = np.column_stack([pd.to_numeric(samples[k], errors='coerce').to_numpy(float) for k in COEFF_KEYS])
    for k, key in enumerate(COEFF_KEYS):
        mu, sd = coef[:, :, k].mean(0), coef[:, :, k].std(0)
        eps = 0.1 * sc['coef_std'][k]
        conf[key] = _conf_q(np.abs(tgt[:, k] - mu) / (sd + eps), alpha)
        conf[f'{key}_eps'] = eps
    gt = pd.to_numeric(samples['g'], errors='coerce').to_numpy(float)
    mu, sd = g.mean(0), g.std(0)
    conf['g_eps'] = 0.1 * sc['g_std']
    conf['g'] = _conf_q(np.abs(gt - mu) / (sd + conf['g_eps']), alpha)
    Y, T, M, _, mask = pad_points(samples, points)
    mm = ens.mpsp_members(outs, Y, T)
    mu, sd = mm.mean(0), mm.std(0)
    with np.errstate(divide='ignore', invalid='ignore'):
        s = np.abs(M / mu - 1.) / (sd / np.abs(mu) + 0.01)
    conf['mpsp'] = _conf_q(s[mask], alpha)
    conf['mpsp_eps_rel'] = 0.01
    ens.conformal = conf
    return conf


# %% CLI

def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description='Train a coefficient-predictor ensemble.')
    ap.add_argument('--processed', required=True, help='processed-data directory (samples/points tables)')
    ap.add_argument('--out', required=True, help='checkpoint directory')
    ap.add_argument('--members', type=int, default=5)
    ap.add_argument('--kind', default='physics', choices=('physics', 'plain'))
    ap.add_argument('--width', type=int, default=128)
    ap.add_argument('--depth', type=int, default=3)
    ap.add_argument('--dropout', type=float, default=0.05)
    ap.add_argument('--epochs', type=int, default=300)
    ap.add_argument('--batch-size', type=int, default=128)
    ap.add_argument('--lr', type=float, default=2e-3)
    ap.add_argument('--weight-decay', type=float, default=1e-4)
    ap.add_argument('--w-phys', type=float, default=0.1)
    ap.add_argument('--calib-frac', type=float, default=0.15)
    ap.add_argument('--group-col', default='config_id')
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args(argv)
    samples, points, _ = load_tables(a.processed)
    mcfg = ModelConfig(kind=a.kind, width=a.width, depth=a.depth, dropout=a.dropout)
    lcfg = LossConfig(w_phys=a.w_phys, physics=a.kind == 'physics', collocation=a.kind == 'physics')
    tcfg = TrainConfig(epochs=a.epochs, batch_size=a.batch_size, lr=a.lr, weight_decay=a.weight_decay,
                       threads=a.threads, seed=a.seed, verbose=True)
    ens = train_ensemble(samples, points, a.members, mcfg, lcfg, tcfg, a.calib_frac,
                         group_col=a.group_col or None, seed=a.seed)
    if ens.constraint_report is not None:
        print(ens.constraint_report.to_string(index=False))
    print('conformal:', {k: round(v, 3) for k, v in ens.conformal.items() if isinstance(v, float)})
    ens.save(a.out)
    print('saved', a.out)


if __name__ == '__main__':
    main()
