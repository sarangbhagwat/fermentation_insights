#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Evaluation (DESIGN.md section 6): metrics, baselines, ablations and a
cross-validation runner producing a tidy results table.

Metrics (natural units)
-----------------------
* per target (a', b', c', d', g): MAE, RMSE, R^2;
* MPSP surface: MAPE of the predicted equation on the held-out simulated TY
  points (``mpsp_mape``), the same for the per-sample OLS coefficients
  (``oracle_mape``, the irreducible error of the 4-coefficient form) and the
  mean of the samples' recorded ``fit_mape``;
* constraint-violation rates of the predictions (fraction of samples):
  ``mono_y``, ``mono_t``, ``positive`` at the sample's points, ``g_le_gmax``,
  ``nonneg_fs``/``nonneg_ng`` (price-explicit models), and ``mono_r`` by a
  finite difference in log r (model-agnostic);
* interval coverage and mean relative width (ensembles with conformal
  calibration).

Methods
-------
``pinn`` (physics-informed ensemble), ``plain`` (plain MLP), ``ridge``,
``hgb`` (sklearn HistGradientBoosting per target), ``knn`` (nearest-neighbour
lookup in the standardized feature space, k = 5, distance-weighted) and
ablations ``pinn-<name>`` with name in ``ABLATIONS``.

CLI::

    python -m fermentation_insights.coeff_predictor.evaluate --processed DIR \\
        --scheme product --methods pinn,plain,ridge,hgb,knn --members 3 \\
        --ablations no_phys,no_channels --out results.csv
"""
import copy
import math
import time
import warnings

import numpy as np
import pandas as pd

from .features import FeatureSpec
from .model import ModelConfig, COEFF_KEYS, mpsp_equation_np
from .train import (LossConfig, TrainConfig, train_ensemble, make_splits, pad_points,
                    constraint_case_values, load_tables, ensemble_splits, SCHEMES, SCHEME_GROUP,
                    setup_torch)

__all__ = ('compute_metrics', 'fit_method', 'run_cv', 'summarize', 'grouped_holdout', 'ABLATIONS',
           'METHODS', 'main')

TARGETS = COEFF_KEYS + ('g',)
METHODS = ('pinn', 'plain', 'ridge', 'hgb', 'knn')

#: ablation name -> (ModelConfig overrides, LossConfig overrides)
ABLATIONS = {
    'no_channels': (dict(channels=False), dict(w_ch=0.)),
    'no_price': (dict(price_explicit=False), {}),
    'no_struct': (dict(structured_fs=False, g_bound=False), {}),
    'no_phys': ({}, dict(physics=False, collocation=False)),
    'no_colloc': ({}, dict(collocation=False)),
    'no_chloss': ({}, dict(w_ch=0.)),
    'no_mpsp_loss': ({}, dict(w_mpsp=0.)),
}


# %% Methods

class _limit_threads:
    """Limit BLAS/OpenMP threads of sklearn (OpenMP spin-waiting is very slow
    on the shared, oversubscribed machine)."""

    def __init__(self, n=1):
        self.n, self.ctx = n, None

    def __enter__(self):
        try:
            from threadpoolctl import threadpool_limits
            self.ctx = threadpool_limits(self.n)
        except Exception:  # pragma: no cover
            self.ctx = None
        return self

    def __exit__(self, *exc):
        if self.ctx is not None:
            self.ctx.restore_original_limits()


def _targets(samples):
    return np.column_stack([pd.to_numeric(samples[k], errors='coerce').to_numpy(float) for k in TARGETS])


def _n_configs(df):
    return int(df['config_id'].astype(str).nunique()) if 'config_id' in df else int(len(df))


class _EnsembleMethod:
    def __init__(self, ens):
        self.ens = ens
        m = ens.meta
        self.fit_info = dict(n_fit=m.get('n_fit', m.get('n_train')), n_fit_configs=m.get('n_fit_configs'),
                             n_calib=m.get('n_calib'), calib_group_used=m.get('calib_group_used'))

    def predict(self, df):
        batch, outs = self.ens.predict_samples(df, strict=False)
        return dict(coef=outs['coef'].mean(0), g=outs['g'].mean(0), members=outs,
                    conformal=self.ens.conformal, g_max=batch.g_max, invalid=batch.invalid,
                    intensities=outs.get('intensities'), price_keys=self.ens.layout['price_keys'])


class _SklearnMethod:
    def __init__(self, kind, samples, seed=0):
        self.spec = FeatureSpec(price_explicit=False).fit(samples)
        X = self.spec.full_matrix(self.spec.transform(samples, strict=False))
        Z = _targets(samples)
        ok = np.isfinite(Z).all(1) & np.isfinite(X).all(1)
        X, Z = X[ok], Z[ok]
        self.mu, self.sd = Z.mean(0), Z.std(0) + 1e-12
        Zs = (Z - self.mu) / self.sd
        self.kind = kind
        with _limit_threads(1):
            self._fit(kind, X, Zs, seed)

    def _fit(self, kind, X, Zs, seed):
        from sklearn.linear_model import RidgeCV
        from sklearn.ensemble import HistGradientBoostingRegressor
        from sklearn.neighbors import KNeighborsRegressor
        if kind == 'ridge':
            self.models = [RidgeCV(alphas=np.logspace(-3, 3, 13)).fit(X, Zs)]
        elif kind == 'knn':
            self.models = [KNeighborsRegressor(n_neighbors=min(5, len(X)), weights='distance').fit(X, Zs)]
        elif kind == 'hgb':
            self.models = [HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_leaf_nodes=15,
                                                         min_samples_leaf=min(10, max(2, len(X) // 20)),
                                                         l2_regularization=1e-3, random_state=seed).fit(X, Zs[:, j])
                           for j in range(Zs.shape[1])]
        else:
            raise ValueError(kind)

    def predict(self, df):
        b = self.spec.transform(df, strict=False)
        X = np.nan_to_num(self.spec.full_matrix(b))
        with _limit_threads(1):
            if len(self.models) == 1:
                Zs = self.models[0].predict(X)
            else:
                Zs = np.column_stack([m.predict(X) for m in self.models])
        Z = Zs * self.sd + self.mu
        return dict(coef=Z[:, :4], g=Z[:, 4], g_max=b.g_max, invalid=b.invalid)


def _configs(method, mcfg_kw=None, lcfg_kw=None):
    base, _, abl = method.partition('-')
    mk = dict(mcfg_kw or {})
    lk = dict(lcfg_kw or {})
    if base == 'plain':
        mk['kind'] = 'plain'
        lk.update(physics=False, collocation=False, w_ch=0.)
    if abl:
        for a in abl.split('+'):
            if a not in ABLATIONS:
                raise ValueError(f'unknown ablation {a!r}; choose from {sorted(ABLATIONS)}')
            m_o, l_o = ABLATIONS[a]
            mk.update(m_o)
            lk.update(l_o)
    return ModelConfig(**mk), LossConfig(**lk)


def fit_method(method, samples, points, n_members=1, mcfg_kw=None, lcfg_kw=None, tcfg_kw=None,
               group_col=None, calib_frac=0.15, seed=0, log=None):
    """Fit ``method`` on training tables; returns an object with
    ``predict(samples_df) -> dict(coef, g, ...)`` and ``fit_info`` (samples /
    configurations actually fitted). The sklearn baselines are fit on exactly
    the samples the network members train on (same calibration and
    early-stopping hold-outs, :func:`train.ensemble_splits`), so that methods
    are compared at equal training data."""
    base = method.partition('-')[0]
    tcfg = TrainConfig(**(tcfg_kw or {}))
    if base in ('ridge', 'hgb', 'knn'):
        samples = samples.reset_index(drop=True)
        sp = ensemble_splits(samples, n_members, calib_frac, tcfg.val_frac, group_col, seed)
        s_fit = samples.iloc[sp['fit']].reset_index(drop=True)
        s_use = s_fit.iloc[sorted({int(i) for tr, _, _ in sp['members'] for i in tr})].reset_index(drop=True)
        m = _SklearnMethod(base, s_use, seed)
        m.fit_info = dict(n_fit=len(s_use), n_fit_configs=_n_configs(s_use), n_calib=len(sp['calib']),
                          calib_group_used=sp['calib_group'])
        return m
    if base not in ('pinn', 'plain'):
        raise ValueError(f'unknown method {method!r}')
    mcfg, lcfg = _configs(method, mcfg_kw, lcfg_kw)
    ens = train_ensemble(samples, points, n_members, mcfg, lcfg, tcfg, calib_frac=calib_frac,
                         group_col=group_col, seed=seed, log=log)
    return _EnsembleMethod(ens)


# %% Metrics

def _r2(y, yhat):
    m = np.isfinite(y) & np.isfinite(yhat)
    if m.sum() < 2:
        return float('nan')
    ss = np.sum((y[m] - y[m].mean()) ** 2)
    return float(1 - np.sum((y[m] - yhat[m]) ** 2) / ss) if ss > 0 else float('nan')


def compute_metrics(pred, samples, points, model=None, fd_delta=0.05, tol=1e-3):
    """Metrics of predictions ``pred`` (output of ``method.predict``) on the
    held-out ``samples``/``points``. ``model`` (the fitted method) enables the
    finite-difference ``mono_r`` check. Returns a list of
    ``(metric, target, value)``."""
    rows = []
    if pred.get('invalid') is not None:
        rows.append(('n_invalid', 'rows', float(np.sum(pred['invalid']))))
    Z = _targets(samples)
    P = np.column_stack([pred['coef'], pred['g']])
    for j, k in enumerate(TARGETS):
        m = np.isfinite(Z[:, j]) & np.isfinite(P[:, j])
        e = P[m, j] - Z[m, j]
        rows += [('MAE', k, float(np.mean(np.abs(e))) if m.any() else np.nan),
                 ('RMSE', k, float(np.sqrt(np.mean(e ** 2))) if m.any() else np.nan),
                 ('R2', k, _r2(Z[:, j], P[:, j]))]
    Y, T, M, _, mask = pad_points(samples, points)
    Mh = mpsp_equation_np(pred['coef'][:, None, :], pred['g'][:, None], Y, T)
    Mo = mpsp_equation_np(Z[:, None, :4], Z[:, 4:5], Y, T)
    with np.errstate(divide='ignore', invalid='ignore'):
        ape = np.abs(Mh - M) / np.abs(M)
        ape_o = np.abs(Mo - M) / np.abs(M)
    v, vo = ape[mask], ape_o[mask]
    rows.append(('mpsp_mape', 'MPSP', float(np.nanmean(v)) if len(v) else np.nan))
    rows.append(('mpsp_mape_median', 'MPSP', float(np.nanmedian(v)) if len(v) else np.nan))
    rows.append(('oracle_mape', 'MPSP', float(np.nanmean(vo)) if len(vo) else np.nan))
    if 'fit_mape' in samples:
        rows.append(('fit_mape', 'MPSP', float(np.nanmean(pd.to_numeric(samples['fit_mape'], errors='coerce')))))
    # constraint violations of the predictions
    qfs = qng = None
    if pred.get('intensities') is not None:
        qi = pred['intensities'].mean(0)
        pk = pred.get('price_keys', [])
        qfs = qi[:, pk.index('fs')] if 'fs' in pk else None
        qng = qi[:, pk.index('ng')] if 'ng' in pk else None
    cv = constraint_case_values(pred['coef'], pred['g'], Y, T, mask, None, qfs, qng, tol)
    for k, val in cv.items():
        rows.append(('violation_rate', k, float(np.nanmean(1 - val)) if np.isfinite(val).any() else np.nan))
    gm = pred.get('g_max')
    if gm is not None:
        ok = np.isfinite(gm) & np.isfinite(pred['g'])
        rows.append(('violation_rate', 'g_le_gmax', float(np.mean(pred['g'][ok] > gm[ok] * (1 + 1e-6)))
                     if ok.any() else np.nan))
    if model is not None and 'r' in samples:
        lo, hi = samples.copy(), samples.copy()
        r = pd.to_numeric(samples['r'], errors='coerce')
        lo['r'], hi['r'] = r * math.exp(-fd_delta), r * math.exp(fd_delta)
        pl, ph = model.predict(lo), model.predict(hi)
        Ml = mpsp_equation_np(pl['coef'][:, None, :], pl['g'][:, None], Y, T)
        Mh2 = mpsp_equation_np(ph['coef'][:, None, :], ph['g'][:, None], Y, T)
        viol = np.where(mask, Mh2 > Ml + tol * np.abs(Ml), False).any(1)
        rows.append(('violation_rate', 'mono_r', float(viol[mask.any(1)].mean()) if mask.any() else np.nan))
    # interval coverage
    conf, outs = pred.get('conformal') or {}, pred.get('members')
    if outs is not None and conf:
        for j, k in enumerate(TARGETS):
            q = conf.get(k)
            if q is None or not np.isfinite(q):
                continue
            arr = outs['coef'][:, :, j] if j < 4 else outs['g']
            mu, sd = arr.mean(0), arr.std(0)
            hw = q * (sd + conf[f'{k}_eps'])
            m = np.isfinite(Z[:, j])
            rows.append(('coverage', k, float(np.mean(np.abs(Z[m, j] - mu[m]) <= hw[m]))))
            rows.append(('interval_halfwidth', k, float(np.mean(hw[m]))))
        q = conf.get('mpsp')
        if q is not None and np.isfinite(q):
            mm = mpsp_equation_np(outs['coef'][:, :, None, :], outs['g'][:, :, None], Y[None], T[None])
            mu, sd = mm.mean(0), mm.std(0)
            hw = q * (sd + conf.get('mpsp_eps_rel', 0.01) * np.abs(mu))
            cov = (np.abs(M - mu) <= hw)[mask]
            rows.append(('coverage', 'MPSP', float(cov.mean()) if len(cov) else np.nan))
            with np.errstate(divide='ignore', invalid='ignore'):
                rows.append(('interval_halfwidth_rel', 'MPSP', float(np.nanmean((hw / np.abs(mu))[mask]))))
    return rows


# %% CV runner

def run_cv(samples, points, scheme='random', methods=METHODS, n_splits=5, max_folds=None, n_members=1,
           mcfg_kw=None, lcfg_kw=None, tcfg_kw=None, seed=0, log=print, calib_frac=0.15):
    """Cross-validate ``methods`` under splitting ``scheme``. Feature and
    target scaling are fit inside each training fold; inner calibration /
    early-stopping splits are grouped like the scheme (``SCHEME_GROUP``).
    Splits with an empty train or test set are skipped (logged). Returns a
    tidy DataFrame with columns ``scheme, fold, label, method, n_train,
    n_test, n_fit, n_fit_configs, n_calib, calib_group, seconds, metric,
    target, value``
    (``n_train`` = fold size; ``n_fit``/``n_fit_configs`` = samples /
    configurations actually fitted; failures are ``metric='error'`` rows).
    Pooled out-of-fold rows (fold -1) use only the folds in which every
    method succeeded."""
    samples = samples.reset_index(drop=True)
    splits = make_splits(samples, scheme, n_splits, seed)
    if max_folds:
        splits = splits[:max_folds]
    group_col = SCHEME_GROUP.get(scheme)
    rows, oof = [], {}
    ok_folds = {m: set() for m in methods}
    for f, (tr, te, label) in enumerate(splits):
        if not len(tr) or not len(te):
            if log:
                log(f'[{scheme} fold {f} ({label})] skipped: empty {"train" if not len(tr) else "test"} set')
            continue
        s_tr, s_te = samples.iloc[tr].reset_index(drop=True), samples.iloc[te].reset_index(drop=True)
        p_tr = points[points['sample_id'].astype(str).isin(set(s_tr['sample_id'].astype(str)))]
        p_te = points[points['sample_id'].astype(str).isin(set(s_te['sample_id'].astype(str)))]
        for meth in methods:
            t0 = time.time()
            info = {}
            try:
                model = fit_method(meth, s_tr, p_tr, n_members, mcfg_kw, lcfg_kw, tcfg_kw,
                                   group_col=group_col, calib_frac=calib_frac, seed=seed + f)
                info = getattr(model, 'fit_info', {}) or {}
                pred = model.predict(s_te)
                mets = compute_metrics(pred, s_te, p_te, model)
                oof.setdefault(meth, {})[f] = (pred, s_te, p_te)
                ok_folds[meth].add(f)
            except Exception as e:  # keep the CV going; record the failure
                mets = [('error', type(e).__name__ + ': ' + str(e)[:200], np.nan)]
                if log:
                    log(f'[{scheme} fold {f} ({label})] {meth} FAILED: {type(e).__name__}: {e}')
            dt = time.time() - t0
            for metric, target, value in mets:
                rows.append(dict(scheme=scheme, fold=f, label=label, method=meth, n_train=len(tr),
                                 n_test=len(te), n_fit=info.get('n_fit', np.nan),
                                 n_fit_configs=info.get('n_fit_configs', np.nan),
                                 n_calib=info.get('n_calib', np.nan),
                                 calib_group=info.get('calib_group_used'), seconds=dt, metric=metric,
                                 target=target, value=value))
            if log:
                mm = {(m, t): v for m, t, v in mets}
                log(f'[{scheme} fold {f} ({label})] {meth:18s} mpsp_mape={mm.get(("mpsp_mape", "MPSP"), np.nan):.4f} '
                    f'oracle={mm.get(("oracle_mape", "MPSP"), np.nan):.4f} R2(a,b,c,d,g)='
                    + ','.join(f'{mm.get(("R2", k), np.nan):.3f}' for k in TARGETS) + f' ({dt:.0f}s)')
    # pooled out-of-fold metrics (R^2 within a single held-out group is not
    # meaningful for leave-one-group-out schemes), over the common folds
    common = set.intersection(*ok_folds.values()) if ok_folds else set()
    if any(ok_folds[m] != common for m in ok_folds):
        warnings.warn(f'run_cv({scheme}): methods failed in different folds; pooled metrics use only the '
                      f'{len(common)} fold(s) in which every method succeeded', UserWarning, stacklevel=2)
    if len(splits) > 1:
        for meth, parts in oof.items():
            parts = [parts[f] for f in sorted(common)]
            if not parts:
                continue
            pred = {k: np.concatenate([p[0][k] for p in parts]) for k in ('coef', 'g', 'g_max')}
            s_all = pd.concat([p[1] for p in parts], ignore_index=True)
            p_all = pd.concat([p[2] for p in parts], ignore_index=True)
            for metric, target, value in compute_metrics(pred, s_all, p_all):
                if metric in ('MAE', 'RMSE', 'R2', 'mpsp_mape', 'mpsp_mape_median', 'oracle_mape', 'fit_mape'):
                    rows.append(dict(scheme=scheme, fold=-1, label='pooled', method=meth, n_train=np.nan,
                                     n_test=len(s_all), seconds=np.nan, metric=metric, target=target,
                                     value=value))
    return pd.DataFrame(rows)


def summarize(results, metrics=('mpsp_mape', 'oracle_mape', 'R2', 'MAE', 'violation_rate', 'coverage'),
              common_folds=True):
    """Mean over folds per (method, metric, target); R^2 and MAE are the
    pooled out-of-fold values for the (scheme, method) pairs that have them.
    With ``common_folds`` the per-fold means use only the folds in which
    every method of the scheme succeeded (a warning lists exclusions). Adds
    the row ``n_folds_ok`` (folds without error) per scheme and method."""
    res = results
    err = res[res['metric'] == 'error']
    folds = res[res['fold'] != -1][['scheme', 'method', 'fold']].drop_duplicates()
    failed = set(map(tuple, err[['scheme', 'method', 'fold']].to_numpy()))
    folds = folds.assign(ok=[(s, m, f) not in failed for s, m, f in folds.to_numpy()])
    n_ok = folds.groupby(['scheme', 'method'])['ok'].sum().astype(float)
    r = res[res['metric'].isin(metrics)]
    pooled = r[r['fold'] == -1]
    per = r[r['fold'] != -1]
    if common_folds and len(failed):
        bad = {(s, f) for s, m, f in failed}
        drop = np.array([(s, f) in bad for s, f in zip(per['scheme'], per['fold'])], bool)
        if drop.any():
            warnings.warn(f'summarize: excluded fold(s) {sorted(bad)} where some method failed', UserWarning,
                          stacklevel=2)
            per = per[~drop]
    if len(pooled):
        keys = set(zip(pooled['scheme'], pooled['method']))
        has = np.array([(s, m) in keys for s, m in zip(per['scheme'], per['method'])], bool)
        per = per[~(per['metric'].isin(['R2', 'MAE']).to_numpy() & has)]
        per = pd.concat([per, pooled[pooled['metric'].isin(['R2', 'MAE'])]])
    out = per.groupby(['scheme', 'method', 'metric', 'target'])['value'].mean()
    nk = pd.Series(n_ok.to_numpy(), index=pd.MultiIndex.from_tuples(
        [(s, m, 'n_folds_ok', 'folds') for s, m in n_ok.index], names=out.index.names))
    return pd.concat([out, nk]).unstack('method')


def grouped_holdout(samples, frac, seed=12345, col='config_id'):
    """Strict grouped final hold-out: ``round(frac * n_groups)`` whole groups
    (at least 1, at least 1 kept). Raises ValueError with fewer than 2
    groups (a random split would not be a hold-out of configurations)."""
    g = samples[col].astype(str).to_numpy()
    u = np.array(sorted(set(g)))
    if len(u) < 2:
        raise ValueError(f'a grouped final test set needs >= 2 distinct {col!r} values (found {len(u)})')
    k = min(max(1, int(round(frac * len(u)))), len(u) - 1)
    held = np.random.RandomState(seed).choice(u, k, replace=False)
    te = np.isin(g, held)
    return np.flatnonzero(~te), np.flatnonzero(te)


# %% CLI

def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description='Cross-validate the coefficient predictor and baselines.')
    ap.add_argument('--processed', required=True, help='processed-data directory')
    ap.add_argument('--scheme', default='random', choices=sorted(SCHEMES))
    ap.add_argument('--methods', default='pinn,plain,ridge,hgb,knn')
    ap.add_argument('--ablations', default='', help=f'comma list from {sorted(ABLATIONS)} (pinn-<name>)')
    ap.add_argument('--folds', type=int, default=5, help='K for the random scheme')
    ap.add_argument('--max-folds', type=int, default=0, help='evaluate only the first N folds (0 = all)')
    ap.add_argument('--members', type=int, default=1)
    ap.add_argument('--width', type=int, default=128)
    ap.add_argument('--depth', type=int, default=3)
    ap.add_argument('--dropout', type=float, default=0.05)
    ap.add_argument('--epochs', type=int, default=300)
    ap.add_argument('--batch-size', type=int, default=128)
    ap.add_argument('--lr', type=float, default=2e-3)
    ap.add_argument('--weight-decay', type=float, default=1e-4)
    ap.add_argument('--patience', type=int, default=40)
    ap.add_argument('--w-phys', type=float, default=0.1)
    ap.add_argument('--w-ch', type=float, default=0.5)
    ap.add_argument('--holdout-frac', type=float, default=0.,
                    help='set aside a final test set of whole configurations before CV '
                         '(needs >= 2 configurations)')
    ap.add_argument('--report-final', action='store_true',
                    help='also train on all CV data and report once on the final test set')
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', default='cv_results.csv')
    a = ap.parse_args(argv)
    if a.report_final and not a.holdout_frac > 0:
        ap.error('--report-final requires --holdout-frac > 0')
    setup_torch(a.threads)
    samples, points, _ = load_tables(a.processed)
    samples = samples.reset_index(drop=True)
    final = None
    if a.holdout_frac > 0:
        try:
            cv_idx, te_idx = grouped_holdout(samples, a.holdout_frac, 12345, 'config_id')
        except ValueError as e:
            ap.error(f'--holdout-frac: {e}')
        final = samples.iloc[te_idx].reset_index(drop=True)
        samples = samples.iloc[cv_idx].reset_index(drop=True)
    methods = [m for m in a.methods.split(',') if m]
    methods += [f'pinn-{x}' for x in a.ablations.split(',') if x]
    mk = dict(width=a.width, depth=a.depth, dropout=a.dropout)
    lk = dict(w_phys=a.w_phys, w_ch=a.w_ch)
    tk = dict(epochs=a.epochs, batch_size=a.batch_size, lr=a.lr, weight_decay=a.weight_decay,
              patience=a.patience, threads=a.threads, seed=a.seed)
    res = run_cv(samples, points, a.scheme, methods, a.folds, a.max_folds or None, a.members, mk, lk, tk, a.seed)
    if final is not None and a.report_final:
        p_tr = points[points['sample_id'].astype(str).isin(set(samples['sample_id'].astype(str)))]
        p_te = points[points['sample_id'].astype(str).isin(set(final['sample_id'].astype(str)))]
        rows = []
        for meth in methods:
            model = fit_method(meth, samples, p_tr, a.members, mk, lk, tk, group_col='config_id', seed=a.seed)
            for metric, target, value in compute_metrics(model.predict(final), final, p_te, model):
                info = getattr(model, 'fit_info', {}) or {}
                rows.append(dict(scheme='final_test', fold=0, label='final', method=meth, n_train=len(samples),
                                 n_test=len(final), n_fit=info.get('n_fit', np.nan),
                                 n_fit_configs=info.get('n_fit_configs', np.nan),
                                 n_calib=info.get('n_calib', np.nan), seconds=np.nan, metric=metric,
                                 target=target, value=value))
        res = pd.concat([res, pd.DataFrame(rows)], ignore_index=True)
    res.to_csv(a.out, index=False)
    with pd.option_context('display.width', 200, 'display.max_columns', 20, 'display.max_rows', 200):
        print(summarize(res).round(4).to_string())
    print('wrote', a.out)


if __name__ == '__main__':
    main()
