#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Report-only statistics of the experiment pipeline (``experiments.py``),
computed from its saved outputs; none of them is used to fit, tune, select
or calibrate anything:

* :func:`ablation_stats` -- paired fold-level statistics of the ablations
  (and of the plain MLP) against the full PI-NN;
* :func:`paired_test_comparison` -- final-test MAPE of every method with a
  configuration-cluster bootstrap CI of the paired difference to the PI-NN;
* :func:`shipped_quantile_coverage` -- coverage of the shipped interval
  quantiles on the out-of-fold scores of each CV scheme (new feedstocks /
  products);
* :func:`ad_validation` -- applicability-domain flag rates and errors in and
  out of the domain on the test sets;
* :func:`raw_sample_coverage` -- planned / completed / retained samples per
  configuration.
"""
import json
import math

import numpy as np
import pandas as pd

__all__ = ('ablation_stats', 'paired_test_comparison', 'shipped_quantile_coverage', 'ad_validation',
           'raw_sample_coverage', 'sample_ape')


def _paired(d):
    """Statistics of paired fold differences ``d`` (Series indexed by fold label)."""
    d = d.dropna().astype(float)
    n = len(d)
    out = dict(n_folds=n, mean=float(d.mean()) if n else math.nan, median=float(d.median()) if n else math.nan)
    se = float(d.std(ddof=1) / math.sqrt(n)) if n > 1 else math.nan
    out.update(se=se, t=out['mean'] / se if se and np.isfinite(se) and se > 0 else math.nan,
               n_variant_better=int((d < 0).sum()), n_variant_worse=int((d > 0).sum()))
    p = math.nan
    if n >= 5 and (d != 0).any():
        try:
            from scipy.stats import wilcoxon
            p = float(wilcoxon(d.to_numpy()).pvalue)
        except Exception:  # pragma: no cover
            pass
    out['p_wilcoxon'] = p
    if n:
        k = d.abs().idxmax()
        out.update(dominant_fold=str(k), dominant_delta=float(d[k]),
                   mean_without_dominant=float(d.drop(k).mean()) if n > 1 else math.nan)
    return out


def ablation_stats(res, schemes=('config', 'product'), ablations=('pinn-no_channels', 'pinn-no_price',
                                                                 'pinn-no_phys'), reference='pinn'):
    """Paired fold-level statistics of (variant - full PI-NN) fold MPSP MAPE
    (fractions) per CV scheme, from the tidy CV results (``cv_results.csv``):
    every ablation on the folds where it ran; the plain MLP on all folds and,
    when the ablations ran on a subset, on that subset too (``subset``).
    Returns (table, per-fold differences)."""
    r = res[(res['fold'] >= 0) & (res['metric'] == 'mpsp_mape') & (res['target'] == 'MPSP')]
    rows, per = [], []
    for sch in schemes:
        g = r[r['scheme'] == sch].pivot_table(index='label', columns='method', values='value', aggfunc='mean')
        if reference not in g:
            continue
        sub = None
        for m in list(ablations) + ['plain']:
            if m not in g:
                continue
            ok = g[[m, reference]].dropna()
            if not len(ok):
                continue
            d = ok[m] - ok[reference]
            rows.append(dict(scheme=sch, variant=m, subset='all folds of the variant', **_paired(d)))
            per += [dict(scheme=sch, variant=m, subset='all folds of the variant', fold=str(k), delta=float(v))
                    for k, v in d.items()]
            if m in ablations and len(ok) < len(g[reference].dropna()):
                sub = ok.index
        if sub is not None and 'plain' in g:
            ok = g.loc[g.index.isin(sub), ['plain', reference]].dropna()
            d = ok['plain'] - ok[reference]
            rows.append(dict(scheme=sch, variant='plain', subset='ablation folds', **_paired(d)))
            per += [dict(scheme=sch, variant='plain', subset='ablation folds', fold=str(k), delta=float(v))
                    for k, v in d.items()]
    return pd.DataFrame(rows), pd.DataFrame(per)


def baseline_fold_stats(res, schemes=('config', 'feedstock', 'product'), methods=('knn', 'plain', 'ridge', 'hgb'),
                        reference='pinn'):
    """Paired fold-level statistics of (baseline - PI-NN) fold MPSP MAPE
    (fractions; > 0: the PI-NN is better) per CV scheme and baseline, over
    the folds where both ran (``cv_results.csv``): mean, s.e., t, median,
    folds where either is better, Wilcoxon signed-rank p (>= 5 folds), the
    two-sided sign-test p and the dominant fold."""
    r = res[(res['fold'] >= 0) & (res['metric'] == 'mpsp_mape') & (res['target'] == 'MPSP')]
    rows = []
    for sch in schemes:
        g = r[r['scheme'] == sch].pivot_table(index='label', columns='method', values='value', aggfunc='mean')
        if reference not in g:
            continue
        for m in methods:
            if m not in g:
                continue
            ok = g[[m, reference]].dropna()
            if not len(ok):
                continue
            st = _paired(ok[m] - ok[reference])
            nb, nw = st.pop('n_variant_better'), st.pop('n_variant_worse')
            p_sign = math.nan
            if nb + nw:
                try:
                    from scipy.stats import binomtest
                    p_sign = float(binomtest(nw, nb + nw, 0.5).pvalue)
                except Exception:  # pragma: no cover
                    pass
            rows.append(dict(scheme=sch, method=m, n_folds=st.pop('n_folds'), mean=st['mean'], se=st['se'],
                             t=st['t'], median=st['median'], n_pinn_better=nw, n_baseline_better=nb,
                             p_wilcoxon=st['p_wilcoxon'], p_sign=p_sign, dominant_fold=st.get('dominant_fold'),
                             dominant_delta=st.get('dominant_delta'),
                             mean_without_dominant=st.get('mean_without_dominant')))
    return pd.DataFrame(rows)


def family_coverage(scores, conf, scheme='config', method='pinn', keys=('g', 'mpsp'), families=None):
    """Coverage of the interval quantiles ``conf`` on the out-of-fold scores
    of one CV scheme per biorefinery family (pooled over its held-out
    samples) and per held-out configuration: (per family, per configuration)
    tables. The pooled coverage can sit at the nominal level while whole
    families are over- or under-covered."""
    from fermentation_insights.coeff_predictor.train import rescore_conformal, CONF_KEYS
    eps = {k: conf.get(f'{k}_eps') for k in CONF_KEYS if k != 'mpsp' and conf.get(f'{k}_eps') is not None}
    eps['mpsp_rel'] = conf.get('mpsp_eps_rel', 0.01)
    g = scores[(scores['scheme'] == scheme) & (scores['method'] == method)]
    if not len(g):
        return pd.DataFrame(), pd.DataFrame()
    sc, _ = rescore_conformal(g, eps)
    cid = sc['config_id'].astype(str)
    fam = cid.map(families) if families else cid.str.split('_').str[0]
    fam = fam.where(fam.notna(), cid)
    fr, cr = [], []
    for k in keys:
        q = conf.get(k)
        if q is None or not np.isfinite(q) or k not in sc:
            continue
        v = sc[k].to_numpy(float)
        ok = np.isfinite(v)
        cov = pd.Series(v[ok] <= q, index=sc.index[ok], dtype=float)
        byc = cov.groupby(sc.loc[ok, 'config_id'].astype(str)).mean()
        fr += [dict(key=k, family=str(f), n=int(len(c)), n_configs=int(sc.loc[c.index, 'config_id'].nunique()),
                    coverage=float(c.mean()))
               for f, c in cov.groupby(fam[ok].astype(str))]
        cr += [dict(key=k, config_id=str(c), coverage=float(v_)) for c, v_ in byc.items()]
    return pd.DataFrame(fr), pd.DataFrame(cr)


def sample_ape(ev, samples_by_id, points_by_id):
    """Per held-out sample of an evaluation dict (``ids``, ``coef``, ``g``):
    ``sample_id``, ``config_id``, sum of the absolute percentage errors of
    MPSP over its simulated TY points (``ape_sum``) and their number (``n``)."""
    from fermentation_insights.coeff_predictor.train import pad_points
    from fermentation_insights.coeff_predictor.model import mpsp_equation_np
    s = samples_by_id.loc[ev['ids']].reset_index(drop=True)
    p = pd.concat([points_by_id[i] for i in ev['ids'] if i in points_by_id], ignore_index=True)
    Y, T, M, _, mask = pad_points(s, p)
    Mh = mpsp_equation_np(np.asarray(ev['coef'], float)[:, None, :], np.asarray(ev['g'], float)[:, None], Y, T)
    with np.errstate(divide='ignore', invalid='ignore'):
        ape = np.where(mask, np.abs(Mh - M) / np.abs(M), 0.)
    ape = np.where(np.isfinite(ape), ape, 0.)
    return pd.DataFrame(dict(sample_id=s['sample_id'].astype(str), config_id=s['config_id'].astype(str),
                             ape_sum=ape.sum(1), n=mask.sum(1)))


def paired_test_comparison(final, samples_by_id, points_by_id, sets=('test_ood', 'test_id'), reference='pinn',
                           n_boot=2000, seed=0, level=0.95):
    """Final-test MPSP MAPE (pooled over the TY points) of every method and
    the paired difference (reference - method) with a percentile bootstrap CI
    that resamples the test CONFIGURATIONS (cluster bootstrap; the same
    draws for every method). ``final``: {method: fit result} of the
    pipeline's final stage. Report only (the test split is never used to
    choose anything)."""
    rng = np.random.default_rng(seed)
    rows = []
    for st in sets:
        per = {}
        for m, r in (final or {}).items():
            if r is None or r.get('error') or st not in r.get('evals', {}):
                continue
            a = sample_ape(r['evals'][st], samples_by_id, points_by_id)
            per[m] = a.groupby('config_id')[['ape_sum', 'n']].sum()
        if reference not in per:
            continue
        cfg = per[reference].index.to_numpy()
        draws = rng.integers(0, len(cfg), size=(n_boot, len(cfg)))

        def boot(t):
            t = t.reindex(cfg)
            s, n = t['ape_sum'].to_numpy(float), t['n'].to_numpy(float)
            return s[draws].sum(1) / n[draws].sum(1), s.sum() / n.sum()

        ref_b, ref = boot(per[reference])
        for m, t in per.items():
            mb, mv = boot(t)
            row = dict(set=st, method=m, n_configs=len(cfg), n_samples=None, mpsp_mape=mv)
            if m != reference:
                d = ref_b - mb
                lo, hi = np.percentile(d, [50 * (1 - level), 50 * (1 + level)])
                row.update(diff_reference_minus_method=ref - mv, ci_lo=float(lo), ci_hi=float(hi),
                           frac_boot_reference_better=float((d < 0).mean()))
            rows.append(row)
    return pd.DataFrame(rows)


def shipped_quantile_coverage(scores, conf, schemes=('config', 'feedstock', 'product'), method='pinn'):
    """Coverage of the interval quantiles ``conf`` (an ensemble's
    ``conformal`` dict, e.g. of ``models/release``) on the out-of-fold scores
    (``cv_oof_scores.csv.gz``) of each CV scheme, re-scored with that
    ensemble's spread floors (``train.rescore_conformal``): per scheme x key,
    pooled over the held-out samples, mean / min over the held-out groups
    (folds) and the worst group. For the LOCO (``config``) scheme the
    quantiles were computed from these very scores (in-sample: about the
    nominal level by construction; the cross-fitted estimate is the honest
    one); for leave-one-feedstock-out and leave-one-product-out it is the
    coverage the shipped intervals would have for a new feedstock / product."""
    from fermentation_insights.coeff_predictor.train import rescore_conformal, CONF_KEYS
    eps = {k: conf.get(f'{k}_eps') for k in CONF_KEYS if k != 'mpsp' and conf.get(f'{k}_eps') is not None}
    eps['mpsp_rel'] = conf.get('mpsp_eps_rel', 0.01)
    rows = []
    for sch in schemes:
        g = scores[(scores['scheme'] == sch) & (scores['method'] == method)]
        if not len(g):
            continue
        sc, _ = rescore_conformal(g, eps)
        for k in CONF_KEYS:
            q = conf.get(k)
            if q is None or not np.isfinite(q) or k not in sc:
                continue
            v = sc[k].to_numpy(float)
            ok = np.isfinite(v)
            cov = pd.Series(v[ok] <= q, index=sc.index[ok])
            by = cov.groupby(sc.loc[ok, 'fold'].astype(str)).mean()
            rows.append(dict(scheme=sch, key=k, q=float(q), n=int(ok.sum()), n_groups=int(len(by)),
                             pooled=float(cov.mean()), group_mean=float(by.mean()), group_min=float(by.min()),
                             worst_group=str(by.idxmin()), in_sample=sch == 'config'))
    return pd.DataFrame(rows)


def ad_validation(P, evals, samples_by_id, points_by_id):
    """Applicability domain of a checkpoint (``CoefficientPredictor``) on
    evaluation sets {name: eval dict of that checkpoint}: fraction of samples
    in the domain, MPSP MAPE in / out of it, and the Spearman correlation
    between ``ad_ratio`` and the sample's MAPE. Features as in training
    (``spec.transform``), so the flags are those of the validated path.
    Returns (summary, per-configuration table)."""
    from scipy.stats import spearmanr
    rows, per = [], []
    for name, ev in evals.items():
        s = samples_by_id.loc[ev['ids']].reset_index(drop=True)
        batch = P.spec.transform(s, strict=False)
        ratio, oor = P._ad(batch)
        ind = np.isfinite(ratio) & (ratio <= 1.) & np.array([not o for o in oor])
        a = sample_ape(ev, samples_by_id, points_by_id)
        sm = (a['ape_sum'] / a['n'].clip(lower=1)).to_numpy()

        def mape(m):
            return float(a.loc[m, 'ape_sum'].sum() / max(a.loc[m, 'n'].sum(), 1)) if m.any() else math.nan

        rho = spearmanr(ratio, sm).statistic if np.isfinite(ratio).all() and len(s) > 2 else math.nan
        rows.append(dict(set=name, n_samples=len(s), frac_in_domain=float(ind.mean()), mape_in_domain=mape(ind),
                         mape_out_of_domain=mape(~ind), spearman_ad_ratio_vs_ape=float(rho),
                         frac_outside_range=float(np.mean([bool(o) for o in oor]))))
        for cid, idx in s.groupby('config_id').indices.items():
            m = np.zeros(len(s), bool)
            m[idx] = True
            per.append(dict(set=name, config_id=cid, n_samples=len(idx), frac_in_domain=float(ind[idx].mean()),
                            mpsp_mape=mape(m), median_ad_ratio=float(np.median(ratio[idx]))))
    return pd.DataFrame(rows), pd.DataFrame(per)


def raw_sample_coverage(files, samples, planned):
    """Per configuration: planned samples (``planned``: {biorefinery: n per
    configuration}, from ``data_generation/scripts/run_all.sh``), samples
    started / completed in the raw files (``sample_end`` not marked retry),
    completed without any simulated point because the system could not be
    restored after a failure (``empty_reset_failed``), and retained in the
    processed tables."""
    started, ended, empty = {}, {}, {}
    for f in files:
        with open(f) as fh:
            for line in fh:
                if '"sample_start"' not in line and '"sample_end"' not in line:
                    continue
                try:
                    rec = json.loads(line)
                except Exception:
                    continue
                sid = str(rec.get('sample_id', ''))
                parts = sid.split('|')
                if len(parts) < 2:
                    continue
                key = (parts[0], parts[1])
                if rec.get('record') == 'sample_start':
                    started.setdefault(key, set()).add(sid)
                elif rec.get('record') == 'sample_end' and not rec.get('retry'):
                    ended.setdefault(key, set()).add(sid)
                    if not rec.get('n_ok') and rec.get('reset_failed'):
                        empty.setdefault(key, set()).add(sid)
                    else:
                        empty.setdefault(key, set()).discard(sid)
    kept = samples.groupby(['biorefinery', 'config_id']).size()
    rows = []
    for key in sorted(set(started) | set(ended) | set(kept.index)):
        bio, cid = key
        rows.append(dict(biorefinery=bio, config_id=cid, planned=planned.get(bio), started=len(started.get(key, ())),
                         completed=len(ended.get(key, ())), empty_reset_failed=len(empty.get(key, ())),
                         retained=int(kept.get(key, 0))))
    return pd.DataFrame(rows)


def price_response_error(P, samples, keys=(('fs', 'feedstock_price'), ('el', 'electricity_price'),
                                           ('ng', 'natural_gas_price'))):
    """Accuracy of the predicted price sensitivity: per sample, the error of
    the predicted MPSP response to a 100 % change of each price,
    (dMPSP/dp_pred - dMPSP/dp_sim) * p / MPSP_sim, at the configuration's
    baseline yield and titer (dMPSP/dp = (i_a + i_b/y + i_c/t + i_d/(y t)) / g
    from the predicted intensities ``int__<key>__<k>`` and g, resp. the
    decomposed simulation data and the fitted g); for the feedstock also the
    relative error of dMPSP/dp. Returns the per-configuration medians."""
    s = samples.reset_index(drop=True)
    out = P.predict(s, warn=False)
    y, t = s['baseline_yield'].to_numpy(float), s['baseline_titer'].to_numpy(float)

    def dm(C, g):
        return (C[:, 0] + C[:, 1] / y + C[:, 2] / t + C[:, 3] / (y * t)) / g

    K = ['a', 'b', 'c', 'd']
    msim = dm(s[K].to_numpy(float), s['g'].to_numpy(float))
    res = pd.DataFrame(dict(config_id=s['config_id'].astype(str)))
    for key, price in keys:
        cs, cp = [f'int__{key}__{k}' for k in K], [f'int__{key}__{k}' for k in K]
        if not all(c in s for c in cs) or not all(c in out for c in cp):
            continue
        dsim = dm(s[cs].to_numpy(float), s['g'].to_numpy(float))
        dpred = dm(out[cp].to_numpy(float), out['g'].to_numpy(float))
        res[f'{key}_response_err'] = (dpred - dsim) * s[price].to_numpy(float) / msim
        if key == 'fs':
            with np.errstate(divide='ignore', invalid='ignore'):
                res['fs_sensitivity_rel_err'] = dpred / dsim - 1
    return res.groupby('config_id').median().reset_index()
