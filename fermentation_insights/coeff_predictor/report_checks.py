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
  configuration;
* :func:`baseline_family_stats` -- the fold-level comparison with the
  baselines at the level of the biorefinery families;
* :func:`loco_subset`, :func:`loco_per_family` -- LOCO over a subset of
  folds (e.g. without the new-product folds) and per family;
* :func:`input_twins`, :func:`final_test_without_configs` -- configurations
  with identical model inputs, and the final test without them;
* :func:`data_violation_rates` -- the constraint-violation rates of the
  data themselves (fitted coefficients), next to those of the predictions.
"""
import json
import os
import math

import numpy as np
import pandas as pd

__all__ = ('ablation_stats', 'paired_test_comparison', 'shipped_quantile_coverage', 'ad_validation',
           'raw_sample_coverage', 'sample_ape', 'baseline_fold_stats', 'baseline_family_stats',
           'single_config_products', 'loco_subset', 'loco_per_family', 'input_twins', 'drop_configs',
           'final_test_without_configs', 'baseline_bias', 'yt_definition_check', 'data_violation_rates', 'data_violation_per_sample',
           'data_violations_by_config')


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


def _fold_mape(res, scheme='config'):
    r = res[(res['scheme'] == scheme) & (res['fold'] >= 0) & (res['metric'] == 'mpsp_mape')
            & (res['target'] == 'MPSP')]
    return r.pivot_table(index='label', columns='method', values='value', aggfunc='mean')


def _sign_p(n_better, n_worse):
    if not n_better + n_worse:
        return math.nan
    try:
        from scipy.stats import binomtest
        return float(binomtest(n_worse, n_better + n_worse, 0.5).pvalue)
    except Exception:  # pragma: no cover
        return math.nan


def baseline_family_stats(res, families, scheme='config', methods=('knn', 'plain', 'ridge', 'hgb'),
                          reference='pinn'):
    """Family-level version of :func:`baseline_fold_stats` for a scheme
    whose folds cluster in families (LOCO: 37 configurations of 7
    biorefinery families, 15 of them HP). The fold-level tests treat the
    configurations as independent; here the paired fold differences
    (baseline - PI-NN) are first averaged per family (``families``:
    {fold label: family}), and the test runs over the families: mean, s.e.,
    families where the PI-NN is better, Wilcoxon signed-rank p (>= 5
    families) and two-sided sign-test p."""
    g = _fold_mape(res, scheme)
    rows = []
    if reference not in g:
        return pd.DataFrame(rows)
    for m in methods:
        if m not in g:
            continue
        d = (g[m] - g[reference]).dropna()
        fam = d.groupby(d.index.map(lambda x: families.get(x, x))).mean()
        st = _paired(fam)
        nb, nw = st['n_variant_better'], st['n_variant_worse']
        rows.append(dict(scheme=scheme, method=m, n_folds=int(len(d)), n_families=st['n_folds'],
                         mean=st['mean'], se=st['se'], median=st['median'], n_pinn_better=nw,
                         n_baseline_better=nb, p_wilcoxon=st['p_wilcoxon'], p_sign=_sign_p(nb, nw),
                         dominant_family=st.get('dominant_fold'), dominant_delta=st.get('dominant_delta')))
    return pd.DataFrame(rows)


def single_config_products(samples, col='fermentation_product'):
    """Configurations of ``samples`` whose fermentation product has no other
    configuration there: their leave-one-configuration-out folds hold out
    the product itself (new-product folds)."""
    c = samples.groupby('config_id')[col].first().astype(str)
    return sorted(c.index[c.map(c.value_counts()) == 1])


def loco_subset(res, exclude, methods=('pinn', 'plain', 'ridge', 'hgb', 'knn'), reference='pinn'):
    """Leave-one-configuration-out fold-mean MPSP MAPE of every method over
    the folds NOT in ``exclude``, with the paired fold statistics of each
    baseline against the PI-NN on those folds (:func:`baseline_fold_stats`)."""
    g = _fold_mape(res, 'config')
    g = g.loc[~g.index.isin(list(exclude))]
    sub = res[(res['scheme'] != 'config') | ~res['label'].isin(list(exclude))]
    bp = baseline_fold_stats(sub, schemes=('config',), methods=[m for m in methods if m != reference],
                             reference=reference).set_index('method')
    rows = []
    for m in methods:
        if m not in g:
            continue
        v = g[m].dropna()
        row = dict(method=m, n_folds=int(len(v)), mean=float(v.mean()),
                   se=float(v.std(ddof=1) / math.sqrt(len(v))) if len(v) > 1 else math.nan)
        if m in bp.index:
            x = bp.loc[m]
            row.update(diff_mean=x['mean'], diff_se=x['se'], n_pinn_better=int(x['n_pinn_better']),
                       p_wilcoxon=x['p_wilcoxon'])
        rows.append(row)
    return pd.DataFrame(rows)


def loco_per_family(res, families, method='pinn'):
    """Per family: the PI-NN's leave-one-configuration-out fold MAPE (mean,
    min, max and the configurations at the extremes)."""
    g = _fold_mape(res, 'config')
    if method not in g:
        return pd.DataFrame()
    v = g[method].dropna()
    fam = v.index.map(lambda x: families.get(x, x))
    rows = []
    for f, x in v.groupby(fam):
        rows.append(dict(family=str(f), n_folds=int(len(x)), mean=float(x.mean()), min=float(x.min()),
                         max=float(x.max()), min_config=str(x.idxmin()), max_config=str(x.idxmax())))
    return pd.DataFrame(rows).sort_values('mean').reset_index(drop=True)


def input_twins(P, samples, sets=None, tol=1e-9):
    """Configurations whose baseline rows encode to the SAME model input
    under checkpoint ``P`` (``spec.transform`` as in training, plus the
    productivity features and the standardized prices: ``spec.full_matrix``).
    Such configurations differ only in settings that are not inputs (e.g.
    the saccharification mode SHF / SSCF of the lactic configurations), so a
    held-out one is not new in input space. ``sets``: {config_id: set name}
    (e.g. 'test_ood' / 'development'). Returns (pairs table, dict of the
    closest pair that is not a twin: its configurations and Euclidean
    distance in the standardized features)."""
    s = samples.copy()
    s['_base'] = ((s['sample_kind'].astype(str) != 'baseline').astype(int) if 'sample_kind' in s else 0)
    keys = ['config_id', '_base'] + [c for c in ('sample_index', 'sample_id') if c in s]
    base = (s.sort_values(keys, kind='mergesort')
            .groupby('config_id', sort=True).head(1).drop(columns='_base').reset_index(drop=True))
    Z = np.asarray(P.spec.full_matrix(P.spec.transform(base, strict=False)), float)
    Z = np.where(np.isfinite(Z), Z, 0.)
    cids = base['config_id'].astype(str).tolist()
    sets = sets or {}
    rows, best = [], dict(config_id=None, other=None, distance=math.inf)
    for i in range(len(cids)):
        for j in range(i + 1, len(cids)):
            diff = np.abs(Z[i] - Z[j])
            dist = float(np.sqrt((diff ** 2).sum()))
            if diff.max() <= tol:
                a_, b_ = cids[i], cids[j]
                rows.append(dict(config_id=a_, twin=b_, set=sets.get(a_, ''), twin_set=sets.get(b_, ''),
                                 n_features=int(Z.shape[1]), n_features_different=int((diff > tol).sum()),
                                 crosses_split=sets.get(a_, '') != sets.get(b_, '')))
            elif dist < best['distance']:
                best = dict(config_id=cids[i], other=cids[j], distance=dist)
    return pd.DataFrame(rows, columns=['config_id', 'twin', 'set', 'twin_set', 'n_features',
                                       'n_features_different', 'crosses_split']), best


def yt_definition_check(samples, points, tol=0.01):
    """Measured vs. nominal yield and titer per configuration (ADAPTERS.md
    checklist item 2): per simulated point the ratios
    F_ferm_product / F_sugar_fed / y and F_ferm_product / V_broth / t
    (fermentation product leaving fermentation per sugar fed, resp. per broth
    volume, over the nominal y resp. t); min / median / max per
    configuration and whether either departs from 1 by more than ``tol``."""
    need = ('F_ferm_product', 'F_sugar_fed', 'V_broth', 'y', 't')
    if not all(c in points for c in need):
        return pd.DataFrame()
    cid = dict(zip(samples['sample_id'].astype(str), samples['config_id'].astype(str)))
    p = points[['sample_id', *need]].copy()
    p['config_id'] = p['sample_id'].astype(str).map(cid)
    with np.errstate(divide='ignore', invalid='ignore'):
        p['y_ratio'] = p['F_ferm_product'] / p['F_sugar_fed'] / p['y']
        p['t_ratio'] = p['F_ferm_product'] / p['V_broth'] / p['t']
    rows = []
    for c, g in p.dropna(subset=['config_id']).groupby('config_id'):
        r = dict(config_id=c, n_points=int(len(g)))
        for k in ('y_ratio', 't_ratio'):
            v = g[k].replace([np.inf, -np.inf], np.nan).dropna()
            r.update({f'{k}_min': float(v.min()), f'{k}_median': float(v.median()), f'{k}_max': float(v.max())})
        r['departs'] = bool(max(abs(r['y_ratio_min'] - 1), abs(r['y_ratio_max'] - 1), abs(r['t_ratio_min'] - 1),
                                abs(r['t_ratio_max'] - 1)) > tol)
        rows.append(r)
    return pd.DataFrame(rows)


def baseline_bias(P, samples, configs):
    """In-sample signed MPSP error of checkpoint ``P`` per configuration of
    ``configs`` (e.g. input twins, which get the same coefficients): the
    predicted MPSP at each sample's baseline yield and titer relative to the
    sample's fitted surface (the training target) - 1, for the baseline
    sample and as the mean over the configuration's samples."""
    from fermentation_insights.coeff_predictor.model import mpsp_equation_np
    s = samples[samples['config_id'].astype(str).isin([str(c) for c in configs])].reset_index(drop=True)
    if not len(s) or not {'baseline_yield', 'baseline_titer'} <= set(s.columns):
        return pd.DataFrame()
    out = P.predict(s, warn=False)
    Y = s['baseline_yield'].to_numpy(float)[:, None]
    T = s['baseline_titer'].to_numpy(float)[:, None]
    K = ['a', 'b', 'c', 'd']
    pred = mpsp_equation_np(out[K].to_numpy(float)[:, None, :], out['g'].to_numpy(float)[:, None], Y, T)[:, 0]
    fit = mpsp_equation_np(s[K].to_numpy(float)[:, None, :], s['g'].to_numpy(float)[:, None], Y, T)[:, 0]
    e = pd.Series(pred / fit - 1, index=s.index)
    rows = []
    for cid, idx in s.groupby('config_id').groups.items():
        b = s.loc[idx]
        bi = b.index[b['sample_kind'].astype(str) == 'baseline']
        rows.append(dict(config_id=str(cid), n_samples=int(len(idx)),
                         err_baseline_sample=float(e[bi[0]]) if len(bi) else math.nan,
                         err_mean=float(e[idx].mean()), coef_a=float(out.loc[bi[0], 'a']) if len(bi) else math.nan))
    return pd.DataFrame(rows)


def drop_configs(ev, samples_by_id, drop):
    """An evaluation dict (``ids``, ``coef``, ``g``, ``coef_std``, ``g_std``,
    ``scores``) without the samples of the configurations ``drop``; other
    (aggregate) entries are left out."""
    drop = {str(c) for c in drop}
    cid = samples_by_id.loc[list(ev['ids']), 'config_id'].astype(str).to_numpy()
    keep = ~np.isin(cid, list(drop))
    out = dict(ids=[i for i, k in zip(ev['ids'], keep) if k])
    for k in ('coef', 'g', 'coef_std', 'g_std'):
        v = ev.get(k)
        if v is not None and np.ndim(v) and len(v) == len(keep):
            out[k] = np.asarray(v, float)[keep]
    if ev.get('scores') is not None:
        sc = pd.DataFrame(ev['scores'])
        out['scores'] = sc[~sc['config_id'].astype(str).isin(drop)].reset_index(drop=True).to_dict('list')
    return out


def final_test_without_configs(final, samples_by_id, points_by_id, drop, set_name='test_ood', seed=0,
                         keys=('mpsp', 'g')):
    """Final-test numbers of one test set without the configurations
    ``drop`` (e.g. a held-out configuration whose model input is identical
    to a development configuration's): the paired comparison of every
    method with the PI-NN (:func:`paired_test_comparison`, configuration-
    cluster bootstrap) and the coverage of each network's intervals (its own
    ``conformal`` quantiles; pooled over samples / mean over configurations).
    Report only. Returns (paired table, coverage table)."""
    sub = {}
    for m, r in (final or {}).items():
        if r is None or r.get('error') or set_name not in r.get('evals', {}):
            continue
        sub[m] = dict(evals={set_name: drop_configs(r['evals'][set_name], samples_by_id, drop)},
                      conformal=r.get('conformal'))
    tp = paired_test_comparison(sub, samples_by_id, points_by_id, sets=(set_name,), seed=seed)
    if len(tp):
        tp['n_samples'] = [len(sub[m]['evals'][set_name]['ids']) for m in tp['method']]
    cov = []
    for m, r in sub.items():
        conf, ev = r.get('conformal') or {}, r['evals'][set_name]
        if not conf or ev.get('scores') is None:
            continue
        sc = pd.DataFrame(ev['scores'])
        for k in keys:
            q = conf.get(k)
            if q is None or k not in sc:
                continue
            v = sc[k].to_numpy(float)
            ok = np.isfinite(v)
            c = pd.Series(v[ok] <= q, dtype=float)
            byc = c.groupby(sc.loc[ok, 'config_id'].astype(str).to_numpy()).mean()
            cov.append(dict(method=m, set=set_name, key=k, n=int(ok.sum()), n_configs=int(byc.size),
                            coverage=float(c.mean()), coverage_cfg_mean=float(byc.mean())))
    return tp, pd.DataFrame(cov)


def data_violation_per_sample(samples_by_id, points_by_id, ids, tol=1e-3):
    """Per sample of ``ids``: whether the DATA (its fitted coefficients)
    violate ``mono_y``, ``mono_t``, ``positive`` (1 / 0 / NaN; definitions
    and tolerance of ``train.constraint_case_values``, the tolerance scaled
    by the spread of the coefficients over ``ids``) and ``g_le_gmax``, and
    g / g_max."""
    from fermentation_insights.coeff_predictor.train import _np_targets, constraint_case_values, G_MAX
    ids = [i for i in ids if i in points_by_id]
    if not ids:
        return pd.DataFrame()
    s = samples_by_id.loc[ids].reset_index(drop=True)
    p = pd.concat([points_by_id[i] for i in ids], ignore_index=True)
    coef, g, Y, T, mask = _np_targets(s, p)
    out = pd.DataFrame(dict(sample_id=s['sample_id'].astype(str), config_id=s['config_id'].astype(str)))
    for k, v in constraint_case_values(coef, g, Y, T, mask, None, None, None, tol).items():
        out[f'viol_{k}'] = 1 - v
    if G_MAX in s:
        gm = pd.to_numeric(s[G_MAX], errors='coerce').to_numpy(float)
        with np.errstate(divide='ignore', invalid='ignore'):
            out['g_over_gmax'] = g / gm
        out['viol_g_le_gmax'] = np.where(np.isfinite(gm) & np.isfinite(g), (g > gm * (1 + 1e-6)).astype(float),
                                         np.nan)
    return out


def data_violations_by_config(samples_by_id, points_by_id, ids, tol=1e-3):
    """Configurations whose DATA violate a constraint (:func:`data_violation_per_sample`
    over ``ids``): config_id, constraint, violating and checkable samples,
    and the maximum g / g_max."""
    d = data_violation_per_sample(samples_by_id, points_by_id, ids, tol)
    rows = []
    for c in [c for c in d.columns if c.startswith('viol_')]:
        v = d[np.isfinite(d[c].astype(float))]
        for cid, g in v.groupby('config_id'):
            n = int(g[c].sum())
            if n:
                rows.append(dict(constraint=c[5:], config_id=cid, n_violating=n, n_samples=int(len(g)),
                                 n_samples_all_configs=int(len(v)),
                                 max_g_over_gmax=float(g['g_over_gmax'].max()) if 'g_over_gmax' in g else math.nan))
    return pd.DataFrame(rows, columns=['constraint', 'config_id', 'n_violating', 'n_samples',
                                       'n_samples_all_configs', 'max_g_over_gmax'])


def data_violation_rates(samples_by_id, points_by_id, ids, tol=1e-3):
    """Constraint-violation rates of the DATA (the per-sample fitted
    coefficients, the targets) over the samples ``ids``, with the same
    definitions and tolerance as the rates of the predictions
    (``evaluate``: ``mono_y``, ``mono_t``, ``positive`` at the sample's
    points, ``g_le_gmax``). A model that reproduces the data exactly has
    these rates; e.g. isobutanol's MPSP rises with the yield (its yield
    deficit is sold as ethanol), so its samples violate ``mono_y`` in the
    data."""
    d = data_violation_per_sample(samples_by_id, points_by_id, ids, tol)
    return {c: float(np.nanmean(d[c].astype(float))) if np.isfinite(d[c].astype(float)).any() else math.nan
            for c in d.columns if c.startswith('viol_')}


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


#: markers of the generated accuracy block in README.md (``readme_summary``)
README_BEGIN = '<!-- BEGIN generated by report_checks.readme_summary (results/tables); do not edit by hand -->'
README_END = '<!-- END generated -->'
_README_SCHEMES = (('random', 'random (by sample)', 'interpolation: a known configuration, new parameter draws'),
                   ('config', 'leave-one-configuration-out', 'a new configuration of known feedstocks and products'),
                   ('feedstock', 'leave-one-feedstock-out', 'a new feedstock (4 folds)'),
                   ('product', 'leave-one-product-out', 'a new product family (7 folds)'),
                   ('productivity', 'productivity extrapolation', 'both tails of the productivity range'))
_README_METHODS = dict(pinn='PI-NN', plain='plain MLP', ridge='ridge', hgb='gradient boosting', knn='kNN')


def readme_summary(tables_dir):
    """Markdown block of the headline accuracy and coverage numbers of a
    pipeline run (from its ``results/tables``), kept in README.md between
    :data:`README_BEGIN` and :data:`README_END`; ``tests/test_readme_numbers.py``
    checks that the README carries the block of the committed tables."""
    t = lambda f: pd.read_csv(os.path.join(tables_dir, f))
    head, bp = t('cv_headline.csv'), t('cv_baselines_paired.csv')
    cf, shc = t('cv_coverage_crossfit.csv'), t('interval_coverage_by_scheme.csv')
    fh = t('final_test_headline.csv')
    pct = lambda v: '–' if not np.isfinite(v) else f'{100 * v:.1f}'
    L = ['| CV scheme (development set) | estimates | PI-NN MPSP MAPE % (± s.e.) | best baseline % | '
         'PI-NN vs. kNN, paired over folds: kNN − PI-NN pp (PI-NN better / folds) | coverage MPSP / g '
         '(nominal 0.90) |', '|---|---|---|---|---|---|']
    for sch, lab, what in _README_SCHEMES:
        h = head[head['scheme'] == sch].set_index('method')
        if 'pinn' not in h.index:
            continue
        base = h.drop(index=[m for m in h.index if m not in _README_METHODS or m == 'pinn'])['mpsp_mape']
        bm = base.idxmin()
        se = h.loc['pinn', 'mpsp_mape_se']
        k = bp[(bp['scheme'] == sch) & (bp['method'] == 'knn')]
        kn = (f'{100 * k["mean"].iloc[0]:+.1f} ± {100 * k["se"].iloc[0]:.1f} '
              f'({int(k["n_pinn_better"].iloc[0])}/{int(k["n_folds"].iloc[0])})' if len(k) else '–')
        if sch == 'config':
            c = cf[(cf['scheme'] == 'config') & (cf['method'] == 'pinn') & (cf['weighting'] == 'sample')]
            c = c.set_index('key')['pooled']
            cov = f'{c["mpsp"]:.2f} / {c["g"]:.2f} (cross-fitted)'
        elif sch in ('feedstock', 'product'):
            c = shc[shc['scheme'] == sch].set_index('key')['pooled']
            cov = f'{c["mpsp"]:.2f} / {c["g"]:.2f} (not calibrated)'
        else:
            cov = '–'
        if sch == 'config' and os.path.exists(os.path.join(tables_dir, 'cv_loco_known_product.csv')):
            k = t('cv_loco_known_product.csv')
            nk = int(k.loc[k['method'] == 'pinn', 'n_folds'].iloc[0])
            nx = len(str(k['excluded_folds'].iloc[0]).split(', '))
            what += f' (all {nk + nx} folds; {nx} of them hold out the only configuration of a product, see below)'
        L.append(f'| {lab} | {what} | {pct(h.loc["pinn", "mpsp_mape"])}'
                 + (f' (± {pct(se)})' if np.isfinite(se) else '') + f' | {_README_METHODS[bm]} '
                 f'{pct(base[bm])} | {kn} | {cov} |')
    L.append('')
    for st, lab in (('test_ood', 'held-out configurations (HP_corn, TAL_sugarcane, succinic_glucose, '
                                 'lactic_neutral_SSCF_cornstover)'), ('test_id', 'in-distribution samples')):
        f = fh[(fh['set'] == st) & fh['method'].isin(list(_README_METHODS))].set_index('method')
        if 'pinn' not in f.index:
            continue
        b = f.drop(index='pinn')['mpsp_mape']
        tw = ''
        if st == 'test_ood' and os.path.exists(os.path.join(tables_dir, 'final_test_without_twins.csv')):
            w = t('final_test_without_twins.csv')
            w = w[w['method'].isin(list(_README_METHODS))].set_index('method')
            pairs = t('input_twins.csv') if os.path.exists(os.path.join(tables_dir, 'input_twins.csv')) else None
            if 'pinn' in w.index and pairs is not None:
                x = pairs[pairs['crosses_split'].astype(bool)]
                held = sorted({r.config_id if r.set == 'test_ood' else r.twin for r in x.itertuples()})
                mate = sorted({r.twin if r.set == 'test_ood' else r.config_id for r in x.itertuples()})
                wb = w.drop(index='pinn')['mpsp_mape']
                tw = (f'; without {", ".join(held)}, whose model input is identical to that of the development '
                      f'configuration {", ".join(mate)} (REPORT.md section 2): PI-NN {pct(w.loc["pinn", "mpsp_mape"])} % '
                      f'(n = {int(w.loc["pinn", "n_samples"])}), best baseline {_README_METHODS[wb.idxmin()]} '
                      f'{pct(wb.min())} %')
        L.append(f'* Final test, {lab}: PI-NN {pct(f.loc["pinn", "mpsp_mape"])} % (n = '
                 f'{int(f.loc["pinn", "n_test"])} samples), best baseline {_README_METHODS[b.idxmin()]} '
                 f'{pct(b.min())} %' + tw + ' (`models/final_dev`; a re-used hold-out, REPORT.md).')
    kp = os.path.join(tables_dir, 'cv_loco_known_product.csv')
    if os.path.exists(kp):
        k = t('cv_loco_known_product.csv').set_index('method')
        pf = t('cv_loco_per_family.csv') if os.path.exists(os.path.join(tables_dir, 'cv_loco_per_family.csv')) \
            else None
        ex = str(k['excluded_folds'].iloc[0]).split(', ')
        rng = lambda x: f'{x.family} {pct(x.mean)}' + (f' ({pct(x.min)}-{pct(x.max)})' if x.n_folds > 1 else '')
        fam = ''
        if pf is not None:
            known = [rng(x) for x in pf.itertuples() if x.min_config not in ex]
            new = [rng(x) for x in pf.itertuples() if x.min_config in ex]
            fam = ('; PI-NN per family (fold mean, range): ' + ', '.join(known) + ' %'
                   + (f'; the new-product folds {", ".join(new)} %' if new else ''))
        L.append(f'* Leave-one-configuration-out over the {int(k.loc["pinn", "n_folds"])} folds of a KNOWN product '
                 f'(the folds {", ".join(ex)} hold out the only configuration of their product, so they are '
                 f'new-product folds): PI-NN {pct(k.loc["pinn", "mean"])} %, '
                 + ', '.join(f'{_README_METHODS[m]} {pct(k.loc[m, "mean"])} %' for m in ('plain', 'knn')
                             if m in k.index) + fam + '.')
    return '\n'.join([README_BEGIN] + L + [README_END])


def update_readme(readme, tables_dir):
    """Replace the generated block of ``readme`` by :func:`readme_summary`
    of ``tables_dir``; returns True if the file changed."""
    with open(readme, encoding='utf-8') as fh:
        txt = fh.read()
    i, j = txt.index(README_BEGIN), txt.index(README_END) + len(README_END)
    new = txt[:i] + readme_summary(tables_dir) + txt[j:]
    if new != txt:
        with open(readme, 'w', encoding='utf-8') as fh:
            fh.write(new)
    return new != txt


if __name__ == '__main__':      # python -m ...report_checks: refresh the README block
    _here = os.path.dirname(os.path.abspath(__file__))
    print('README.md updated:', update_readme(os.path.join(_here, 'README.md'),
                                               os.path.join(_here, 'results', 'tables')))
