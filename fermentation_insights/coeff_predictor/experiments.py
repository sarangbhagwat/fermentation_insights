#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Validation experiments of the coefficient predictor: ONE reproducible script
that produces every number, table and figure of the validation report.

    PYTHONPATH=/home/user/fermentation_insights python -m \\
        fermentation_insights.coeff_predictor.experiments [--quick] [--threads 3]

Stages (all deterministic; every model fit is cached, so a crashed or
interrupted run resumes where it stopped)
------------------------------------------------------------------------
1. ``data``    processed tables from the raw harness files (``dataset.build``),
               cached in ``data/processed/`` and rebuilt only when the raw files
               change (size/mtime fingerprint);
2. ``split``   final test set, FROZEN in ``results/split.json`` (committed) before
               any modelling and never used for selection: one whole
               configuration from every biorefinery with >= ``min_configs``
               PLANNED configurations (``PLANNED_CONFIGS``; seeded hash order,
               preferring feedstocks not yet held out) = out-of-distribution
               test, plus every sample of the other configurations whose seeded
               sample-id hash is < 10 % (in-distribution test). The choice depends
               only on the planned list, the seed and the sample ids, never on
               which data exist; the file is loaded, never regenerated, and the
               final stage refuses to run while a test configuration has no
               data. The rest = development set (``development_set``);
3. ``hpo``     random search (candidate 0 = package defaults) of width, depth,
               dropout, weight decay, learning rate, w_mpsp, w_ch, w_phys,
               collocation ratio, early-stopping criterion and input clipping by
               K-fold leave-configurations-out CV on the development set (folds
               stratified by biorefinery; ``hpo_members`` networks per fit),
               selecting by the fold mean of the configuration-balanced held-out
               MPSP-surface MAPE (``HPO_METRIC``). Full mode sizes the search
               (and, if needed, the number of leave-one-configuration-out folds)
               to ``--budget-hours`` from the measured cost of the first
               candidate (``plan.json``);
4. ``cv``      generalization CV on the development set with the selected
               hyperparameters: random-by-sample, leave-one-configuration-out,
               leave-one-feedstock-out, leave-one-product-out and productivity
               extrapolation; methods PI-NN ensemble, plain MLP, ridge, HGB, kNN
               and the PI-NN ablations no_channels / no_price / no_phys (config
               and product schemes);
5. ``final``   selected PI-NN ensemble (and the baselines) trained on the whole
               development set, evaluated ONCE on the final test set
               (``models/final_dev``); then the release ensemble trained on ALL
               data (``models/release``, ``train.Ensemble`` checkpoint format);
6. ``checks``  publication baseline coefficients vs. ours, predicted vs.
               simulated price sensitivities, channel breakdown of the test
               configurations;
7. ``figures`` / ``report``  PNGs (``figures/``), CSV tables (``tables/``) and
               ``REPORT.md``.

``--quick`` runs the same pipeline with a tiny search, 2 members and a few
epochs (outputs go to ``results/quick`` and ``models/quick``) for debugging.
``--threads N`` runs N single-threaded worker processes (one torch/BLAS thread
per fit was much faster than multi-threaded fits for these small networks).
Caches are keyed by the data hash, the split, the job specification and
``CODE_VERSION``; ``--fresh`` clears them.
"""
import argparse
import glob
import hashlib
import json
import os
import pickle
import platform
import shutil
import sys
import time
import traceback
import warnings

PKG_DIR = os.path.dirname(os.path.abspath(__file__))
MODULE = 'fermentation_insights.coeff_predictor.experiments'
CODE_VERSION = '3'

DEFAULT_RAW = '/home/user/work/data/raw'
DEFAULT_PROCESSED = os.path.join(PKG_DIR, 'data', 'processed')
DEFAULT_RESULTS = os.path.join(PKG_DIR, 'results')
DEFAULT_MODELS = os.path.join(PKG_DIR, 'models')
PUB_DIR = os.path.join(os.path.dirname(PKG_DIR), 'results', 'data')

TARGETS = ('a', 'b', 'c', 'd', 'g')
SCHEMES = ('random', 'config', 'feedstock', 'product', 'productivity')
#: inner calibration / early-stopping grouping per scheme. Like the package
#: default (train.SCHEME_GROUP) except leave-one-feedstock-out: with only four
#: feedstocks, grouping the calibration AND validation hold-outs by feedstock
#: would leave a single feedstock to train on, so they hold out whole
#: configurations instead (intervals are then calibrated for new
#: configurations, not for new feedstocks).
INNER_GROUP = {'random': None, 'config': 'config_id', 'feedstock': 'config_id',
               'product': 'fermentation_product', 'productivity': 'config_id'}
BASE_METHODS = ('pinn', 'plain', 'ridge', 'hgb', 'knn')
ABLATIONS = ('pinn-no_channels', 'pinn-no_price', 'pinn-no_phys')
ABLATION_SCHEMES = ('config', 'product')
NN_BASES = ('pinn', 'plain')

#: PI-NN extrapolation safeguards of the pipeline (need the downstream
#: FeatureSpec; DESIGN.md sections 3 and 6, round 4): the sugar-availability
#: prior of the structured feedstock intensity, the neutralization gate, and
#: clipping of the feedstock-composition inputs to the training range
EXTRAPOLATION_SAFEGUARDS = dict(fs_prior=True, neutralization_gate=True, clip_groups=['feedstock', 'pretreatment'])

HP_DEFAULT = dict(width=128, depth=3, dropout=0.05, weight_decay=1e-4, lr=2e-3,
                  w_mpsp=1.0, w_ch=0.5, w_phys=0.1, colloc_ratio=1.0, early_stop='data', clip_inputs=-1.0)
#: HPO selection criterion: fold mean of the configuration-balanced held-out
#: MPSP MAPE (every held-out configuration counts once). The point-weighted
#: MAPE of a fold is dominated by its largest configuration (e.g. oxalic: 189
#: of 911 development samples), so selecting on it selects on one
#: extrapolation.
HPO_METRIC = 'mpsp_mape_cfg'

#: colour-blind-safe categorical palette (fixed order; colour follows entity)
PALETTE = ('#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948')
INK, INK2, INK3, GRID = '#0b0b0b', '#52514e', '#8a8984', '#e4e3df'
METHOD_ORDER = BASE_METHODS + ABLATIONS
METHOD_COLOR = dict(zip(METHOD_ORDER, PALETTE))
METHOD_LABEL = {'pinn': 'PI-NN ensemble', 'plain': 'Plain MLP', 'ridge': 'Ridge', 'hgb': 'Gradient boosting',
                'knn': 'kNN (k = 5)', 'pinn-no_channels': 'PI-NN w/o channels',
                'pinn-no_price': 'PI-NN w/o price-explicit', 'pinn-no_phys': 'PI-NN w/o physics/colloc.'}
SCHEME_LABEL = {'random': 'Random\n(by sample)', 'config': 'Leave-one-\nconfiguration-out',
                'feedstock': 'Leave-one-\nfeedstock-out', 'product': 'Leave-one-\nproduct-out',
                'productivity': 'Productivity\nextrapolation'}
BIO_ORDER = ('HP', 'TAL', 'succinic', 'lactic', 'ethanol', 'oxalic', 'isobutanol')
BIO_MARKER = dict(zip(BIO_ORDER, ('o', 's', '^', 'D', 'v', 'P', 'X')))
BIO_COLOR = dict(zip(BIO_ORDER, PALETTE))
TARGET_LABEL = {'a': "a′ [USD/kg]", 'b': "b′ [USD/kg·(g/g)]", 'c': "c′ [USD/kg·(g/L)]",
                'd': "d′ [USD/kg·(g/g)·(g/L)]", 'g': 'g [kg/kg]'}
PRICE_INFO = {'fs': ('feedstock_price', 'feedstock', 'USD/kg wet feedstock'),
              'el': ('electricity_price', 'electricity', 'USD/kWh'),
              'ng': ('natural_gas_price', 'natural gas', 'USD/kg')}

SENS_STEP = 0.5   # relative price step of the sensitivity secants
_W = {}       # worker state: samples / points tables
_LOG = {}     # log file handle


# %% Utilities

def log(msg):
    line = f'[{time.strftime("%H:%M:%S")}] {msg}'
    print(line, flush=True)
    fh = _LOG.get('fh')
    if fh is not None:
        fh.write(line + '\n')
        fh.flush()


def _sha(obj, n=16):
    s = obj if isinstance(obj, str) else json.dumps(obj, sort_keys=True, default=str)
    return hashlib.sha1(s.encode()).hexdigest()[:n]


def _h(*parts):
    """Deterministic hash rank (independent of PYTHONHASHSEED)."""
    return int(hashlib.sha1('|'.join(map(str, parts)).encode()).hexdigest()[:12], 16)


def _dump_json(obj, path):
    tmp = path + '.tmp'
    with open(tmp, 'w') as fh:
        json.dump(obj, fh, indent=1, default=_json_default)
    os.replace(tmp, path)


def _json_default(o):
    import numpy as np
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (set, tuple)):
        return list(o)
    return str(o)


def _set_thread_env():
    for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
              'VECLIB_MAXIMUM_THREADS', 'NUMBA_NUM_THREADS'):
        os.environ[k] = '1'
    os.environ.setdefault('MPLBACKEND', 'Agg')


# %% Stage 1: data

def raw_fingerprint(raw_dirs):
    files = []
    for d in raw_dirs:
        for p in sorted(glob.glob(os.path.join(d, '**', '*.jsonl'), recursive=True)):
            st = os.stat(p)
            files.append([os.path.abspath(p), int(st.st_size), int(st.st_mtime)])
    return _sha(files), files


def data_hash(samples, points):
    cols = [c for c in ('sample_id', *TARGETS, 'fit_mape', 'r') if c in samples]
    return _sha(samples[cols].to_csv(index=False) + f'|{len(points)}|{points["MPSP"].sum():.10g}')


def stage_data(a):
    from fermentation_insights.coeff_predictor import dataset
    from fermentation_insights.coeff_predictor.train import load_tables
    proc = a.processed
    fp_path = os.path.join(proc, 'experiments_fingerprint.json')
    fp, files = raw_fingerprint(a.raw)
    have = all(os.path.exists(os.path.join(proc, f'{n}.csv.gz')) for n in ('samples', 'points', 'configs'))
    old = {}
    if os.path.exists(fp_path):
        with open(fp_path) as fh:
            old = json.load(fh)
    if have and (a.no_rebuild or (old.get('raw_fingerprint') == fp and not a.rebuild)):
        log(f'data: cached processed tables in {proc} '
            f'({"--no-rebuild" if a.no_rebuild else "raw files unchanged"})')
    else:
        if not files:
            raise SystemExit(f'no raw *.jsonl files under {a.raw}')
        t0 = time.time()
        log(f'data: building processed tables from {len(files)} raw files ...')
        dataset.build(a.raw, out_dir=proc)
        _dump_json(dict(raw_fingerprint=fp, raw_dirs=a.raw, files=files,
                        built=time.strftime('%Y-%m-%d %H:%M:%S')), fp_path)
        log(f'data: built in {time.time() - t0:.0f}s -> {proc}')
    # processed tables built before the downstream specification (patch S1):
    # add it on disk (idempotent), so that worker processes read the same tables
    from fermentation_insights.coeff_predictor import downstream
    PD = dataset.ProcessedData.load(proc)
    s_new, n_s = downstream.apply_to_frame(PD.samples)
    c_new, n_c = downstream.apply_to_frame(PD.configs)
    if n_s or n_c:
        dataset.ProcessedData(samples=s_new, points=PD.points, configs=c_new,
                              report=dict(PD.report, downstream_applied=n_s)).save(proc)
        log(f'data: added the downstream specification (S1, version {downstream.SPEC_VERSION}) to {n_s} '
            f'processed samples')
    samples, points, configs = load_tables(proc)
    samples = samples.sort_values('sample_id', kind='mergesort').reset_index(drop=True)
    rep = {}
    if os.path.exists(os.path.join(proc, 'report.json')):
        with open(os.path.join(proc, 'report.json')) as fh:
            rep = json.load(fh)
    now = time.time()
    recent = [f[0] for f in files if now - f[2] < 900]
    if recent:
        log(f'data: WARNING {len(recent)} raw file(s) modified in the last 15 min '
            '(campaigns may still be running; results describe a partial data set)')
    dh = data_hash(samples, points)
    log(f'data: {len(samples)} samples, {len(points)} points, {samples["config_id"].nunique()} configurations, '
        f'{samples["biorefinery"].nunique()} biorefineries (hash {dh})')
    return dict(samples=samples, points=points, configs=configs, report=rep, hash=dh,
                recent_raw_files=recent, n_raw_files=len(files))


def data_summary(samples, points):
    import numpy as np
    npts = points.groupby('sample_id').size()
    s = samples[['biorefinery', 'config_id', 'sample_id', 'sample_kind', 'feedstock', 'fermentation_product',
                 'final_product', 'R2', 'fit_mape', *TARGETS]].copy()
    s['n_pts'] = s['sample_id'].map(npts).fillna(0)
    g = s.groupby(['biorefinery', 'config_id'])
    out = g.agg(feedstock=('feedstock', 'first'), fermentation_product=('fermentation_product', 'first'),
                final_product=('final_product', 'first'), n_samples=('sample_id', 'size'),
                n_points=('n_pts', 'sum'), median_R2=('R2', 'median'), median_fit_mape=('fit_mape', 'median'),
                median_a=('a', 'median'), median_b=('b', 'median'), median_c=('c', 'median'),
                median_d=('d', 'median'), median_g=('g', 'median')).reset_index()
    kinds = s.groupby(['config_id', 'sample_kind']).size().unstack(fill_value=0)
    for k in ('baseline', 'productivity', 'random'):
        out[f'n_{k}'] = out['config_id'].map(kinds[k] if k in kinds else {}).fillna(0).astype(int)
    out['n_points'] = out['n_points'].astype(int)
    return out.sort_values(['biorefinery', 'config_id'], key=lambda c: c.map(_bio_rank) if c.name == 'biorefinery'
                           else c).reset_index(drop=True)


def data_flags(samples, max_fit_mape=0.05):
    """Configurations whose data look physically odd: median fit MAPE above
    ``max_fit_mape`` (4-coefficient form fits poorly) or a negative median
    feedstock-intensity b' (feedstock use per kg product should grow as 1/y;
    the PI-NN structures b'_fs = exp(z)/w_sugar > 0)."""
    import numpy as np
    import pandas as pd
    rows = []
    for cid, g in samples.groupby('config_id'):
        fm = float(np.nanmedian(pd.to_numeric(g['fit_mape'], errors='coerce')))
        if fm > max_fit_mape:
            rows.append(dict(config_id=cid, flag=f'median fit MAPE {100 * fm:.1f} %'))
        if 'int__fs__b' in g:
            b = float(np.nanmedian(pd.to_numeric(g['int__fs__b'], errors='coerce')))
            if b < 0:
                d = float(np.nanmedian(pd.to_numeric(g.get('int__fs__d'), errors='coerce')))
                rows.append(dict(config_id=cid, flag=f"negative feedstock-intensity b' (median {b:.3g}; d' {d:.3g})"))
    return pd.DataFrame(rows, columns=['config_id', 'flag'])


def _bio_rank(b):
    return BIO_ORDER.index(b) if b in BIO_ORDER else len(BIO_ORDER)


# %% Stage 2: final test split (frozen)

def _planned(*rows):
    return tuple(dict(config_id=c, biorefinery=b, feedstock=f) for c, b, f in rows)


_FS4 = ('glucose', 'corn', 'sugarcane', 'cornstover')
#: every configuration the data-generation campaigns plan to produce (the
#: adapters' ``configurations()``). The held-out test configurations are drawn
#: from THIS list and the split seed only -- never from whichever data happen
#: to exist when the split is made -- so the split cannot drift as campaigns
#: finish.
PLANNED_CONFIGS = _planned(
    *[(f'{v}_{f}', 'HP', f) for v in ('HP', 'HP_neutral', 'HP_hexanol', 'HP_neutral_hexanol') for f in _FS4],
    *[(f'{v}_{f}', 'TAL', f) for v in ('TAL', 'TAL_SA') for f in _FS4],
    *[(f'{v}_{f}', 'succinic', f) for v in ('succinic', 'succinic_neutral') for f in _FS4],
    *[(f'lactic_{r}_{k}_cornstover', 'lactic', 'cornstover') for r in ('neutral', 'lowpH') for k in ('SSCF', 'SHF')],
    *[(f'ethanol_{f}', 'ethanol', f) for f in ('cornstover', 'sugarcane', 'corn')],
    ('oxalic_neutral_sugarcane', 'oxalic', 'sugarcane'),
    ('IBO_EtOH_corn', 'isobutanol', 'corn'))
SPLIT_FILE = os.path.join(DEFAULT_RESULTS, 'split.json')
SPLIT_RULE_VERSION = 2


def _unit_hash(*parts):
    """Deterministic uniform number in [0, 1) from ``parts``."""
    return _h(*parts) / float(16 ** 12)


def make_frozen_split(seed, planned=PLANNED_CONFIGS, min_configs=4, id_frac=0.10):
    """The final test split as a pure function of the PLANNED configuration
    list and the seed.

    * out-of-distribution test: one whole configuration from every
      biorefinery with >= ``min_configs`` planned configurations
      (most-constrained biorefinery first; within it the configurations in
      seeded hash order, preferring a feedstock not held out yet);
    * in-distribution test: every sample of the remaining configurations whose
      seeded hash of its ``sample_id`` falls below ``id_frac`` (whole
      samples, i.e. grouped by sample). Membership depends on the sample id
      only, so it is stable as more samples arrive (a sample never moves
      between development and test).
    """
    cfg = {c['config_id']: c for c in planned}
    if len(cfg) != len(planned):
        raise ValueError('duplicate config_id in the planned configuration list')
    by_bio = {}
    for c in planned:
        by_bio.setdefault(c['biorefinery'], []).append(c['config_id'])
    eligible = {b: sorted(c) for b, c in by_bio.items() if len(c) >= max(2, min_configs)}
    used_fs, test_cfgs = set(), []
    for b in sorted(eligible, key=lambda b: (len({cfg[c]['feedstock'] for c in eligible[b]}), len(eligible[b]), b)):
        cands = sorted(eligible[b], key=lambda c: _h(seed, 'test-config', c))
        pick = next((c for c in cands if cfg[c]['feedstock'] not in used_fs), cands[0])
        used_fs.add(cfg[pick]['feedstock'])
        test_cfgs.append(pick)
    test_cfgs = sorted(test_cfgs, key=lambda c: (_bio_rank(cfg[c]['biorefinery']), c))
    return dict(
        rule_version=SPLIT_RULE_VERSION, seed=int(seed), min_configs=int(min_configs), id_frac=float(id_frac),
        rule=(f'OOD: one whole configuration from every biorefinery with >= {max(2, min_configs)} PLANNED '
              f'configurations (most-constrained biorefinery first; sha1 rank of "seed|test-config|config_id", '
              f'preferring a feedstock not held out yet). ID: every sample of the other configurations with '
              f'sha1("seed|id-test|sample_id")[:12] / 16^12 < {id_frac:g}. Depends only on the planned list, the '
              f'seed and the sample id; never on which data exist.'),
        planned_configs=[dict(c) for c in planned],
        test_configs=[dict(cfg[c]) for c in test_cfgs],
        biorefineries_without_held_out_config=sorted(set(by_bio) - set(eligible)))


def is_id_test(sample_ids, seed, id_frac):
    import numpy as np
    return np.array([_unit_hash(seed, 'id-test', str(s)) < id_frac for s in sample_ids], bool)


def split_hash(frozen):
    keys = ('rule_version', 'seed', 'min_configs', 'id_frac', 'planned_configs', 'test_configs')
    return _sha({k: frozen[k] for k in keys})


def apply_split(samples, frozen, require_test_configs=False):
    """(split info with the current test sample ids, development table).
    ``require_test_configs``: raise if a held-out configuration has no data
    (final evaluation)."""
    test_cfgs = [c['config_id'] for c in frozen['test_configs']]
    planned = {c['config_id'] for c in frozen['planned_configs']}
    cid = samples['config_id'].astype(str)
    present = set(cid)
    missing = [c for c in test_cfgs if c not in present]
    if missing and require_test_configs:
        raise SystemExit(f'final evaluation: held-out test configuration(s) {missing} have no samples in the data; '
                         f'the frozen split ({SPLIT_FILE}) must not be changed -- generate their data first')
    unplanned = sorted(present - planned)
    sid = samples['sample_id'].astype(str)
    is_ood = cid.isin(test_cfgs).to_numpy()
    is_id = ~is_ood & is_id_test(sid, frozen['seed'], frozen['id_frac'])
    dev = samples[~is_ood & ~is_id].reset_index(drop=True)
    sp = dict(frozen, id_test_sample_ids=sorted(sid[is_id]), ood_test_sample_ids=sorted(sid[is_ood]),
              n_samples=int(len(samples)), n_dev=int(len(dev)), n_test_ood=int(is_ood.sum()),
              n_test_id=int(is_id.sum()), n_dev_configs=int(dev['config_id'].nunique()),
              test_configs_missing=missing, unplanned_configs=unplanned,
              test_configs=[dict(c, n_samples=int((cid == c['config_id']).sum())) for c in frozen['test_configs']])
    sp['hash'] = split_hash(frozen)
    return sp, dev


def load_frozen_split(path=SPLIT_FILE):
    with open(path) as fh:
        fr = json.load(fh)
    if fr.get('rule_version') != SPLIT_RULE_VERSION:
        raise SystemExit(f'{path}: split rule version {fr.get("rule_version")} != {SPLIT_RULE_VERSION}')
    regen = make_frozen_split(fr['seed'], _planned(*[(c['config_id'], c['biorefinery'], c['feedstock'])
                                                     for c in fr['planned_configs']]),
                              fr['min_configs'], fr['id_frac'])
    if regen['test_configs'] != fr['test_configs']:
        raise SystemExit(f'{path}: test configurations do not follow from its own rule (edited by hand?)')
    return fr


def development_set(samples, path=SPLIT_FILE):
    """Development table only (for diagnostics / model development)."""
    return apply_split(samples, load_frozen_split(path))[1]


def stage_split(a, D):
    path = a.split_file
    if os.path.exists(path):
        frozen = load_frozen_split(path)
        asked = (a.split_seed, a.min_test_configs, a.id_test_frac)
        have = (frozen['seed'], frozen['min_configs'], frozen['id_frac'])
        if a.split_seed_given and asked != have:
            raise SystemExit(f'split: {path} is frozen with (seed, min_configs, id_frac) = {have}; '
                             f'refusing to use {asked}')
        if [dict(c) for c in PLANNED_CONFIGS] != frozen['planned_configs']:
            log(f'split: NOTE the planned configuration list in the code differs from {path}; the frozen '
                'split is used unchanged')
    else:
        frozen = make_frozen_split(a.split_seed, PLANNED_CONFIGS, a.min_test_configs, a.id_test_frac)
        frozen['created'] = time.strftime('%Y-%m-%d %H:%M:%S')
        _dump_json(frozen, path)
        log(f'split: froze a new split -> {path}')
    sp, dev = apply_split(D['samples'], frozen)
    if sp['test_configs_missing']:
        log(f"split: WARNING held-out configuration(s) without data yet: {sp['test_configs_missing']} "
            '(the final evaluation refuses to run until they have data)')
    if sp['unplanned_configs']:
        log(f"split: WARNING configurations not in the planned list (kept in development): {sp['unplanned_configs']}")
    _dump_json({k: v for k, v in sp.items() if k not in ('planned_configs',)}, os.path.join(a.out, 'split_applied.json'))
    log(f"split: test configurations {[c['config_id'] for c in sp['test_configs']]} ({sp['n_test_ood']} samples) + "
        f"{sp['n_test_id']} in-distribution test samples; development set {len(dev)} samples / "
        f"{dev['config_id'].nunique()} configurations (frozen split {path})")
    return sp, dev


# %% Folds

def stratified_group_kfold(df, k, seed, group='config_id', strata='biorefinery'):
    """K folds of whole groups, groups of each stratum dealt round-robin in
    hash order (continuing across strata, so fold sizes balance)."""
    import numpy as np
    g = df[group].astype(str).to_numpy()
    st = df.groupby(group)[strata].first()
    fold_of, i = {}, 0
    for s in sorted(set(st), key=lambda b: (_bio_rank(b), b)):
        for gg in sorted(st.index[st == s], key=lambda x: _h(seed, 'kfold', x)):
            fold_of[str(gg)] = i % k
            i += 1
    f = np.array([fold_of[x] for x in g])
    return [(np.flatnonzero(f != j), np.flatnonzero(f == j), f'fold{j}') for j in range(k) if (f == j).any()]


def stratified_subset(df, n, seed, group='config_id', strata='biorefinery'):
    """``n`` groups, allocated to strata proportionally (largest remainder,
    at least one per stratum), hash order within a stratum."""
    import numpy as np
    st = df.groupby(group)[strata].first()
    by = {s: sorted(map(str, st.index[st == s]), key=lambda x: _h(seed, 'subset', x)) for s in set(st)}
    if n >= len(st):
        return sorted(map(str, st.index))
    keys = sorted(by, key=lambda b: (_bio_rank(b), b))
    alloc = {b: 1 for b in keys} if n >= len(keys) else {b: 0 for b in keys}
    if n < len(keys):
        for b in sorted(keys, key=lambda b: -len(by[b]))[:n]:
            alloc[b] = 1
    left = n - sum(alloc.values())
    if left > 0:
        cap = {b: len(by[b]) - alloc[b] for b in keys}
        tot = sum(cap.values())
        quota = {b: left * cap[b] / tot for b in keys}
        add = {b: int(np.floor(quota[b])) for b in keys}
        rem = left - sum(add.values())
        for b in sorted(keys, key=lambda b: -(quota[b] - add[b]))[:rem]:
            add[b] += 1
        for b in keys:
            alloc[b] += min(add[b], cap[b])
    return sorted(x for b in keys for x in by[b][:alloc[b]])


def scheme_splits(dev, scheme, plan, seed):
    """[(fold_index, train_idx, test_idx, label)] on the development table."""
    from fermentation_insights.coeff_predictor.train import split_random, split_group, split_productivity
    if scheme == 'random':
        sp = split_random(dev, plan['random_k'], seed)
    elif scheme == 'config':
        sp = split_group(dev, 'config_id')
    elif scheme == 'feedstock':
        sp = split_group(dev, 'feedstock')
    elif scheme == 'product':
        sp = split_group(dev, 'fermentation_product')
    elif scheme == 'productivity':
        sp = split_productivity(dev)
    else:
        raise ValueError(scheme)
    out = [(f, tr, te, lab) for f, (tr, te, lab) in enumerate(sp) if len(tr) and len(te)]
    if scheme == 'config' and plan.get('loco_folds') and plan['loco_folds'] < len(out):
        keep = set(stratified_subset(dev, plan['loco_folds'], seed))
        out = [x for x in out if x[3] in keep]
    return out


# %% Jobs (fits), worker processes and the cache

def _worker_init(processed_dir, expected_hash):
    _set_thread_env()
    warnings.filterwarnings('ignore')
    import torch
    torch.set_num_threads(1)
    from fermentation_insights.coeff_predictor.train import load_tables
    s, p, _ = load_tables(processed_dir)
    s = s.sort_values('sample_id', kind='mergesort').reset_index(drop=True)
    h = data_hash(s, p)
    if expected_hash and h != expected_hash:
        raise RuntimeError(f'worker data hash {h} != {expected_hash} (processed tables changed?)')
    s['sample_id'] = s['sample_id'].astype(str)
    p = p.copy()
    p['sample_id'] = p['sample_id'].astype(str)
    _W['samples'] = s.set_index('sample_id', drop=False)
    _W['points'] = p
    _W['pts'] = {k: v for k, v in p.groupby('sample_id', sort=False)}


def _subset(ids):
    import pandas as pd
    S = _W['samples']
    s = S.loc[list(ids)].reset_index(drop=True)
    parts = [_W['pts'][i] for i in ids if i in _W['pts']]
    p = pd.concat(parts, ignore_index=True) if parts else _W['points'].iloc[:0]
    return s, p


def _price_sensitivities(model, df):
    """Central secant slopes (price +-50 %; exact for price-linear models,
    and wide enough to see the response of piecewise-constant trees and to
    keep float32 round-off small) of a', b', c', d', g w.r.t. each price
    column: {key: (n, 5)}."""
    import numpy as np
    import pandas as pd
    out = {}
    for key, (col, _, _) in PRICE_INFO.items():
        if col not in df:
            continue
        p = pd.to_numeric(df[col], errors='coerce').to_numpy(float)
        h = SENS_STEP * np.abs(p)
        ok = np.isfinite(p) & (h > 0)
        if not ok.any():
            continue
        lo, hi = df.copy(), df.copy()
        lo[col] = np.where(ok, p - h, p)
        hi[col] = np.where(ok, p + h, p)
        a, b = model.predict(lo), model.predict(hi)
        d = np.column_stack([b['coef'] - a['coef'], b['g'] - a['g']])
        out[key] = d / np.where(ok, 2 * h, np.nan)[:, None]
    return out


def _run_job(job):
    """Fit one method on the job's training samples and evaluate it on each
    of its evaluation sets. Never raises (errors are returned)."""
    import numpy as np
    from fermentation_insights.coeff_predictor.evaluate import fit_method, compute_metrics
    t0 = time.time()
    out = dict(key=job['key'], stage=job['stage'], tag=job['tag'], method=job['method'], error=None,
               evals={}, info={}, epochs=[], n_train=len(job['train_ids']))
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter('always')
        try:
            s_tr, p_tr = _subset(job['train_ids'])
            model = fit_method(job['method'], s_tr, p_tr, job['members'], job['mk'], job['lk'], job['tk'],
                               group_col=job['group_col'], calib_frac=job['calib_frac'], seed=job['seed'],
                               spec_kw=job.get('spec_kw') or None)
            out['fit_seconds'] = time.time() - t0
            out['info'] = {k: v for k, v in (getattr(model, 'fit_info', {}) or {}).items()}
            ens = getattr(model, 'ens', None)
            if ens is not None:
                out['epochs'] = [int(e) for e in ens.meta.get('epochs', [])]
                out['conformal'] = dict(ens.conformal)
                out['n_fit_rows'] = int(ens.meta.get('n_fit', 0))
            for name, ids in job['evals'].items():
                if not len(ids):
                    continue
                s_te, p_te = _subset(ids)
                pred = model.predict(s_te)
                ev = dict(ids=list(ids), rows=compute_metrics(pred, s_te, p_te, model),
                          coef=np.asarray(pred['coef'], float), g=np.asarray(pred['g'], float),
                          g_max=None if pred.get('g_max') is None else np.asarray(pred['g_max'], float))
                mem = pred.get('members')
                if mem is not None:
                    ev['coef_std'], ev['g_std'] = mem['coef'].std(0), mem['g'].std(0)
                    if 'channels' in job['extras']:
                        if 'channels' in mem:
                            ev['channels'] = mem['channels'].mean(0)
                            ev['channel_names'] = list(ens.channel_names())
                        if 'intensities' in mem:
                            ev['intensities'] = mem['intensities'].mean(0)
                            ev['price_keys'] = list(ens.layout['price_keys'])
                if 'sens' in job['extras']:
                    ev['sens'] = _price_sensitivities(model, s_te)
                out['evals'][name] = ev
            if job.get('save_model'):
                if ens is None:
                    raise ValueError('save_model needs an ensemble method')
                path = job['save_model']
                if os.path.isdir(path):
                    shutil.rmtree(path)
                ens.save(path)
                _dump_json(dict(job_key=job['key'], tag=job['tag'], method=job['method'], members=job['members'],
                                mk=job['mk'], lk=job['lk'], tk=job['tk'], seed=job['seed'],
                                group_col=job['group_col'], calib_frac=job['calib_frac'],
                                n_train_samples=len(job['train_ids']), data_hash=job.get('data_hash'),
                                saved=time.strftime('%Y-%m-%d %H:%M:%S')),
                           os.path.join(path, 'experiment_info.json'))
                out['saved_model'] = path
        except Exception as e:  # recorded, the pipeline continues
            out['error'] = f'{type(e).__name__}: {str(e)[:300]}'
            out['traceback'] = traceback.format_exc()[-3000:]
    out['warnings'] = sorted({f'{w.category.__name__}: {str(w.message)[:160]}' for w in wl})[:12]
    out['seconds'] = time.time() - t0
    return out


def make_job(stage, tag, method, train_ids, evals, members, mk, lk, tk, group_col, calib_frac, seed, ctx,
             extras=(), save_model=None, meta=None, spec_kw=None):
    """One fit + evaluation. ``spec_kw``: FeatureSpec options of every method
    (e.g. ``{'downstream': True}``; default ``ctx['spec_kw']``, off with CLI
    ``--no-downstream``; part of the cache key only when non-empty).
    ``ctx['mk_extra']`` (water gate, off with CLI ``--no-water-gate``; channel-wise
    input masking, off with ``--no-masking``; EXTRAPOLATION_SAFEGUARDS, off with
    ``--no-extrapolation-safeguards``) is merged into ``mk``. A masked
    PI-NN without explicit sub-network sizes gets ``masked_width = width // 2``
    and ``masked_depth = depth - 1`` (64 / 2 for the defaults, as in the
    development CV), so that the HPO dimensions width and depth still apply."""
    if spec_kw is None:
        spec_kw = ctx.get('spec_kw')
    if ctx.get('mk_extra'):
        mk = {**mk, **ctx['mk_extra']}
    if mk.get('channel_inputs') and method.split('-')[0] == 'pinn':
        mk = dict(mk)
        mk.setdefault('masked_width', max(8, int(mk.get('width', 128)) // 2))
        mk.setdefault('masked_depth', max(1, int(mk.get('depth', 3)) - 1))
    spec = dict(v=CODE_VERSION, data=ctx['data_hash'], split=ctx['split_hash'], stage=stage, tag=tag,
                method=method, members=members, mk=mk, lk=lk, tk=tk, group_col=group_col, calib_frac=calib_frac,
                seed=seed, train=_sha(list(train_ids)), evals={k: _sha(list(v)) for k, v in evals.items()},
                extras=sorted(extras), save=bool(save_model))
    if spec_kw:
        spec['spec_kw'] = dict(spec_kw)
    return dict(key=_sha(spec, 24), stage=stage, tag=tag, method=method, train_ids=list(train_ids),
                evals={k: list(v) for k, v in evals.items()}, members=members, mk=mk, lk=lk, tk=tk,
                spec_kw=dict(spec_kw or {}),
                group_col=group_col, calib_frac=calib_frac, seed=seed, extras=tuple(extras),
                save_model=save_model, meta=meta or {}, data_hash=ctx['data_hash'],
                cost=(members if method.split('-')[0] in NN_BASES else 0.3) * len(train_ids))


class Runner:
    """Runs jobs in-process (``workers == 1``) or in a pool of spawned
    single-threaded worker processes; caches each finished job."""

    def __init__(self, workers, processed_dir, dhash, cache_dir, force_pool=False):
        self.workers, self.processed, self.dhash, self.cache = int(workers), processed_dir, dhash, cache_dir
        self.force_pool = bool(force_pool)
        self.pool = None
        self.M = None

    def _mod(self):
        if self.M is None:
            import importlib
            self.M = importlib.import_module(MODULE)
        return self.M

    def path(self, job):
        return os.path.join(self.cache, job['stage'], job['key'] + '.pkl')

    def cached(self, job):
        p = self.path(job)
        if not os.path.exists(p):
            return None
        if job.get('save_model') and not os.path.exists(os.path.join(job['save_model'], 'meta.json')):
            return None
        try:
            with open(p, 'rb') as fh:
                return pickle.load(fh)
        except Exception:
            return None

    def _store(self, job, res):
        p = self.path(job)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        if res.get('error') is None:      # failures are retried on the next run
            tmp = p + '.tmp'
            with open(tmp, 'wb') as fh:
                pickle.dump(res, fh, protocol=4)
            os.replace(tmp, p)

    def run(self, jobs, name):
        out, todo = {}, []
        for j in jobs:
            r = self.cached(j)
            if r is None:
                todo.append(j)
            else:
                out[j['key']] = r
        log(f'{name}: {len(jobs)} fits ({len(jobs) - len(todo)} cached, {len(todo)} to run, '
            f'{self.workers} worker(s))')
        if not todo:
            return out
        todo.sort(key=lambda j: -j['cost'])        # longest first
        total = sum(j['cost'] for j in todo)
        done_cost, t0 = 0., time.time()
        M = self._mod()

        def _done(j, r):
            nonlocal done_cost
            self._store(j, r)
            out[j['key']] = r
            done_cost += j['cost']
            el = time.time() - t0
            eta = el * (total - done_cost) / max(done_cost, 1e-9)
            ev = r['evals'].get('test') or next(iter(r['evals'].values()), None)
            mm = {(m, t): v for m, t, v in ev['rows']} if ev else {}
            msg = (f'mape={mm.get(("mpsp_mape", "MPSP"), float("nan")):.4f} fit_mape='
                   f'{mm.get(("fit_mape", "MPSP"), float("nan")):.4f}' if ev else '')
            if r['error']:
                msg = 'FAILED ' + r['error']
            log(f'  [{name} {j["tag"]}] {j["method"]:<17s} {msg} ({r["seconds"]:.0f}s, epochs {r["epochs"]}) '
                f'[{len([k for k in out if k in {x["key"] for x in todo}])}/{len(todo)}, '
                f'elapsed {el / 60:.1f} min, ETA {eta / 60:.1f} min]')

        if self.workers <= 1 and not self.force_pool:
            if not M._W:
                M._worker_init(self.processed, self.dhash)
            for j in todo:
                _done(j, M._run_job(j))
        else:
            from concurrent.futures import ProcessPoolExecutor, as_completed
            import multiprocessing as mp
            if self.pool is None:
                self.pool = ProcessPoolExecutor(max(1, self.workers), mp_context=mp.get_context('spawn'),
                                                initializer=M._worker_init, initargs=(self.processed, self.dhash))
            futs = {self.pool.submit(M._run_job, j): j for j in todo}
            for f in as_completed(futs):
                j = futs[f]
                try:
                    r = f.result()
                except Exception as e:  # worker crash
                    r = dict(key=j['key'], tag=j['tag'], method=j['method'], error=f'worker: {e!r}', evals={},
                             epochs=[], seconds=0., info={}, warnings=[])
                _done(j, r)
        return out

    def close(self):
        if self.pool is not None:
            self.pool.shutdown()
            self.pool = None


# %% Plan (budget)

def base_plan(a):
    if a.quick:
        return dict(mode='quick', hpo_k=3, hpo_candidates=3, hpo_members=1, cv_members=2, final_members=2,
                    release_members=2, random_k=3, loco_folds=4, epochs=12, patience=6, min_epochs=4, batch_size=128)
    # hpo_members: single-member fits differ between seeds by up to ~50 % (relative)
    # in held-out MAPE on the development folds, so one member per candidate
    # and fold mostly ranks seeds, not hyperparameters
    return dict(mode='full', hpo_k=5, hpo_candidates=None, hpo_members=a.hpo_members, cv_members=a.members,
                final_members=5, release_members=5, random_k=5, loco_folds=None, epochs=300, patience=40,
                min_epochs=20, batch_size=128)


def hpo_candidates(n, seed):
    """Fixed seeded sequence (candidate 0 = package defaults); the first n."""
    import numpy as np
    rng = np.random.RandomState(seed + 2027)
    out = [dict(HP_DEFAULT)]
    sig = lambda v: float(f'{v:.2g}')
    while len(out) < max(n, 1):
        c = dict(width=int(rng.choice([64, 128, 256])), depth=int(rng.choice([2, 3, 4])),
                 dropout=float(rng.choice([0., 0.05, 0.1, 0.2])),
                 weight_decay=sig(10 ** rng.uniform(-5, -3)), lr=sig(10 ** rng.uniform(np.log10(5e-4), np.log10(5e-3))),
                 w_mpsp=float(rng.choice([0.3, 1., 3.])), w_ch=float(rng.choice([0.1, 0.3, 0.5, 1.])),
                 w_phys=float(rng.choice([0.03, 0.1, 0.3, 1.])), colloc_ratio=float(rng.choice([0.5, 1., 2.])),
                 early_stop=str(rng.choice(['data', 'mpsp'])), clip_inputs=float(rng.choice([-1., 0.])))
        if c not in out:
            out.append(c)
    return out[:n]


def hp_kwargs(hp, plan):
    mk = dict(width=int(hp['width']), depth=int(hp['depth']), dropout=float(hp['dropout']),
              clip_inputs=float(hp.get('clip_inputs', -1.0)))
    lk = dict(w_mpsp=float(hp['w_mpsp']), w_ch=float(hp['w_ch']), w_phys=float(hp['w_phys']),
              colloc_ratio=float(hp['colloc_ratio']))
    tk = dict(lr=float(hp['lr']), weight_decay=float(hp['weight_decay']), epochs=plan['epochs'],
              patience=plan['patience'], min_epochs=plan['min_epochs'], batch_size=plan['batch_size'], threads=1,
              early_stop=str(hp.get('early_stop', 'data')))
    return mk, lk, tk


def plan_budget(plan, a, pilot, dev, test_n, n_all, seed):
    """Full mode: size the HPO search and the number of LOCO folds so that
    the estimated wall time fits ``--budget-hours`` with ``--threads``
    workers. Cost model: seconds per member fit proportional to the number
    of training samples (measured on the HPO pilot = candidate 0)."""
    import numpy as np
    per = pilot['seconds_per_member_sample']
    hgb = pilot.get('hgb_seconds_per_sample', 0.004)

    def fit(n, members, method='pinn'):
        return members * per * n if method.split('-')[0] in NN_BASES else hgb * n

    def cv_cost(loco_folds):
        tot = 0.
        for sch in SCHEMES:
            sp = scheme_splits(dev, sch, dict(plan, loco_folds=loco_folds), seed)
            meths = list(BASE_METHODS) + (list(ABLATIONS) if sch in ABLATION_SCHEMES else [])
            for _, tr, _, _ in sp:
                tot += sum(fit(len(tr), plan['cv_members'], m) for m in meths if m != 'ridge' and m != 'knn')
        return tot

    n_dev = len(dev)
    final = fit(n_dev, plan['final_members']) * 2 + fit(n_all, plan['release_members']) + 3 * fit(n_dev, 1, 'hgb')
    hpo_one = plan['hpo_k'] * fit(n_dev * (plan['hpo_k'] - 1) / plan['hpo_k'], plan.get('hpo_members', 1))
    budget = a.budget_hours * 3600 * max(1, a.threads)
    n_loco = dev['config_id'].nunique()
    n_cand = a.hpo_candidates or 12
    loco = a.loco_folds or None
    total = lambda nc, lf: nc * hpo_one + cv_cost(lf) + final
    if not a.hpo_candidates and not a.loco_folds:
        # shrink LOCO folds (min 12), then candidates (min 8); grow candidates (max 16) with slack
        while total(n_cand, loco) > budget and (loco or n_loco) > 12:
            loco = max(12, (loco or n_loco) - 2)
        while total(n_cand, loco) > budget and n_cand > 8:
            n_cand -= 1
        while n_cand < 16 and total(n_cand + 1, loco) <= budget:
            n_cand += 1
    est = total(n_cand, loco)
    plan.update(hpo_candidates=int(n_cand), loco_folds=None if (loco is None or loco >= n_loco) else int(loco),
                estimate=dict(seconds_per_member_sample=per, hpo_hours=n_cand * hpo_one / 3600 / max(1, a.threads),
                              cv_hours=cv_cost(loco) / 3600 / max(1, a.threads),
                              final_hours=final / 3600 / max(1, a.threads),
                              total_hours=est / 3600 / max(1, a.threads), budget_hours=a.budget_hours,
                              threads=a.threads, fits_budget=bool(est <= budget)))
    return plan


# %% Results assembly

def tidy_rows(res, scheme, fold, label, method, n_train, n_test):
    rows = []
    info = res.get('info') or {}
    if res.get('error'):
        mets = [('error', res['error'][:200], float('nan'))]
    else:
        mets = res['evals']['test']['rows']
    for metric, target, value in mets:
        rows.append(dict(scheme=scheme, fold=fold, label=label, method=method, n_train=n_train, n_test=n_test,
                         n_fit=info.get('n_fit', float('nan')), n_fit_configs=info.get('n_fit_configs', float('nan')),
                         n_calib=info.get('n_calib', float('nan')), calib_group=info.get('calib_group_used'),
                         seconds=res.get('seconds', float('nan')),
                         epochs=float(sum(res.get('epochs') or [])) or float('nan'),
                         metric=metric, target=target, value=value))
    return rows


def pooled_rows(scheme, parts, samples_by_id, points_by_id):
    """Pooled out-of-fold metrics over the common (all-methods-ok) folds."""
    import numpy as np
    import pandas as pd
    from fermentation_insights.coeff_predictor.evaluate import compute_metrics
    rows, oof = [], []
    for meth, evs in parts.items():
        if not evs:
            continue
        ids = [i for ev in evs for i in ev['ids']]
        pred = dict(coef=np.concatenate([ev['coef'] for ev in evs]), g=np.concatenate([ev['g'] for ev in evs]),
                    g_max=np.concatenate([ev['g_max'] if ev['g_max'] is not None else np.full(len(ev['g']), np.nan)
                                          for ev in evs]))
        s = samples_by_id.loc[ids].reset_index(drop=True)
        p = pd.concat([points_by_id[i] for i in ids if i in points_by_id], ignore_index=True)
        for metric, target, value in compute_metrics(pred, s, p):
            if metric in ('MAE', 'RMSE', 'R2', 'mpsp_mape', 'mpsp_mape_median', 'mpsp_mape_cfg',
                          'mpsp_mape_cfg_gmean', 'oracle_mape', 'fit_mape'):
                rows.append(dict(scheme=scheme, fold=-1, label='pooled', method=meth, n_train=np.nan, n_test=len(s),
                                 seconds=np.nan, metric=metric, target=target, value=value))
        fold = np.concatenate([np.full(len(ev['ids']), ev['_fold']) for ev in evs])
        oof.append(pd.DataFrame(dict(scheme=scheme, method=meth, fold=fold, sample_id=ids,
                                     **{f'{k}_pred': pred['coef'][:, j] for j, k in enumerate(TARGETS[:4])},
                                     g_pred=pred['g'])))
    return rows, oof


def headline(res, methods=None):
    """One row per (scheme, method): fold-mean MPSP MAPE (+- s.e.), pooled
    MAPE, fit MAPE, pooled R2 / MAE per target, violation and coverage."""
    import numpy as np
    import pandas as pd
    out = []
    err = res[res['metric'] == 'error']
    bad = set(zip(err['scheme'], err['fold']))
    for (sch, meth), g in res.groupby(['scheme', 'method'], sort=False):
        per = g[(g['fold'] >= 0) & ~np.array([(s, f) in bad for s, f in zip(g['scheme'], g['fold'])], bool)]
        pool = g[g['fold'] == -1]
        v = lambda d, m, t: d[(d['metric'] == m) & (d['target'] == t)]['value']
        mp = v(per, 'mpsp_mape', 'MPSP')
        row = dict(scheme=sch, method=meth, n_folds=int(g[g['fold'] >= 0]['fold'].nunique()),
                   n_folds_ok=int(g[(g['fold'] >= 0) & (g['metric'] != 'error')]['fold'].nunique()),
                   n_folds_common=int(per['fold'].nunique()),
                   mpsp_mape=mp.mean(), mpsp_mape_se=mp.std(ddof=1) / np.sqrt(len(mp)) if len(mp) > 1 else np.nan,
                   mpsp_mape_median=v(per, 'mpsp_mape_median', 'MPSP').mean(),
                   mpsp_mape_cfg=v(per, 'mpsp_mape_cfg', 'MPSP').mean(),
                   mpsp_mape_cfg_gmean=v(per, 'mpsp_mape_cfg_gmean', 'MPSP').mean(),
                   mpsp_mape_pooled=v(pool, 'mpsp_mape', 'MPSP').mean(),
                   fit_mape=v(per, 'fit_mape', 'MPSP').mean(), oracle_mape=v(per, 'oracle_mape', 'MPSP').mean(),
                   # spec-determined (bp_* at baseline) vs random samples; Gram-norm coefficient error
                   mpsp_mape_cfg_gmean_spec=v(per, 'mpsp_mape_cfg_gmean_spec', 'MPSP').mean(),
                   mpsp_mape_cfg_gmean_random=v(per, 'mpsp_mape_cfg_gmean_random', 'MPSP').mean(),
                   coef_gram_rel=v(per, 'coef_gram_rel', 'abcd').mean(),
                   coef_gram_rel_cfg=v(per, 'coef_gram_rel_cfg', 'abcd').mean())
        for k in TARGETS:
            src = pool if len(v(pool, 'R2', k)) else per
            row[f'R2_{k}'] = v(src, 'R2', k).mean()
            row[f'MAE_{k}'] = v(src, 'MAE', k).mean()
        for k in ('mono_y', 'mono_t', 'positive', 'mono_r', 'g_le_gmax', 'nonneg_fs', 'nonneg_ng'):
            row[f'viol_{k}'] = v(per, 'violation_rate', k).mean()
        for k in ('MPSP',) + TARGETS:
            row[f'coverage_{k}'] = v(per, 'coverage', k).mean()
        out.append(row)
    df = pd.DataFrame(out)
    if methods:
        df['_m'] = df['method'].map({m: i for i, m in enumerate(methods)})
        df['_s'] = df['scheme'].map({s: i for i, s in enumerate(SCHEMES + ('final_test',))})
        df = df.sort_values(['_s', '_m']).drop(columns=['_m', '_s'])
    return df.reset_index(drop=True)


# %% Stages 3-5

def stage_hpo(a, plan, runner, dev, ctx, n_test, n_all):
    folds = stratified_group_kfold(dev, plan['hpo_k'], a.seed)
    ids = dev['sample_id'].astype(str).to_numpy()
    n_max = 16 if plan['hpo_candidates'] is None else plan['hpo_candidates']
    cands = hpo_candidates(max(n_max, 1), a.seed)

    def jobs_for(ci):
        mk, lk, tk = hp_kwargs(cands[ci], plan)
        return [make_job('hpo', f'cand{ci} {lab}', 'pinn', ids[tr], {'test': ids[te]}, plan.get('hpo_members', 1),
                         mk, lk, tk, 'config_id',
                         0.0, a.seed + 101 * ci + f, ctx, meta=dict(cand=ci, fold=f))
                for f, (tr, te, lab) in enumerate(folds)]

    results = runner.run(jobs_for(0), 'hpo pilot (candidate 0 = defaults)')
    if plan['hpo_candidates'] is None:      # full mode: size the search from the pilot
        pl_path = os.path.join(a.out, 'plan.json')
        old = {}
        if os.path.exists(pl_path) and not a.replan:
            with open(pl_path) as fh:
                old = json.load(fh)
        if old.get('data_hash') == ctx['data_hash'] and old.get('split_hash') == ctx['split_hash'] \
                and old.get('mode') == 'full' and old.get('threads') == a.threads:
            plan.update({k: old[k] for k in ('hpo_candidates', 'loco_folds', 'estimate')})
            log(f'plan: reusing {pl_path}')
        else:
            ok = [r for r in results.values() if not r['error'] and r.get('n_fit_rows')]
            per = (sum(r['fit_seconds'] for r in ok) / sum(r['n_fit_rows'] * len(r['epochs'] or [1]) for r in ok)
                   if ok else 0.05)
            plan_budget(plan, a, dict(seconds_per_member_sample=per), dev, n_test, n_all, a.seed)
        e = plan['estimate']
        log(f"plan: {plan['hpo_candidates']} HPO candidates, LOCO folds = {plan['loco_folds'] or 'all'}; "
            f"estimated {e['total_hours']:.2f} h on {e['threads']} worker(s) (HPO {e['hpo_hours']:.2f}, "
            f"CV {e['cv_hours']:.2f}, final {e['final_hours']:.2f}; budget {e['budget_hours']} h"
            f"{'' if e['fits_budget'] else ' -- EXCEEDS BUDGET'})")
    if a.plan_only:
        return None
    cands = cands[:plan['hpo_candidates']]
    jobs = [j for ci in range(len(cands)) for j in jobs_for(ci)]
    results.update(runner.run([j for j in jobs if j['meta']['cand'] > 0], 'hpo search'))
    import numpy as np
    import pandas as pd
    rows, fold_rows = [], []
    for ci, hp in enumerate(cands):
        js = [j for j in jobs if j['meta']['cand'] == ci]
        rs = [results[j['key']] for j in js]
        mm = []
        for j, r in zip(js, rs):
            m = {} if r['error'] else {(x, t): v for x, t, v in r['evals']['test']['rows']}
            mm.append(m)
            fold_rows.append(dict(candidate=ci, fold=j['meta']['fold'], error=r['error'],
                                  mpsp_mape=m.get(('mpsp_mape', 'MPSP'), np.nan),
                                  mpsp_mape_cfg=m.get(('mpsp_mape_cfg', 'MPSP'), np.nan),
                                  mpsp_mape_cfg_gmean=m.get(('mpsp_mape_cfg_gmean', 'MPSP'), np.nan),
                                  fit_mape=m.get(('fit_mape', 'MPSP'), np.nan), seconds=r['seconds'],
                                  epochs=np.mean(r['epochs']) if r['epochs'] else np.nan,
                                  **{f'R2_{k}': m.get(('R2', k), np.nan) for k in TARGETS}))
        mape = np.array([m.get(('mpsp_mape', 'MPSP'), np.nan) for m in mm])
        crit = np.array([m.get((HPO_METRIC, 'MPSP'), np.nan) for m in mm])
        n_ok = int(np.isfinite(crit).sum())
        rows.append(dict(candidate=ci, **hp, n_folds_ok=n_ok,
                         criterion=float(np.nanmean(crit)) if n_ok == len(crit) else np.inf,
                         criterion_se=float(np.nanstd(crit, ddof=1) / np.sqrt(n_ok)) if n_ok > 1 else np.nan,
                         mpsp_mape=float(np.nanmean(mape)) if np.isfinite(mape).all() else np.inf,
                         mpsp_mape_se=float(np.nanstd(mape, ddof=1) / np.sqrt(len(mape))) if len(mape) > 1 else np.nan,
                         mpsp_mape_cfg_gmean=float(np.nanmean([m.get(('mpsp_mape_cfg_gmean', 'MPSP'), np.nan)
                                                               for m in mm])),
                         mpsp_mape_median=float(np.nanmean([m.get(('mpsp_mape_median', 'MPSP'), np.nan) for m in mm])),
                         fit_mape=float(np.nanmean([m.get(('fit_mape', 'MPSP'), np.nan) for m in mm])),
                         seconds_per_fit=float(np.mean([r['seconds'] for r in rs])),
                         epochs=float(np.nanmean([np.mean(r['epochs']) if r['epochs'] else np.nan for r in rs]))))
    tab = pd.DataFrame(rows)
    best = int(tab.sort_values(['criterion', 'criterion_se', 'candidate']).iloc[0]['candidate'])
    tab['selected'] = tab['candidate'] == best
    sel = dict(cands[best]) if np.isfinite(tab.loc[best, 'criterion']) else dict(HP_DEFAULT)
    if not np.isfinite(tab.loc[best, 'criterion']):
        log('hpo: WARNING every candidate failed in some fold; falling back to the package defaults')
    fold_info = [dict(fold=lab, n_test_configs=int(len(set(dev['config_id'].iloc[te]))), n_test=int(len(te)))
                 for tr, te, lab in folds]
    log(f'hpo: selected candidate {best}: {sel} (fold-mean configuration-balanced held-out MPSP MAPE '
        f'{tab.loc[best, "criterion"]:.4f}; point-weighted {tab.loc[best, "mpsp_mape"]:.4f})')
    return dict(table=tab, folds=pd.DataFrame(fold_rows), selected=sel, best=best, fold_info=fold_info)


def stage_cv(a, plan, runner, dev, ctx, hp, D):
    import pandas as pd
    mk, lk, tk = hp_kwargs(hp, plan)
    ids = dev['sample_id'].astype(str).to_numpy()
    jobs, meta = [], {}
    for sch in SCHEMES:
        meths = list(BASE_METHODS) + (list(ABLATIONS) if sch in ABLATION_SCHEMES else [])
        for f, tr, te, lab in scheme_splits(dev, sch, plan, a.seed):
            for m in meths:
                j = make_job('cv', f'{sch} {lab}', m, ids[tr], {'test': ids[te]}, plan['cv_members'], mk, lk, tk,
                             INNER_GROUP[sch], 0.15, a.seed + f, ctx, meta=dict(scheme=sch, fold=f, label=lab))
                jobs.append(j)
    results = runner.run(jobs, 'cv')
    S = D['samples'].assign(sample_id=D['samples']['sample_id'].astype(str)).set_index('sample_id', drop=False)
    P = {k: v for k, v in D['points'].assign(sample_id=D['points']['sample_id'].astype(str)).groupby('sample_id')}
    rows, oofs = [], []
    for sch in SCHEMES:
        js = [j for j in jobs if j['meta']['scheme'] == sch]
        failed = {j['meta']['fold'] for j in js if results[j['key']]['error']}
        parts = {}
        for j in js:
            r = results[j['key']]
            rows += tidy_rows(r, sch, j['meta']['fold'], j['meta']['label'], j['method'], len(j['train_ids']),
                              len(j['evals']['test']))
            if not r['error'] and j['meta']['fold'] not in failed:
                ev = dict(r['evals']['test'], _fold=j['meta']['fold'])
                parts.setdefault(j['method'], []).append(ev)
        if failed:
            log(f'cv {sch}: pooled metrics exclude fold(s) {sorted(failed)} (a method failed there)')
        if len({j['meta']['fold'] for j in js}) > 1:
            pr, oof = pooled_rows(sch, parts, S, P)
            rows += pr
            oofs += oof
    res = pd.DataFrame(rows)
    oof = pd.concat(oofs, ignore_index=True) if oofs else pd.DataFrame()
    return dict(results=res, oof=oof, jobs=jobs, raw=results)


def stage_final(a, plan, runner, dev, ctx, hp, sp, D, models_dir):
    if sp['test_configs_missing']:
        msg = (f"final: held-out test configuration(s) {sp['test_configs_missing']} have no samples; the frozen "
               f"split ({a.split_file}) is never changed -- generate their data before the final evaluation")
        if not a.quick:
            raise SystemExit(msg)
        log('WARNING (quick mode only) ' + msg)
    mk, lk, tk = hp_kwargs(hp, plan)
    s = D['samples']
    dev_ids = dev['sample_id'].astype(str).tolist()
    evals = dict(test_ood=list(sp['ood_test_sample_ids']), test_id=list(sp['id_test_sample_ids']),
                 test_all=sorted(set(sp['ood_test_sample_ids']) | set(sp['id_test_sample_ids'])))
    all_ids = s['sample_id'].astype(str).tolist()
    jobs = []
    for m in BASE_METHODS:
        nn = m in NN_BASES
        jobs.append(make_job('final', 'final(dev->test)', m, dev_ids, evals, plan['final_members'], mk, lk, tk,
                             'config_id', 0.15, a.seed, ctx,
                             extras=('channels', 'sens') if m == 'pinn' else ('sens',),
                             save_model=os.path.join(models_dir, 'final_dev') if m == 'pinn' else None,
                             meta=dict(kind='final', nn=nn)))
    # release: ALL samples are fitted (calib_frac = 0; the members' early-stopping
    # hold-outs rotate), intervals are transferred from the dev model (below)
    jobs.append(make_job('final', 'release(all data)', 'pinn', all_ids, dict(train_all=all_ids),
                         plan['release_members'], mk, lk, tk, 'config_id', 0.0, a.seed, ctx,
                         save_model=os.path.join(models_dir, 'release'), meta=dict(kind='release')))
    results = runner.run(jobs, 'final + release')
    import pandas as pd
    rows = []
    for j in jobs:
        r = results[j['key']]
        if r['error']:
            log(f"final: {j['tag']} {j['method']} FAILED: {r['error']}")
            rows.append(dict(set='error', method=j['method'], model=j['tag'], metric='error', target=r['error'],
                             value=float('nan')))
            continue
        for name, ev in r['evals'].items():
            for metric, target, value in ev['rows']:
                rows.append(dict(set=name, method=j['method'], model=j['tag'], n_train=len(j['train_ids']),
                                 n_test=len(ev['ids']), metric=metric, target=target, value=value))
    out = dict(jobs=jobs, raw=results, table=pd.DataFrame(rows))
    out['final'] = {j['method']: results[j['key']] for j in jobs if j['meta']['kind'] == 'final'}
    out['release'] = next(results[j['key']] for j in jobs if j['meta']['kind'] == 'release')
    rel = os.path.join(models_dir, 'release')
    out['release_conformal'] = None
    if not out['release']['error'] and not out['final']['pinn']['error']:
        try:
            out['release_conformal'] = transfer_conformal(os.path.join(models_dir, 'final_dev'), rel, dev_ids,
                                                          sp, plan, D, a.seed)
            c = out['release_conformal']
            log(f"release: intervals calibrated with the dev model on {c['n']} samples of {c['n_calib_groups']} "
                f"held-out configurations (level {c.get('mpsp_level', float('nan')):.3f}) -> {rel}")
        except Exception as e:
            log(f'release: WARNING conformal transfer failed: {type(e).__name__}: {e}')
    # release smoke test through the user-facing API
    out['release_check'] = None
    if not out['release']['error']:
        try:
            from fermentation_insights.coeff_predictor.predict import CoefficientPredictor
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                P = CoefficientPredictor.load(rel)
                base = s[s['sample_kind'].astype(str) == 'baseline'].head(3)
                pr = P.predict(base, warn=False)
            out['release_check'] = dict(ok=True, n=len(pr), columns=int(pr.shape[1]),
                                        example={k: float(pr[k].iloc[0]) for k in TARGETS},
                                        intervals=all(f'{k}_lo' in pr for k in TARGETS),
                                        example_config=str(base['config_id'].iloc[0]))
            log(f'release: {rel} loads with CoefficientPredictor; prediction for {base["config_id"].iloc[0]} '
                f'baseline: ' + ', '.join(f'{k}={pr[k].iloc[0]:.4g}' for k in TARGETS)
                + f'; intervals: {out["release_check"]["intervals"]}')
        except Exception as e:
            out['release_check'] = dict(ok=False, error=f'{type(e).__name__}: {e}')
            log(f'release: WARNING CoefficientPredictor check failed: {type(e).__name__}: {e}')
    return out


def transfer_conformal(dev_path, rel_path, dev_ids, sp, plan, D, seed):
    """Split-conformal quantiles for the release ensemble (which fits ALL
    data, so nothing is left to calibrate it on): scores of the dev-trained
    ensemble (same architecture and hyperparameters, slightly less data) on
    every configuration it never saw -- its own calibration configurations
    plus the held-out test configurations -- grouped by configuration. Valid
    for the dev ensemble; for the release ensemble an approximation (expected
    to be conservative), hence ``coverage_guaranteed = False``."""
    import pandas as pd
    from fermentation_insights.coeff_predictor.train import Ensemble, ensemble_splits, calibrate_conformal
    ens_dev = Ensemble.load(dev_path)
    S = D['samples'].assign(sample_id=D['samples']['sample_id'].astype(str)).set_index('sample_id', drop=False)
    dev_s = S.loc[dev_ids].reset_index(drop=True)
    spl = ensemble_splits(dev_s, plan['final_members'], 0.15, ens_dev.tcfg.val_frac, 'config_id', seed)
    if len(spl['calib']) != int(ens_dev.meta.get('n_calib', -1)):
        raise RuntimeError('cannot reproduce the calibration split of the dev ensemble')
    ids = list(dev_s['sample_id'].iloc[spl['calib']]) + list(sp['ood_test_sample_ids'])
    s_cal = S.loc[ids].reset_index(drop=True)
    pts = D['points'][D['points']['sample_id'].astype(str).isin(set(ids))]
    conf = dict(calibrate_conformal(ens_dev, s_cal, pts, alpha=0.1, group_col='config_id'))
    conf.update(transferred=True, coverage_guaranteed=False,
                calibration_source='dev-trained ensemble (models/final_dev) scored on its calibration '
                                   'configurations + the held-out test configurations',
                calibration_configs=sorted(set(s_cal['config_id'].astype(str))))
    ens_rel = Ensemble.load(rel_path)
    ens_rel.conformal = conf
    ens_rel.save(rel_path)
    return conf


def final_headline(table):
    import numpy as np
    import pandas as pd
    t = table[table['set'] != 'error']
    out = []
    for (st, meth), g in t.groupby(['set', 'method'], sort=False):
        v = lambda m, k: g[(g['metric'] == m) & (g['target'] == k)]['value'].mean()
        row = dict(set=st, method=meth, model=g['model'].iloc[0], n_train=int(g['n_train'].iloc[0]),
                   n_test=int(g['n_test'].iloc[0]), mpsp_mape=v('mpsp_mape', 'MPSP'),
                   mpsp_mape_median=v('mpsp_mape_median', 'MPSP'), fit_mape=v('fit_mape', 'MPSP'))
        for k in TARGETS:
            row[f'R2_{k}'], row[f'MAE_{k}'] = v('R2', k), v('MAE', k)
        for k in ('mono_y', 'mono_t', 'positive', 'mono_r', 'g_le_gmax', 'nonneg_fs', 'nonneg_ng'):
            row[f'viol_{k}'] = v('violation_rate', k)
        for k in ('MPSP',) + TARGETS:
            row[f'coverage_{k}'] = v('coverage', k)
        out.append(row)
    df = pd.DataFrame(out)
    if len(df):
        df['_s'] = df['set'].map({'test_ood': 0, 'test_id': 1, 'test_all': 2, 'train_all': 3})
        df['_m'] = df['method'].map({m: i for i, m in enumerate(METHOD_ORDER)})
        df = df.sort_values(['_s', '_m']).drop(columns=['_s', '_m']).reset_index(drop=True)
    return df


# %% Stage 6: checks

def pub_comparison(samples, test_cfgs):
    import numpy as np
    import pandas as pd
    rows = []
    base = samples[samples['sample_kind'].astype(str) == 'baseline'].sort_values('sample_id')
    for cid, g in base.groupby('config_id'):
        f = os.path.join(PUB_DIR, f'{cid}_coefficients.npy')
        if not os.path.exists(f):
            continue
        pub = np.load(f).astype(float)
        r = g.iloc[0]
        row = dict(config_id=cid, biorefinery=r['biorefinery'], feedstock=r['feedstock'],
                   in_final_test=cid in test_cfgs, sample_id=r['sample_id'])
        for j, k in enumerate(TARGETS + ('R2',)):
            ours = float(r[k])
            row[f'{k}_ours'], row[f'{k}_pub'] = ours, float(pub[j])
            row[f'{k}_ratio'] = ours / pub[j] if pub[j] != 0 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def sensitivity_table(final, samples_by_id):
    import numpy as np
    import pandas as pd
    rows, pts = [], []
    for meth, r in final.items():
        if r['error'] or 'test_all' not in r['evals']:
            continue
        ev = r['evals']['test_all']
        s = samples_by_id.loc[ev['ids']]
        for key, d in (ev.get('sens') or {}).items():
            for j, k in enumerate(TARGETS[:4]):
                col = f'int__{key}__{k}'
                if col not in s:
                    continue
                y = pd.to_numeric(s[col], errors='coerce').to_numpy(float)
                x = d[:, j]
                m = np.isfinite(x) & np.isfinite(y)
                if m.sum() < 2:
                    continue
                ss = np.sum((y[m] - y[m].mean()) ** 2)
                # exact linearity of price-explicit models: FD derivative == predicted intensity
                lin = np.nan
                pk = ev.get('price_keys') or []
                if ev.get('intensities') is not None and key in pk:
                    q = ev['intensities'][:, pk.index(key), j]
                    lin = float(np.nanmax(np.abs(x - q)) / max(np.nanmax(np.abs(q)), 1e-12))
                rows.append(dict(method=meth, price=key, coefficient=k, n=int(m.sum()), fd_vs_intensity_rel=lin,
                                 R2=float(1 - np.sum((x[m] - y[m]) ** 2) / ss) if ss > 0 else np.nan,
                                 MAE=float(np.mean(np.abs(x[m] - y[m]))), mean_simulated=float(np.mean(y[m])),
                                 mean_predicted=float(np.mean(x[m])),
                                 slope=float(np.sum(x[m] * y[m]) / np.sum(y[m] ** 2)) if np.sum(y[m] ** 2) > 0
                                 else np.nan))
                pts.append(pd.DataFrame(dict(method=meth, price=key, coefficient=k, sample_id=np.asarray(ev['ids'])[m],
                                             simulated=y[m], predicted=x[m])))
    return pd.DataFrame(rows), (pd.concat(pts, ignore_index=True) if pts else pd.DataFrame())


def channel_table(final_pinn, samples_by_id, test_cfgs):
    import numpy as np
    import pandas as pd
    if final_pinn is None or final_pinn['error'] or 'test_ood' not in final_pinn['evals']:
        return pd.DataFrame()
    ev = final_pinn['evals']['test_ood']
    if 'channels' not in ev:
        return pd.DataFrame()
    s = samples_by_id.loc[ev['ids']].reset_index(drop=True)
    names = ev['channel_names']
    rows = []
    for cid in test_cfgs:
        m = (s['config_id'] == cid).to_numpy()
        if not m.any():
            continue
        for ci, ch in enumerate(names):
            for j, k in enumerate(TARGETS[:4]):
                col = f'ch__{ch}__{k}'
                sim = pd.to_numeric(s.loc[m, col], errors='coerce').to_numpy(float) if col in s else np.full(m.sum(),
                                                                                                            np.nan)
                pred = ev['channels'][m, ci, j]
                rows.append(dict(config_id=cid, channel=ch, coefficient=k, n=int(m.sum()),
                                 simulated_mean=float(np.nanmean(sim)) if np.isfinite(sim).any() else np.nan,
                                 predicted_mean=float(np.mean(pred)),
                                 MAE=float(np.nanmean(np.abs(pred - sim))) if np.isfinite(sim).any() else np.nan))
    return pd.DataFrame(rows)


def per_group_mape(ev, samples_by_id, points_by_id, col):
    """MPSP MAPE / fit MAPE per group (e.g. config_id) of an evaluation set."""
    import numpy as np
    import pandas as pd
    from fermentation_insights.coeff_predictor.train import pad_points
    from fermentation_insights.coeff_predictor.model import mpsp_equation_np
    s = samples_by_id.loc[ev['ids']].reset_index(drop=True)
    p = pd.concat([points_by_id[i] for i in ev['ids'] if i in points_by_id], ignore_index=True)
    Y, T, M, _, mask = pad_points(s, p)
    Mh = mpsp_equation_np(ev['coef'][:, None, :], ev['g'][:, None], Y, T)
    with np.errstate(divide='ignore', invalid='ignore'):
        ape = np.where(mask, np.abs(Mh - M) / np.abs(M), np.nan)
    rows = []
    for gval, idx in s.groupby(col).indices.items():
        a = ape[idx][mask[idx]]
        rows.append(dict(**{col: gval}, n_samples=len(idx), mpsp_mape=float(np.nanmean(a)),
                         mpsp_mape_median=float(np.nanmedian(a)),
                         fit_mape=float(pd.to_numeric(s.loc[idx, 'fit_mape'], errors='coerce').mean())))
    return pd.DataFrame(rows)


# %% Stage 7: figures

def _style():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
        'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'legend.fontsize': 8, 'axes.edgecolor': INK3, 'axes.labelcolor': INK, 'xtick.color': INK2,
        'ytick.color': INK2, 'text.color': INK, 'axes.spines.top': False, 'axes.spines.right': False,
        'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.axisbelow': True,
        'legend.frameon': False, 'lines.linewidth': 2, 'savefig.dpi': 200, 'figure.dpi': 100})
    return plt


def _save(fig, path, made):
    fig.savefig(path, dpi=200, bbox_inches='tight')
    import matplotlib.pyplot as plt
    plt.close(fig)
    made.append(path)


def _parity_ax(ax, x, y, groups, colors, markers, k, legend_labels=None):
    import numpy as np
    m = np.isfinite(x) & np.isfinite(y)
    for gname in [gg for gg in colors if (groups == gg).any()]:
        mm = m & (groups == gname)
        ax.scatter(x[mm], y[mm], s=16, c=colors[gname], marker=markers.get(gname, 'o'), edgecolors='white',
                   linewidths=0.4, alpha=0.85, label=(legend_labels or {}).get(gname, gname), zorder=3)
    if m.any():
        lo, hi = min(x[m].min(), y[m].min()), max(x[m].max(), y[m].max())
        pad = 0.05 * (hi - lo if hi > lo else 1.)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=INK3, lw=1, ls='--', zorder=2)
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lo - pad, hi + pad)
        ss = np.sum((x[m] - x[m].mean()) ** 2)
        r2 = 1 - np.sum((y[m] - x[m]) ** 2) / ss if ss > 0 else np.nan
        ax.text(0.04, 0.96, f'R² = {r2:.3f}\nMAE = {np.mean(np.abs(y[m] - x[m])):.3g}', transform=ax.transAxes,
                va='top', ha='left', fontsize=8, color=INK2, zorder=5,
                bbox=dict(facecolor='white', alpha=0.85, edgecolor='none', pad=1.5))
    ax.set_xlabel(f'simulated {TARGET_LABEL[k]}')
    ax.set_ylabel(f'predicted {TARGET_LABEL[k]}')


def fig_parity(path, Z, P, groups, colors, markers, title, made, legend_labels=None):
    plt = _style()
    fig, axs = plt.subplots(1, 5, figsize=(17, 3.7))
    for j, k in enumerate(TARGETS):
        _parity_ax(axs[j], Z[:, j], P[:, j], groups, colors, markers, k, legend_labels)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=min(len(l), 7), bbox_to_anchor=(0.5, -0.06))
    fig.suptitle(title, x=0.01, ha='left', fontsize=11)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    _save(fig, path, made)


def _points_pred(ev, samples_by_id, points_by_id):
    import numpy as np
    import pandas as pd
    from fermentation_insights.coeff_predictor.train import pad_points
    from fermentation_insights.coeff_predictor.model import mpsp_equation_np
    s = samples_by_id.loc[ev['ids']].reset_index(drop=True)
    p = pd.concat([points_by_id[i] for i in ev['ids'] if i in points_by_id], ignore_index=True)
    Y, T, M, _, mask = pad_points(s, p)
    Mh = mpsp_equation_np(ev['coef'][:, None, :], ev['g'][:, None], Y, T)
    grp = np.repeat(s['biorefinery'].astype(str).to_numpy()[:, None], Y.shape[1], 1)
    return M[mask], Mh[mask], grp[mask]


def fig_mpsp_parity(path, panels, made):
    import numpy as np
    plt = _style()
    fig, axs = plt.subplots(1, len(panels), figsize=(5.2 * len(panels), 5.0), squeeze=False)
    seen = {}
    for ax, (title, (M, Mh, grp)) in zip(axs[0], panels):
        ok = np.isfinite(M) & np.isfinite(Mh) & (M > 0)
        if not ok.any():
            continue
        # axis range = simulated range (x/÷ 2); predictions outside are drawn on the border
        lo, hi = M[ok].min() / 2, M[ok].max() * 2
        out = ok & ((Mh < lo) | (Mh > hi))
        Mc = np.clip(np.where(Mh > 0, Mh, lo), lo, hi)
        for b in [b for b in BIO_ORDER if (grp == b).any()]:
            mm = ok & (grp == b)
            h = ax.scatter(M[mm], Mc[mm], s=7, c=BIO_COLOR[b], marker=BIO_MARKER[b], alpha=0.55, linewidths=0,
                           zorder=3)
            seen.setdefault(b, h)
        xs = np.array([lo, hi])
        ax.plot(xs, xs, color=INK3, lw=1, ls='--', zorder=2)
        for f in (0.8, 1.2):
            ax.plot(xs, f * xs, color=INK3, lw=0.8, ls=':', zorder=2)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ape = np.abs(Mh[ok] - M[ok]) / M[ok]
        ax.text(0.04, 0.96, f'MAPE = {100 * ape.mean():.1f} %\nmedian APE = {100 * np.median(ape):.1f} %\n'
                f'within ±20 % (dotted): {100 * np.mean(ape <= 0.2):.0f} %\n{ok.sum()} TY points'
                + (f'\n{out.sum()} predictions beyond the axes\n(drawn on the border)' if out.any() else ''),
                transform=ax.transAxes, va='top', fontsize=8, color=INK2, zorder=5,
                bbox=dict(facecolor='white', alpha=0.85, edgecolor='none', pad=1.5))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('simulated MPSP [USD/kg]')
        ax.set_ylabel('predicted MPSP [USD/kg]')
        ax.set_title(title, loc='left')
    if seen:
        leg = fig.legend(list(seen.values()), list(seen), loc='lower center', ncol=len(seen),
                         bbox_to_anchor=(0.5, -0.02), markerscale=2.5)
        for lh in leg.legend_handles:
            lh.set_alpha(1)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    _save(fig, path, made)


def fig_cv_bars(path, head, made):
    import numpy as np
    plt = _style()
    schemes = [s for s in SCHEMES if (head['scheme'] == s).any()]
    meths = [m for m in BASE_METHODS if (head['method'] == m).any()]
    fig, ax = plt.subplots(figsize=(11, 4.4))
    w = 0.8 / max(len(meths), 1)
    for i, m in enumerate(meths):
        vals, errs = [], []
        for s in schemes:
            r = head[(head['scheme'] == s) & (head['method'] == m)]
            vals.append(100 * r['mpsp_mape'].iloc[0] if len(r) else np.nan)
            errs.append(100 * r['mpsp_mape_se'].iloc[0] if len(r) else np.nan)
        x = np.arange(len(schemes)) + (i - (len(meths) - 1) / 2) * w
        ax.bar(x, vals, w * 0.92, color=METHOD_COLOR[m], label=METHOD_LABEL[m], zorder=3)
        ax.errorbar(x, vals, yerr=np.nan_to_num(errs), fmt='none', ecolor=INK2, elinewidth=0.8, capsize=2, zorder=4)
        if m == 'pinn':
            for xx, v in zip(x, vals):
                if np.isfinite(v):
                    ax.text(xx, v * 1.08, f'{v:.1f}', ha='center', va='bottom', fontsize=7, color=INK)
    for i, s in enumerate(schemes):
        r = head[(head['scheme'] == s) & (head['method'] == 'pinn')]
        fm = 100 * r['fit_mape'].iloc[0] if len(r) else np.nan
        ax.hlines(fm, i - 0.45, i + 0.45, colors=INK, linestyles='--', lw=1.2, zorder=5,
                  label='irreducible (per-sample fit MAPE)' if i == 0 else None)
    ax.set_yscale('log')
    ax.set_xticks(np.arange(len(schemes)))
    ax.set_xticklabels([SCHEME_LABEL[s] for s in schemes])
    ax.set_ylabel('MPSP MAPE on held-out TY points [%]\n(mean over folds ± s.e.)')
    ax.legend(ncol=3, loc='upper left', bbox_to_anchor=(0, 1.22))
    fig.tight_layout()
    _save(fig, path, made)


def fig_ablation(path, head, made):
    import numpy as np
    plt = _style()
    meths = ['pinn'] + list(ABLATIONS) + ['plain']
    schemes = [s for s in ABLATION_SCHEMES if (head['scheme'] == s).any()]
    fig, axs = plt.subplots(1, max(len(schemes), 1), figsize=(5.6 * max(len(schemes), 1), 3.8), squeeze=False)
    for ax, s in zip(axs[0], schemes):
        vals, errs, cols, labs = [], [], [], []
        for m in meths:
            r = head[(head['scheme'] == s) & (head['method'] == m)]
            if not len(r):
                continue
            vals.append(100 * r['mpsp_mape'].iloc[0])
            errs.append(100 * r['mpsp_mape_se'].iloc[0])
            cols.append(METHOD_COLOR[m])
            labs.append(METHOD_LABEL[m].replace('PI-NN ', 'PI-NN\n').replace('Plain MLP', 'Plain\nMLP'))
        y = np.arange(len(vals))[::-1]
        ax.barh(y, vals, 0.7, color=cols, zorder=3)
        e = np.nan_to_num(errs)
        ax.errorbar(vals, y, xerr=e, fmt='none', ecolor=INK2, elinewidth=0.8, capsize=2, zorder=4)
        for yy, v, ee in zip(y, vals, e):
            if np.isfinite(v):
                ax.text(v + ee, yy, f'  {v:.1f} %', va='center', fontsize=8, color=INK)
        r = head[(head['scheme'] == s) & (head['method'] == 'pinn')]
        if len(r) and np.isfinite(r['fit_mape'].iloc[0]):
            ax.axvline(100 * r['fit_mape'].iloc[0], color=INK, ls='--', lw=1.2, label='irreducible fit MAPE')
            ax.legend(loc='lower left', bbox_to_anchor=(0, 1.06))
        ax.set_yticks(y)
        ax.set_yticklabels(labs)
        ax.set_xlabel('MPSP MAPE on held-out TY points [%] (mean over folds ± s.e.)')
        ax.set_title(_slabel(s), loc='right')
        top = max([v + ee for v, ee in zip(vals, e) if np.isfinite(v)] or [1])
        ax.set_xlim(0, top * 1.25)
        ax.grid(axis='y', visible=False)
    fig.tight_layout()
    _save(fig, path, made)


def fig_contours(path, cases, made):
    """cases: [(title, sample_row, pred_coef(4), pred_g, points_df)]"""
    import numpy as np
    from matplotlib.colors import LogNorm, LinearSegmentedColormap
    from fermentation_insights.coeff_predictor.model import mpsp_equation_np
    plt = _style()
    cmap = LinearSegmentedColormap.from_list('blues', ['#cde2fb', '#86b6ef', '#3987e5', '#256abf', '#104281',
                                                       '#0d366b'])
    fig, axs = plt.subplots(len(cases), 2, figsize=(10.5, 4.4 * len(cases)), squeeze=False, layout='constrained')
    fig.get_layout_engine().set(h_pad=0.12, hspace=0.06)
    for row, (title, s, pc, pg, pts) in zip(axs, cases):
        # fitted (feasible) domain = range of the sample's feasible simulated points
        y0, y1 = pts['y'].min(), pts['y'].max()
        t0, t1 = pts['t'].min(), pts['t'].max()
        Y, T = np.meshgrid(np.linspace(y0, y1, 120), np.linspace(t0, t1, 120))
        Ms = mpsp_equation_np(np.array([float(s[k]) for k in 'abcd']), float(s['g']), Y, T)
        Mp = mpsp_equation_np(np.asarray(pc, float), float(pg), Y, T)
        vals = np.concatenate([Ms.ravel(), pts['MPSP'].to_numpy(float)])
        vals = vals[np.isfinite(vals) & (vals > 0)]
        lo, hi = np.quantile(vals, 0.01), np.quantile(vals, 0.99)
        levels = np.geomspace(lo, hi, 13)
        norm = LogNorm(lo, hi)
        Mhat_pts = mpsp_equation_np(np.asarray(pc, float), float(pg), pts['y'].to_numpy(), pts['t'].to_numpy())
        ape = np.abs(Mhat_pts - pts['MPSP'].to_numpy()) / pts['MPSP'].to_numpy()
        fit_ape = float(s['fit_mape'])
        for ax, Z, lab in ((row[0], Ms, f'simulated (per-sample fit; fit MAPE {100 * fit_ape:.1f} %)'),
                           (row[1], Mp, f'predicted (PI-NN; MAPE {100 * np.mean(ape):.1f} % at TY points)')):
            cf = ax.contourf(Y, T, np.clip(Z, lo, hi), levels=levels, cmap=cmap, norm=norm, extend='neither')
            cs = ax.contour(Y, T, Z, levels=levels[::2], colors='white', linewidths=0.7)
            ax.clabel(cs, fmt=lambda v: f'{v:.2f}', fontsize=7, colors='white')
            if ax is row[1]:
                ax.contour(Y, T, Ms, levels=levels[::2], colors=INK, linewidths=0.8, linestyles='--')
            ax.scatter(pts['y'], pts['t'], c=pts['MPSP'], cmap=cmap, norm=norm, s=30, edgecolors='white',
                       linewidths=1.0, zorder=4)
            ax.set_xlabel('yield y [g/g glucose-eq.]')
            ax.set_ylabel('titer t [g/L]')
            ax.set_title(f'{title}\n{lab}', loc='left', fontsize=9)
            ax.grid(False)
        cb = fig.colorbar(cf, ax=list(row), shrink=0.9, pad=0.02)
        cb.minorticks_off()
        cb.set_ticks(levels[::2])
        cb.set_ticklabels([f'{v:.2f}' for v in levels[::2]])
        cb.set_label('MPSP [USD/kg]')
    fig.suptitle('Held-out test configurations: MPSP over the feasible TY domain. Dots = simulated TY points '
                 '(colour = simulated MPSP);\ndashed black lines (right) = simulated surface for comparison.',
                 x=0.01, ha='left', fontsize=9, color=INK2)
    _save(fig, path, made)


def fig_pub(path, pub, made):
    import numpy as np
    plt = _style()
    fig, axs = plt.subplots(1, 5, figsize=(17, 3.7))
    for ax, k in zip(axs, TARGETS):
        x, y = pub[f'{k}_pub'].to_numpy(float), pub[f'{k}_ours'].to_numpy(float)
        groups = pub['biorefinery'].astype(str).to_numpy()
        for b in [b for b in BIO_ORDER if (groups == b).any()]:
            m = groups == b
            ax.scatter(x[m], y[m], s=26, c=BIO_COLOR[b], marker=BIO_MARKER[b], edgecolors='white', linewidths=0.5,
                       label=b, zorder=3)
        f = np.isfinite(x) & np.isfinite(y)
        if f.any():
            lo, hi = min(x[f].min(), y[f].min()), max(x[f].max(), y[f].max())
            pad = 0.05 * (hi - lo if hi > lo else 1)
            ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=INK3, ls='--', lw=1)
            med = np.nanmedian(pub[f'{k}_ratio'])
            ax.text(0.04, 0.96, f'median ratio = {med:.3f}\nn = {f.sum()}', transform=ax.transAxes, va='top',
                    fontsize=8, color=INK2)
        ax.set_xlabel(f'publication {TARGET_LABEL[k]}')
        ax.set_ylabel(f'this work {TARGET_LABEL[k]}')
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=len(l), bbox_to_anchor=(0.5, -0.06))
    fig.suptitle('Baseline coefficients: this work (current BioSTEAM, harness v2) vs. publication',
                 x=0.01, ha='left', fontsize=11)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    _save(fig, path, made)


def fig_sens(path, pts, made):
    import numpy as np
    plt = _style()
    sel = [('fs', 'b'), ('el', 'a')]
    meths = [m for m in ('pinn', 'plain', 'hgb') if (pts['method'] == m).any()] if len(pts) else []
    fig, axs = plt.subplots(len(sel), max(len(meths), 1), figsize=(4.4 * max(len(meths), 1), 3.9 * len(sel)),
                            squeeze=False)
    for i, (key, k) in enumerate(sel):
        for jm, m in enumerate(meths):
            ax = axs[i, jm]
            d = pts[(pts['method'] == m) & (pts['price'] == key) & (pts['coefficient'] == k)]
            if not len(d):
                ax.set_visible(False)
                continue
            x, y = d['simulated'].to_numpy(float), d['predicted'].to_numpy(float)
            ax.scatter(x, y, s=16, c=METHOD_COLOR[m], edgecolors='white', linewidths=0.4, zorder=3)
            lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
            pad = 0.05 * (hi - lo if hi > lo else 1)
            ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=INK3, ls='--', lw=1)
            ss = np.sum((x - x.mean()) ** 2)
            r2 = 1 - np.sum((y - x) ** 2) / ss if ss > 0 else np.nan
            ax.text(0.04, 0.96, f'R² = {r2:.3f}', transform=ax.transAxes, va='top', fontsize=8, color=INK2)
            unit = PRICE_INFO[key][2]
            ax.set_xlabel(f"simulated intensity int_{key},{k}′\n(= ∂{k}′/∂p, {k}′ units per {unit})")
            ax.set_ylabel(f"predicted ∂{k}′/∂p_{key}")
            ax.set_title(f'{METHOD_LABEL[m]}: {PRICE_INFO[key][1]} price', loc='left')
    fig.tight_layout()
    _save(fig, path, made)


def fig_channels(path, ch, made):
    import numpy as np
    plt = _style()
    cfgs = list(dict.fromkeys(ch['config_id']))
    fig, axs = plt.subplots(len(cfgs), 4, figsize=(17, 3.6 * len(cfgs)), squeeze=False)
    for i, cid in enumerate(cfgs):
        for j, k in enumerate(TARGETS[:4]):
            ax = axs[i, j]
            d = ch[(ch['config_id'] == cid) & (ch['coefficient'] == k)]
            y = np.arange(len(d))[::-1]
            ax.barh(y + 0.19, d['simulated_mean'], 0.36, color=PALETTE[0], label='simulated', zorder=3)
            ax.barh(y - 0.19, d['predicted_mean'], 0.36, color=PALETTE[1], label='predicted (PI-NN)', zorder=3)
            ax.axvline(0, color=INK3, lw=0.8)
            ax.set_yticks(y)
            ax.set_yticklabels(d['channel'] if j == 0 else [''] * len(d))
            ax.set_title(f'{cid}: {TARGET_LABEL[k]}' if j == 0 else TARGET_LABEL[k], loc='left', fontsize=9)
            ax.grid(axis='y', visible=False)
            if i == 0 and j == 0:
                ax.legend(loc='lower right')
    fig.suptitle('Channel breakdown of the coefficients, held-out test configurations (mean over samples)',
                 x=0.01, ha='left', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    _save(fig, path, made)


def fig_hpo(path, tab, made):
    import numpy as np
    plt = _style()
    k, kse = ('criterion', 'criterion_se') if 'criterion' in tab else ('mpsp_mape', 'mpsp_mape_se')
    t = tab[np.isfinite(tab[k])].sort_values(k)
    fig, ax = plt.subplots(figsize=(7, 0.32 * len(t) + 1.4))
    y = np.arange(len(t))[::-1]
    cols = [PALETTE[0] if s else '#b8b7b1' for s in t['selected']]
    ax.barh(y, 100 * t[k], 0.7, color=cols, zorder=3)
    ax.errorbar(100 * t[k], y, xerr=100 * np.nan_to_num(t[kse]), fmt='none', ecolor=INK2,
                elinewidth=0.8, capsize=2, zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels([f"#{int(r.candidate)}  w{int(r.width)} d{int(r.depth)} lr {r.lr:.1e} wd {r.weight_decay:.0e} "
                        f"w_ch {r.w_ch:g} w_phys {r.w_phys:g} ES {r.early_stop}" for r in t.itertuples()], fontsize=7)
    ax.set_xlabel('configuration-balanced held-out MPSP MAPE [%] (fold mean, grouped CV by configuration)')
    ax.set_title('Hyperparameter search (selected in blue)', loc='left')
    ax.grid(axis='y', visible=False)
    fig.tight_layout()
    _save(fig, path, made)


# %% Stage 8: report

def _slabel(s):
    return SCHEME_LABEL[s].replace('-\n', '-').replace('\n', ' ')


def _fmt(v, pct=False):
    import numpy as np
    if isinstance(v, (bool, np.bool_)):
        return 'yes' if v else 'no'
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        if not np.isfinite(v):
            return '–' if np.isnan(v) else ('∞' if v > 0 else '−∞')
        if pct:
            return f'{100 * v:.1f}'
        av = abs(v)
        return f'{v:.3g}' if (av >= 1e-3 and av < 1e5) or v == 0 else f'{v:.2e}'
    return str(v)


def md_table(df, pct_cols=()):
    cols = list(df.columns)
    out = ['| ' + ' | '.join(map(str, cols)) + ' |', '|' + '|'.join('---' for _ in cols) + '|']
    for r in df.itertuples(index=False):
        out.append('| ' + ' | '.join(_fmt(v, c in pct_cols) for c, v in zip(cols, r)) + ' |')
    return '\n'.join(out)


def write_report(a, ctx, D, sp, dev, plan, hpo, cv, fin, checks, figs, timings):
    import numpy as np
    import pandas as pd
    L = []
    P = L.append
    mode = 'QUICK (debug: tiny search, 2 members, few epochs -- numbers are NOT results)' if a.quick else 'full'
    P('# Coefficient predictor: validation report\n')
    P(f'Generated by `experiments.py` ({mode}) on {time.strftime("%Y-%m-%d %H:%M")}; '
      f'data hash `{ctx["data_hash"]}`, split hash `{ctx["split_hash"]}`, seed {a.seed}, '
      f'{a.threads} worker(s). Command: `{" ".join(ctx["argv"])}`.\n')
    P('Features: ' + ('structured downstream specification (`FeatureSpec(downstream=True)`, DESIGN.md section 5)'
                      if ctx.get('spec_kw', {}).get('downstream') else 'without the downstream specification')
      + ('; PI-NN with the water-removal gate (`ModelConfig.water_gate`)'
         if ctx.get('mk_extra', {}).get('water_gate') else '')
      + ('; PI-NN with channel-wise input masking (`ModelConfig.channel_inputs='
         f"{ctx['mk_extra']['channel_inputs']!r}`, DESIGN.md section 3)"
         if ctx.get('mk_extra', {}).get('channel_inputs') else '')
      + ('; PI-NN with the hard productivity structure (`ModelConfig.r_structure='
         f"{ctx['mk_extra']['r_structure']!r}`, DESIGN.md section 3)"
         if ctx.get('mk_extra', {}).get('r_structure') else '')
      + ('; extrapolation safeguards (sugar-availability prior `fs_prior`, neutralization gate, '
         'feedstock-input clipping `clip_groups`)' if ctx.get('mk_extra', {}).get('fs_prior') else '') + '.\n')
    if D['recent_raw_files']:
        P(f'> **Partial data:** {len(D["recent_raw_files"])} raw file(s) were still being written when the tables '
          'were built; the numbers describe the data available at that time.\n')
    P('**Where numbers come from.** (1) *CV estimates*: models re-trained inside each fold on parts of the '
      'development set only (selected hyperparameters, '
      f'{plan["cv_members"]}-member ensembles for the networks). (2) *Final-test estimate*: the selected PI-NN '
      f'ensemble ({plan["final_members"]} members, `models/final_dev`) and the baselines, trained once on the whole '
      'development set and evaluated once on the final test set, which played no role in any choice. (3) *Release '
      f'model*: {plan["release_members"]}-member PI-NN ensemble with the same hyperparameters trained on ALL data '
      '(`models/release`); it has no independent test estimate -- the CV and final-test numbers are the estimates '
      'of its accuracy (slightly pessimistic: it is trained on more data). MAPE = mean absolute percentage error of MPSP over the held-out samples\' simulated TY points; '
      '"fit MAPE" is the irreducible error of the per-sample 4-coefficient fit on the same points.\n')
    # summary
    head, fh = cv['headline'], fin['headline']
    P('## Summary\n')

    def _m(sch, meth, col='mpsp_mape'):
        r = head[(head['scheme'] == sch) & (head['method'] == meth)]
        return float(r[col].iloc[0]) if len(r) else float('nan')

    for sch in ('product', 'config', 'feedstock', 'random', 'productivity'):
        if not (head['scheme'] == sch).any():
            continue
        base = {m: _m(sch, m) for m in BASE_METHODS if m != 'pinn'}
        bm = min(base, key=lambda m: base[m] if np.isfinite(base[m]) else np.inf)
        se = _m(sch, 'pinn', 'mpsp_mape_se')
        P(f'* {_slabel(sch)} (CV): PI-NN MPSP MAPE {100 * _m(sch, "pinn"):.1f} %'
          + (f' (± {100 * se:.1f} s.e.)' if np.isfinite(se) else ' (single fold)')
          + f' vs. best baseline {METHOD_LABEL[bm]} '
          f'{100 * base[bm]:.1f} %; irreducible fit MAPE {100 * _m(sch, "pinn", "fit_mape"):.1f} %.')
    for sch in ABLATION_SCHEMES:
        d = {m: _m(sch, m) - _m(sch, 'pinn') for m in ABLATIONS if (head['method'] == m).any()}
        if d:
            P(f'* Ablations ({_slabel(sch)}), change in MAPE vs. full PI-NN: '
              + ', '.join(f'{METHOD_LABEL[m].replace("PI-NN ", "")} {100 * v:+.1f} pp' for m, v in d.items()) + '.')
    if len(fh):
        for st, lab in (('test_ood', 'held-out configurations'), ('test_id', 'in-distribution samples')):
            r = fh[(fh['set'] == st) & (fh['method'] == 'pinn')]
            if len(r):
                P(f'* Final test, {lab} (PI-NN trained on dev): MPSP MAPE {100 * r["mpsp_mape"].iloc[0]:.1f} % '
                  f'(median {100 * r["mpsp_mape_median"].iloc[0]:.1f} %; fit MAPE {100 * r["fit_mape"].iloc[0]:.1f} %; '
                  f'n = {int(r["n_test"].iloc[0])} samples).')
    P('')
    # data
    ds = checks['data_summary']
    rep = D['report']
    P('## 1. Data\n')
    P(f'{len(D["samples"])} samples ({len(D["points"])} feasible TY points) of {ds["config_id"].nunique()} '
      f'configurations / {ds["biorefinery"].nunique()} biorefineries from {D["n_raw_files"]} raw files. '
      f'Rejected samples: {rep.get("rejections", {})}. Median fit R² {np.nanmedian(D["samples"]["R2"]):.4f}, '
      f'median fit MAPE {100 * np.nanmedian(D["samples"]["fit_mape"]):.2f} %. Per configuration: '
      '`tables/data_summary.csv`.\n')
    bio = ds.groupby('biorefinery').agg(configs=('config_id', 'size'), samples=('n_samples', 'sum'),
                                        points=('n_points', 'sum'), median_fit_mape=('median_fit_mape', 'median'))
    bio = bio.reindex([b for b in BIO_ORDER if b in bio.index] + [b for b in bio.index if b not in BIO_ORDER])
    P(md_table(bio.reset_index(), pct_cols=('median_fit_mape',)) + '\n\n(median_fit_mape in %)\n')
    nf = checks.get('noise_floor') or {}
    if nf.get('n_configs'):
        P(f'Noise floor of random samples (unobserved biorefinery-specific `bp_*` draws; '
          f'`evaluate.noise_floor`, within-configuration 5-fold ridge on the observable inputs, '
          f'{nf["n_configs"]} development configurations with >= 10 random samples): MPSP MAPE '
          f'{100 * nf["floor_mape_cfg_gmean"]:.1f} % (configuration geo-mean; mean '
          f'{100 * nf["floor_mape_cfg_mean"]:.1f} %), vs. {100 * nf["const_mape_cfg_gmean"]:.1f} % for the '
          f'configuration mean and {100 * nf["fit_mape_cfg_gmean"]:.1f} % for the 4-coefficient fit itself. '
          'Spec-determined samples (baseline / productivity: every `bp_*` at baseline) have no such noise; '
          'the CV tables report both kinds separately (`*_spec`, `*_random`). `tables/noise_floor.csv`.\n')
    flags = checks.get('data_flags')
    if flags is not None and len(flags):
        P('Data-consistency flags (configurations to inspect; they are kept): '
          + '; '.join(f'`{r.config_id}`: {r.flag}' for r in flags.itertuples()) + '.\n')
    # split
    P('## 2. Final test set (fixed before modelling)\n')
    P(f'Rule: {sp["rule"]} (seed {sp["seed"]}). Held-out configurations (out-of-distribution test, '
      f'{sp["n_test_ood"]} samples): ' + ', '.join(f'`{c["config_id"]}` ({c["n_samples"]})' for c in sp['test_configs'])
      + f'; in-distribution test: {sp["n_test_id"]} samples of the other configurations. Development set: '
      f'{sp["n_dev"]} samples, {sp["n_dev_configs"]} configurations. '
      + (f'Biorefineries with too few configurations for a held-out one: '
         f'{", ".join(sp["biorefineries_without_held_out_config"])}. ' if sp['biorefineries_without_held_out_config']
         else '') + '`split.json`.\n')
    # hpo
    t = hpo['table']
    P('## 3. Hyperparameter selection (development set only)\n')
    P(f'Random search over {len(t)} candidates (candidate 0 = package defaults) with {plan["hpo_k"]}-fold '
      f'leave-configurations-out CV (folds stratified by biorefinery; {plan.get("hpo_members", 1)} network(s) per '
      'fit, early stopping on held-out configurations), criterion = fold mean of the configuration-balanced '
      'held-out MPSP MAPE (`criterion`; each held-out configuration counts once). '
      + (f'Search size from the measured cost of candidate 0 and the {a.budget_hours} h budget '
         f'(estimate {plan["estimate"]["total_hours"]:.2f} h for the whole run on {a.threads} workers). '
         if plan.get('estimate') else '') + '\n')
    show = t[['candidate', 'width', 'depth', 'dropout', 'lr', 'weight_decay', 'w_mpsp', 'w_ch', 'w_phys',
              'colloc_ratio', 'early_stop', 'clip_inputs', 'criterion', 'criterion_se', 'mpsp_mape',
              'mpsp_mape_cfg_gmean', 'fit_mape', 'epochs', 'seconds_per_fit', 'selected']]
    P(md_table(show.sort_values('criterion'), pct_cols=('criterion', 'criterion_se', 'mpsp_mape',
                                                        'mpsp_mape_cfg_gmean', 'fit_mape'))
      + '\n\n(MAPE columns in %.) **Selected:** ' + ', '.join(f'{k}={v}' for k, v in hpo['selected'].items())
      + '. `tables/hpo_candidates.csv`, `figures/hpo_search.png`.\n')
    # cv
    head = cv['headline']
    P('## 4. Generalization CV estimates (development set)\n')
    P(f'Selected hyperparameters; inner calibration / early-stopping hold-outs grouped by {INNER_GROUP}. '
      f'Random K = {plan["random_k"]}; leave-one-configuration-out over '
      + (f'a stratified subset of {plan["loco_folds"]} configurations (budget); ' if plan.get('loco_folds') else
         'all development configurations; ')
      + 'productivity: train on the central 20-80 % quantile range of log(r / baseline r), test both tails. '
      'MAPE = mean over folds (± s.e.); pooled = over all out-of-fold samples; R² pooled out-of-fold.\n')
    cols = ['scheme', 'method', 'n_folds_common', 'mpsp_mape', 'mpsp_mape_se', 'mpsp_mape_median', 'mpsp_mape_pooled',
            'fit_mape'] + [f'R2_{k}' for k in TARGETS]
    P(md_table(head[cols], pct_cols=('mpsp_mape', 'mpsp_mape_se', 'mpsp_mape_median', 'mpsp_mape_pooled',
                                     'fit_mape')) + '\n\n(MAPE columns in %.)\n')
    P('Per-coefficient MAE (pooled out-of-fold, natural units) and constraint-violation / coverage rates:\n')
    vcols = ['scheme', 'method'] + [f'MAE_{k}' for k in TARGETS] + ['viol_mono_y', 'viol_mono_t', 'viol_positive',
                                                                     'viol_mono_r', 'viol_g_le_gmax', 'coverage_MPSP']
    P(md_table(head[vcols]) + '\n')
    lp = checks.get('lopo_per_product')
    if lp is not None and len(lp):
        P('Leave-one-product-out, per held-out product (MAPE %):\n')
        P(md_table(lp, pct_cols=tuple(c for c in lp.columns if c not in ('product',))) + '\n')
    P('Full tidy results: `tables/cv_results.csv`; summary: `tables/cv_headline.csv`; out-of-fold predictions: '
      '`tables/cv_oof_predictions.csv.gz`.\n')
    # final
    fh = fin['headline']
    P('## 5. Final test estimate (models trained on the development set, evaluated once)\n')
    if len(fh):
        cols = ['set', 'method', 'n_train', 'n_test', 'mpsp_mape', 'mpsp_mape_median', 'fit_mape'] + \
               [f'R2_{k}' for k in TARGETS] + ['coverage_MPSP', 'viol_mono_y', 'viol_mono_t']
        P(md_table(fh[fh['set'] != 'train_all'][cols], pct_cols=('mpsp_mape', 'mpsp_mape_median', 'fit_mape'))
          + '\n\n(test_ood = held-out configurations, test_id = in-distribution samples; MAPE in %.)\n')
        conf = (fin['final'].get('pinn') or {}).get('conformal') or {}
        if conf:
            P(f'PI-NN conformal calibration: {conf.get("n")} samples / {conf.get("n_calib_groups")} held-out '
              f'configurations, level {1 - conf.get("alpha", 0.1):.2f}, finite-sample valid: '
              f'{conf.get("finite_sample_valid")}, coverage guaranteed (>= 10 groups): '
              f'{conf.get("coverage_guaranteed")}.\n')
        pc = checks.get('test_per_config')
        if pc is not None and len(pc):
            P('PI-NN per held-out test configuration:\n')
            P(md_table(pc, pct_cols=('mpsp_mape', 'mpsp_mape_median', 'fit_mape')) + '\n')
    else:
        P('Final models failed (see log).\n')
    # release
    P('## 6. Release model\n')
    rel = fin['release']
    if rel['error']:
        P(f'Release training FAILED: {rel["error"]}\n')
    else:
        conf = fin.get('release_conformal') or {}
        rc = fin.get('release_check') or {}
        ins = fh[(fh['set'] == 'train_all')]
        ins_txt = f'{100 * ins["mpsp_mape"].iloc[0]:.1f} %' if len(ins) else 'n/a'
        P(f'`{os.path.relpath(rel["saved_model"], PKG_DIR)}`: {plan["release_members"]} members trained on all '
          f'{len(D["samples"])} samples (development + test; each member early-stops on its own held-out '
          f'configurations), epochs {rel["epochs"]}. In-sample MPSP MAPE (NOT a generalization estimate): '
          f'{ins_txt}.\n')
        if conf:
            P(f'Intervals: split-conformal quantiles transferred from the dev-trained ensemble, scored on '
              f'{conf.get("n")} samples of {conf.get("n_calib_groups")} configurations it never saw (its calibration '
              f'configurations + the test configurations); level {1 - conf.get("alpha", 0.1):.2f} '
              f'(finite-sample level reached for MPSP: {_fmt(conf.get("mpsp_level", float("nan")))}). Valid for the dev '
              'ensemble; for the release ensemble (more training data) an approximation, flagged '
              '`coverage_guaranteed = False` in its checkpoint.\n')
        else:
            P('Intervals: none (conformal transfer failed; see log).\n')
        P(f'Loads with `CoefficientPredictor.load`: {rc.get("ok")} '
          + (f'(e.g. {rc.get("example_config")} baseline: ' + ', '.join(f'{k}={v:.4g}' for k, v in
                                                                         rc.get('example', {}).items())
             + f'; intervals reported: {rc.get("intervals")}).' if rc.get('ok') else f'({rc.get("error")}).') + '\n')
    # checks
    P('## 7. Physical-consistency and interpretability checks\n')
    pub = checks['pub']
    P('**(i) Baseline coefficients vs. publication** (data only, no model; `sample_kind == baseline` of each '
      'configuration vs. `results/data/<config>_coefficients.npy`; differences are expected from the newer '
      f'BioSTEAM / harness). {len(pub)} of 32 publication configurations present.\n')
    if len(pub):
        med = {k: np.nanmedian(pub[f'{k}_ratio']) for k in TARGETS}
        P('Median ratio ours/publication: ' + ', '.join(f"{k}′ {v:.3f}" if k != 'g' else f'g {v:.3f}'
                                                     for k, v in med.items()) + '.\n')
        cols = ['config_id'] + [f'{k}_{s}' for k in TARGETS for s in ('ours', 'pub', 'ratio')]
        P(md_table(pub[cols]) + '\n')
    sens = checks['sens']
    P('**(ii) Price sensitivities** (final dev-trained models on the test set): central secant slopes '
      f'(±{100 * SENS_STEP:.0f} % price) of the predicted coefficients vs. the simulated price intensities `int__<fs|el|ng>__<k>` '
      '(exact ∂θ/∂p at fixed design). The PI-NN is linear in price by construction.\n')
    if len(sens):
        P('`fd_vs_intensity_rel` (PI-NN) = max |secant slope − predicted intensity| / max |intensity|: ≈ 0 up to '
          'float32 round-off of the total coefficient = exactly linear in price, as designed.\n')
        P(md_table(sens[sens['coefficient'].isin(['a', 'b', 'c'])].reset_index(drop=True)) + '\n')
    ch = checks['channels']
    P('**(iii) Channel breakdown** of the held-out test configurations (final PI-NN; mean over each configuration\'s '
      'test samples): `tables/channel_breakdown.csv`, `figures/channels_test.png`.\n')
    if len(ch):
        agg = []
        for k, d in ch.groupby('coefficient'):
            x, y = d['simulated_mean'].to_numpy(float), d['predicted_mean'].to_numpy(float)
            f = np.isfinite(x) & np.isfinite(y)
            agg.append(dict(coefficient=k, channels=d['channel'].nunique(), configs=d['config_id'].nunique(),
                            mean_abs_err=float(np.nanmean(d['MAE'])),
                            corr_sim_pred=float(np.corrcoef(x[f], y[f])[0, 1]) if f.sum() > 2 else np.nan))
        P('Agreement over (configuration × channel) pairs:\n')
        P(md_table(pd.DataFrame(agg)) + '\n')
    # figures & timings
    P('## 8. Figures\n')
    P('\n'.join(f'* `{os.path.relpath(f, a.out)}`' for f in figs) + '\n')
    P('## 9. Reproduction and timings\n')
    P(f'`{" ".join(ctx["argv"])}` (PYTHONPATH=/home/user/fermentation_insights). Python {platform.python_version()}, '
      f'{ctx["versions"]}. Every fit is cached under `cache/` (keyed by data hash, split, job and code version).\n')
    P(md_table(timings) + '\n')
    path = os.path.join(a.out, 'REPORT.md')
    with open(path, 'w') as fh_:
        fh_.write('\n'.join(L))
    return path


# %% Main

def parse_args(argv=None):
    ap = argparse.ArgumentParser(description='Reproducible validation experiments of the coefficient predictor.')
    ap.add_argument('--quick', action='store_true', help='debug run: tiny search, 2 members, few epochs')
    ap.add_argument('--raw', nargs='+', default=[DEFAULT_RAW], help='raw harness data dir(s)')
    ap.add_argument('--processed', default=DEFAULT_PROCESSED, help='processed-data cache dir')
    ap.add_argument('--rebuild', action='store_true', help='rebuild the processed tables even if raw is unchanged')
    ap.add_argument('--no-rebuild', action='store_true', help='use existing processed tables even if raw changed')
    ap.add_argument('--out', default=None, help='results dir (default results/, quick: results/quick)')
    ap.add_argument('--models', default=None, help='models dir (default models/, quick: models/quick)')
    ap.add_argument('--threads', type=int, default=1, help='parallel single-threaded worker processes')
    ap.add_argument('--budget-hours', type=float, default=2.5, help='full mode: wall-time budget of the whole run')
    ap.add_argument('--hpo-candidates', type=int, default=0, help='override the budget-derived search size')
    ap.add_argument('--loco-folds', type=int, default=0, help='override: number of LOCO folds (0 = budget / all)')
    ap.add_argument('--members', type=int, default=3, help='full mode: ensemble members in CV (3-5)')
    ap.add_argument('--hpo-members', type=int, default=2, help='full mode: ensemble members per HPO fit')
    ap.add_argument('--split-file', default=SPLIT_FILE,
                    help='frozen final test split (created from PLANNED_CONFIGS if missing, then never changed)')
    ap.add_argument('--min-test-configs', type=int, default=None,
                    help='only when the split is created: biorefineries with >= this many planned configurations '
                         'get a held-out test configuration (default 4)')
    ap.add_argument('--id-test-frac', type=float, default=None, help='only when the split is created (default 0.10)')
    ap.add_argument('--seed', type=int, default=0)
    # on by default since the 2026-10-17 development CV (DESIGN.md section 6)
    ap.add_argument('--no-downstream', action='store_true',
                    help='do not use the structured downstream specification features '
                         '(FeatureSpec(downstream=True); implies --no-water-gate)')
    ap.add_argument('--no-water-gate', action='store_true',
                    help='PI-NN without the water-removal gate (ModelConfig.water_gate)')
    # on by default since the 2026-10-17 masking CV (DESIGN.md section 3)
    ap.add_argument('--no-masking', action='store_true',
                    help="PI-NN with one shared encoder instead of the channel-wise input masking "
                         "(ModelConfig.channel_inputs='physics_v2')")
    # opt-in (2026-10-17 round 4 development CV, DESIGN.md section 3)
    ap.add_argument('--r-structure', action='store_true',
                    help="PI-NN with the hard productivity structure (ModelConfig.r_structure='physics': "
                         "fermentation capital and electricity intensity = theta0(x) + theta1(x) r_ref/r)")
    # on by default since the round-4 extrapolation CV (DESIGN.md section 6)
    ap.add_argument('--no-extrapolation-safeguards', action='store_true',
                    help='PI-NN without the sugar-availability prior of the feedstock intensity '
                         '(ModelConfig.fs_prior), the neutralization gate (ModelConfig.neutralization_gate) '
                         "and the clipping of the feedstock inputs (ModelConfig.clip_groups)")
    ap.add_argument('--split-seed', type=int, default=None, help='only when the split is created (default 2026)')
    ap.add_argument('--replan', action='store_true', help='recompute plan.json from the pilot timings')
    ap.add_argument('--fresh', action='store_true', help='delete the fit cache first')
    ap.add_argument('--plan-only', action='store_true',
                    help='run data, split and the HPO pilot (candidate 0), write plan.json with the budget '
                         'estimate and stop')
    ap.add_argument('--force-pool', action='store_true', help=argparse.SUPPRESS)  # debug: pool even with 1 worker
    ap.add_argument('--nice', type=int, default=0, help='increase the process niceness')
    a = ap.parse_args(argv)
    if a.out is None:
        a.out = os.path.join(DEFAULT_RESULTS, 'quick') if a.quick else DEFAULT_RESULTS
    if a.models is None:
        a.models = os.path.join(DEFAULT_MODELS, 'quick') if a.quick else DEFAULT_MODELS
    a.split_seed_given = any(v is not None for v in (a.split_seed, a.min_test_configs, a.id_test_frac))
    a.split_seed = 2026 if a.split_seed is None else a.split_seed
    a.min_test_configs = 4 if a.min_test_configs is None else a.min_test_configs
    a.id_test_frac = 0.10 if a.id_test_frac is None else a.id_test_frac
    return a


def main(argv=None):
    a = parse_args(argv)
    _set_thread_env()
    if a.nice:
        os.nice(a.nice)
    os.makedirs(a.out, exist_ok=True)
    for d in ('tables', 'figures', 'cache'):
        os.makedirs(os.path.join(a.out, d), exist_ok=True)
    os.makedirs(a.models, exist_ok=True)
    if a.fresh:
        shutil.rmtree(os.path.join(a.out, 'cache'), ignore_errors=True)
        os.makedirs(os.path.join(a.out, 'cache'))
    _LOG['fh'] = open(os.path.join(a.out, 'experiments.log'), 'a')
    t_start = time.time()
    stage_t = {}
    log(f'===== experiments.py {"QUICK" if a.quick else "FULL"} mode; out={a.out}; threads={a.threads}')
    warnings.filterwarnings('ignore', category=UserWarning)
    import numpy as np
    import pandas as pd
    import torch
    import sklearn
    import matplotlib
    torch.set_num_threads(1)
    ctx = dict(argv=['python -m ' + MODULE] + (sys.argv[1:] if argv is None else list(argv)),
               versions=f'torch {torch.__version__}, numpy {np.__version__}, pandas {pd.__version__}, '
                        f'scikit-learn {sklearn.__version__}, matplotlib {matplotlib.__version__}')
    mk_extra = {}
    if not a.no_downstream:
        ctx['spec_kw'] = {'downstream': True}
        if not a.no_water_gate:
            mk_extra['water_gate'] = True
        if not a.no_extrapolation_safeguards:
            mk_extra.update(EXTRAPOLATION_SAFEGUARDS)
    if not a.no_masking:
        mk_extra['channel_inputs'] = 'physics_v2'
    if getattr(a, 'r_structure', False):
        mk_extra['r_structure'] = 'physics'
    if mk_extra:
        ctx['mk_extra'] = mk_extra
    T = lambda name, t0: stage_t.__setitem__(name, time.time() - t0)

    t0 = time.time()
    D = stage_data(a)
    T('data', t0)
    ctx['data_hash'] = D['hash']
    t0 = time.time()
    sp, dev = stage_split(a, D)
    ctx['split_hash'] = sp['hash']
    T('split', t0)
    plan = base_plan(a)
    if a.hpo_candidates:
        plan['hpo_candidates'] = a.hpo_candidates
    if a.loco_folds:
        plan['loco_folds'] = a.loco_folds
    runner = Runner(a.threads, a.processed, D['hash'], os.path.join(a.out, 'cache'), a.force_pool)
    try:
        t0 = time.time()
        hpo = stage_hpo(a, plan, runner, dev, ctx, sp['n_test_ood'] + sp['n_test_id'], len(D['samples']))
        T('hpo', t0)
        _dump_json(dict(plan, data_hash=ctx['data_hash'], split_hash=ctx['split_hash'], threads=a.threads,
                        selected=hpo['selected'] if hpo else None), os.path.join(a.out, 'plan.json'))
        if hpo is None:
            log(f'plan-only: plan written to {os.path.join(a.out, "plan.json")}; stopping')
            return 0
        t0 = time.time()
        cv = stage_cv(a, plan, runner, dev, ctx, hpo['selected'], D)
        T('cv', t0)
        t0 = time.time()
        fin = stage_final(a, plan, runner, dev, ctx, hpo['selected'], sp, D, a.models)
        T('final+release', t0)
    finally:
        runner.close()

    # ---- tables
    t0 = time.time()
    tdir, fdir = os.path.join(a.out, 'tables'), os.path.join(a.out, 'figures')
    S = D['samples'].assign(sample_id=D['samples']['sample_id'].astype(str)).set_index('sample_id', drop=False)
    Pts = {k: v for k, v in D['points'].assign(sample_id=D['points']['sample_id'].astype(str)).groupby('sample_id')}
    checks = dict(data_summary=data_summary(D['samples'], D['points']), data_flags=data_flags(D['samples']))
    checks['data_summary'].to_csv(os.path.join(tdir, 'data_summary.csv'), index=False)
    hpo['table'].to_csv(os.path.join(tdir, 'hpo_candidates.csv'), index=False)
    hpo['folds'].to_csv(os.path.join(tdir, 'hpo_folds.csv'), index=False)
    _dump_json(dict(selected=hpo['selected'], candidate=hpo['best'], hpo_folds=hpo['fold_info'], plan=plan),
               os.path.join(tdir, 'selected_hyperparameters.json'))
    cv['results'].to_csv(os.path.join(tdir, 'cv_results.csv'), index=False)
    # noise floor of random samples (unobserved bp_* parameters; development set only)
    from fermentation_insights.coeff_predictor.evaluate import noise_floor
    nf, nf_summ = noise_floor(dev, D['points'])
    nf.to_csv(os.path.join(tdir, 'noise_floor.csv'), index=False)
    _dump_json(nf_summ, os.path.join(tdir, 'noise_floor_summary.json'))
    checks['noise_floor'] = nf_summ
    cv['headline'] = headline(cv['results'], METHOD_ORDER)
    cv['headline'].to_csv(os.path.join(tdir, 'cv_headline.csv'), index=False)
    if len(cv['oof']):
        cv['oof'].to_csv(os.path.join(tdir, 'cv_oof_predictions.csv.gz'), index=False)
    fin['table'].to_csv(os.path.join(tdir, 'final_test_results.csv'), index=False)
    fin['headline'] = final_headline(fin['table'])
    fin['headline'].to_csv(os.path.join(tdir, 'final_test_headline.csv'), index=False)
    res = cv['results']
    lp = res[(res['scheme'] == 'product') & (res['fold'] >= 0) & (res['metric'] == 'mpsp_mape')]
    if len(lp):
        lpt = lp.pivot_table(index='label', columns='method', values='value', aggfunc='mean')
        fm = res[(res['scheme'] == 'product') & (res['fold'] >= 0) & (res['metric'] == 'fit_mape')]
        lpt['fit_mape'] = fm.groupby('label')['value'].mean()
        lpt = lpt[[m for m in METHOD_ORDER if m in lpt] + ['fit_mape']].reset_index().rename(columns={'label': 'product'})
        checks['lopo_per_product'] = lpt
        lpt.to_csv(os.path.join(tdir, 'cv_lopo_per_product.csv'), index=False)
    test_cfgs = [c['config_id'] for c in sp['test_configs']]
    fp = fin['final'].get('pinn')
    if fp and not fp['error'] and 'test_ood' in fp['evals']:
        checks['test_per_config'] = per_group_mape(fp['evals']['test_ood'], S, Pts, 'config_id')
        checks['test_per_config'].to_csv(os.path.join(tdir, 'final_test_per_config.csv'), index=False)
    checks['pub'] = pub_comparison(D['samples'], set(test_cfgs))
    checks['pub'].to_csv(os.path.join(tdir, 'publication_comparison.csv'), index=False)
    checks['sens'], sens_pts = sensitivity_table(fin['final'], S)
    checks['sens'].to_csv(os.path.join(tdir, 'price_sensitivity.csv'), index=False)
    checks['channels'] = channel_table(fp, S, test_cfgs)
    checks['channels'].to_csv(os.path.join(tdir, 'channel_breakdown.csv'), index=False)
    T('tables', t0)

    # ---- figures
    t0 = time.time()
    made = []
    figure_errors = []

    def _fig(name, fn, *args):
        try:
            fn(os.path.join(fdir, name), *args, made)
        except Exception as e:
            figure_errors.append(f'{name}: {type(e).__name__}: {e}')
            log(f'figure {name} FAILED: {type(e).__name__}: {e}\n{traceback.format_exc()[-1500:]}')

    if fp and not fp['error'] and 'test_all' in fp['evals']:
        ev = fp['evals']['test_all']
        st = S.loc[ev['ids']]
        Z = st[list(TARGETS)].to_numpy(float)
        Pp = np.column_stack([ev['coef'], ev['g']])
        grp = np.where(st['config_id'].isin(test_cfgs), 'held-out configurations', 'in-distribution samples')
        _fig('parity_test.png', fig_parity, Z, Pp, grp,
             {'held-out configurations': PALETTE[0], 'in-distribution samples': PALETTE[1]},
             {'held-out configurations': 'o', 'in-distribution samples': 's'},
             'Final test set: PI-NN ensemble trained on the development set (models/final_dev)')
    oof = cv['oof']
    lopo = oof[(oof['scheme'] == 'product') & (oof['method'] == 'pinn')] if len(oof) else oof
    if len(lopo):
        st = S.loc[lopo['sample_id']]
        Z = st[list(TARGETS)].to_numpy(float)
        Pp = lopo[[f'{k}_pred' for k in TARGETS]].to_numpy(float)
        grp = st['biorefinery'].astype(str).to_numpy()
        _fig('parity_lopo.png', fig_parity, Z, Pp, grp, BIO_COLOR, BIO_MARKER,
             'Leave-one-product-out CV, out-of-fold: PI-NN ensemble (colour = held-out product family)')
    panels = []
    if fp and not fp['error'] and 'test_all' in fp['evals']:
        panels.append(('Final test set (PI-NN, models/final_dev)', _points_pred(fp['evals']['test_all'], S, Pts)))
    if len(lopo):
        ev = dict(ids=lopo['sample_id'].tolist(), coef=lopo[[f'{k}_pred' for k in 'abcd']].to_numpy(float),
                  g=lopo['g_pred'].to_numpy(float))
        panels.append(('Leave-one-product-out CV, out-of-fold (PI-NN)', _points_pred(ev, S, Pts)))
    if panels:
        _fig('mpsp_parity.png', fig_mpsp_parity, panels)
    if len(cv['headline']):
        _fig('cv_mape_by_scheme.png', fig_cv_bars, cv['headline'])
        _fig('ablations.png', fig_ablation, cv['headline'])
    if fp and not fp['error'] and 'test_ood' in fp['evals']:
        ev = fp['evals']['test_ood']
        st = S.loc[ev['ids']].reset_index(drop=True)
        cases = []
        for cid in test_cfgs[:2]:
            idx = np.flatnonzero((st['config_id'] == cid).to_numpy())
            if not len(idx):
                continue
            kinds = st.loc[idx, 'sample_kind'].astype(str).to_numpy()
            i = idx[np.flatnonzero(kinds == 'baseline')[0]] if (kinds == 'baseline').any() else idx[0]
            row = st.iloc[i]
            pts = Pts.get(row['sample_id'])
            if pts is None or len(pts) < 3:
                continue
            cases.append((f'{cid} ({row["sample_kind"]} sample, r = {float(row["r"]):.3g} g/L/h)', row,
                          ev['coef'][i], ev['g'][i], pts))
        if cases:
            _fig('mpsp_contours_test.png', fig_contours, cases)
    if len(checks['pub']):
        _fig('publication_comparison.png', fig_pub, checks['pub'])
    if len(sens_pts):
        _fig('price_sensitivity.png', fig_sens, sens_pts)
    if len(checks['channels']):
        _fig('channels_test.png', fig_channels, checks['channels'])
    _fig('hpo_search.png', fig_hpo, hpo['table'])
    T('figures', t0)

    # ---- timings
    rows = []
    allres = {**cv['raw'], **fin['raw']}
    for j in cv['jobs'] + fin['jobs']:
        r = allres.get(j['key'])
        if r is None:
            continue
        rows.append(dict(stage=j['stage'], method=j['method'], seconds=r.get('seconds', np.nan),
                         n_train=len(j['train_ids']), members=j['members'],
                         member_epochs=float(sum(r.get('epochs') or [])),
                         n_fit=(r.get('info') or {}).get('n_fit', np.nan)))
    tj = pd.DataFrame(rows)
    nn = tj[tj['member_epochs'] > 0]
    sec_per = float((nn['seconds'] / (nn['member_epochs'] * nn['n_fit'].astype(float))).median()) if len(nn) else np.nan
    timings = pd.DataFrame([dict(stage=k, minutes=v / 60) for k, v in stage_t.items()]
                           + [dict(stage='total (this run, incl. cached)', minutes=(time.time() - t_start) / 60)])
    tim_fit = tj.groupby(['stage', 'method']).agg(fits=('seconds', 'size'), mean_seconds=('seconds', 'mean'),
                                                  mean_member_epochs=('member_epochs', 'mean')).reset_index()
    tim_fit.to_csv(os.path.join(tdir, 'timings_fits.csv'), index=False)
    timings['note'] = ''
    timings.loc[len(timings)] = dict(stage='median s per (member-epoch x training sample)', minutes=np.nan,
                                     note=f'{sec_per:.3g}')
    timings.to_csv(os.path.join(tdir, 'timings.csv'), index=False)

    rpath = write_report(a, ctx, D, sp, dev, plan, hpo, cv, fin, checks, made, timings)
    log(f'figures: {len(made)} written{"; errors: " + str(figure_errors) if figure_errors else ""}')
    log(f'report: {rpath}')
    log(f'===== done in {(time.time() - t_start) / 60:.1f} min (stages: '
        + ', '.join(f'{k} {v / 60:.1f}' for k, v in stage_t.items()) + ' min)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
