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
   Every network fit of stages 3-5 uses the same training procedure
   (``PROTOCOL``; v1 behaviour with ``--v1-protocol``): no calibration
   hold-out, every ensemble member early-stops on its own random 15 %
   sample-level hold-out, no refit -- so the CV estimates the training
   procedure of the final models with the selected hyperparameters (the HPO
   is not nested in the CV, so slightly optimistic for HPO + fit);
3. ``hpo``     random search (candidate 0 = package defaults) of width, depth,
               dropout, weight decay, learning rate, w_mpsp, w_ch, w_phys,
               collocation ratio, early-stopping criterion and input clipping by
               K-fold leave-configurations-out CV on the development set (folds
               stratified by biorefinery; ``hpo_members`` networks per fit),
               selecting by the fold mean of the configuration-balanced held-out
               MPSP-surface MAPE (``HPO_METRIC``). Full mode sizes the search
               (and, if needed, the number of LOCO folds of the PI-NN
               ablations) to ``--budget-hours`` from the measured cost of the
               first candidate (``plan.json``);
4. ``cv``      generalization CV on the development set with the selected
               hyperparameters: random-by-sample, leave-one-configuration-out
               (every development configuration for the PI-NN, plain MLP and
               baselines; the ablations on a stratified subset), leave-one-
               feedstock-out, leave-one-product-out and productivity
               extrapolation; methods PI-NN ensemble, plain MLP, ridge, HGB, kNN
               and the PI-NN ablations no_channels / no_price / no_phys (config
               and product schemes). The network fits also return the
               per-sample conformal scores of their held-out samples
               (``tables/cv_oof_scores.csv.gz``) and a cross-fitted coverage
               estimate (``tables/cv_coverage_crossfit.csv``);
5. ``final``   selected PI-NN ensemble (and the plain MLP and the baselines)
               trained on the whole development set and evaluated on the
               final test set (``models/final_dev``), its intervals calibrated
               from the out-of-fold scores of the development LOCO CV
               (``train.calibrate_conformal_oof``; never from test data); then
               the release ensemble trained on ALL data (``models/release``,
               ``train.Ensemble`` checkpoint format), calibrated from the same
               scores. CV, final_dev and release ensembles have the same number
               of members (``--members``);
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
CODE_VERSION = '5'

DEFAULT_RAW = '/home/user/work/data/raw'
DEFAULT_PROCESSED = os.path.join(PKG_DIR, 'data', 'processed')
DEFAULT_RESULTS = os.path.join(PKG_DIR, 'results')
DEFAULT_MODELS = os.path.join(PKG_DIR, 'models')
PUB_DIR = os.path.join(os.path.dirname(PKG_DIR), 'results', 'data')

TARGETS = ('a', 'b', 'c', 'd', 'g')
SCHEMES = ('random', 'config', 'feedstock', 'product', 'productivity')
#: Training protocol of EVERY network fit of the pipeline -- HPO, all CV
#: schemes (PI-NN, plain MLP, ablations), ``models/final_dev`` and
#: ``models/release`` -- so that the CV estimates the training procedure of the
#: final models with the selected hyperparameters (the HPO is not nested in
#: the CV; round 4, DESIGN.md section 6): no calibration hold-out
#: (``calib_frac`` 0); every member early-stops on its own random 15 %
#: SAMPLE-level hold-out (``es_group`` None; the hold-outs differ between
#: members, so the ensemble sees every sample); no refit. Intervals of the
#: final models are calibrated from the out-of-fold conformal scores of the
#: development leave-one-configuration-out CV of this same procedure
#: (``train.calibrate_conformal_oof``), never from test data. On the dev proxy
#: of the frozen test it beat the v1 protocol (configuration-grouped
#: calibration and early-stopping hold-outs): OOD MAPE 7.2 vs 12.0 %, ID 3.9 vs
#: 13.9 % (6 draws, paired).
PROTOCOL = dict(name='v2', calib_frac=0.0, es_group=None, refit='none')
#: the v1 protocol (CLI ``--v1-protocol``; commit fd3f386): CV calibration and
#: early-stopping hold-outs grouped per scheme by ``INNER_GROUP_V1``,
#: ``final_dev`` with a 15 % configuration-grouped calibration hold-out,
#: HPO and release without one, members early-stopping on held-out
#: configurations. Kept to reproduce v1; the release intervals then come
#: from the dev LOCO out-of-fold scores as well (v1 transferred them from the
#: dev model scored on the test configurations, which is no longer done).
PROTOCOL_V1 = dict(name='v1', calib_frac=0.15, es_group='inner', refit='none')
#: v1 inner calibration / early-stopping grouping per scheme. Like the package
#: default (train.SCHEME_GROUP) except leave-one-feedstock-out: with only four
#: feedstocks, grouping the calibration AND validation hold-outs by feedstock
#: would leave a single feedstock to train on, so they hold out whole
#: configurations instead (intervals are then calibrated for new
#: configurations, not for new feedstocks). Round 4 found that grouped
#: early stopping on a few held-out groups is noisy and often stops far too
#: early (dev proxy: best epochs 6-226; LOPO: members stop within ~25 epochs).
INNER_GROUP_V1 = {'random': None, 'config': 'config_id', 'feedstock': 'config_id',
                  'product': 'fermentation_product', 'productivity': 'config_id'}
#: early-stopping grouping per scheme under PROTOCOL (report alias: sample
#: level everywhere)
INNER_GROUP = {s: None for s in SCHEMES}
BASE_METHODS = ('pinn', 'plain', 'ridge', 'hgb', 'knn')
ABLATIONS = ('pinn-no_channels', 'pinn-no_price', 'pinn-no_phys')
ABLATION_SCHEMES = ('config', 'product')
NN_BASES = ('pinn', 'plain')

#: PI-NN extrapolation safeguards of the pipeline (need the downstream
#: FeatureSpec; DESIGN.md sections 3 and 6, round 4): the sugar-availability
#: prior of the structured feedstock intensity, the neutralization gate, and
#: clipping of the feedstock-composition inputs to the training range
EXTRAPOLATION_SAFEGUARDS = dict(fs_prior=True, neutralization_gate=True, clip_groups=['feedstock', 'pretreatment'])
#: hard productivity structure of the pipeline PI-NN (DESIGN.md section 3,
#: round 4; CLI ``--no-r-structure``)
R_STRUCTURE = 'physics'
#: round 5 (DESIGN.md section 6, round 5; each has a CLI flag restoring v2):
#: hard power gate of the pipeline PI-NN (ModelConfig.power_gate: net
#: electricity intensity >= 0 where electricity is not a coproduct; needs the
#: downstream FeatureSpec; CLI ``--no-power-gate``)
POWER_GATE = True
#: pure-b' feedstock intensity (ModelConfig.fs_b_only: kg feedstock per kg
#: fermentation product = 1 / (eta w_sugar y)): NOT adopted (round 5: better
#: new-feedstock sensitivities, but the recovery-shaped feedstock terms the
#: data do have move into other channels, LOPO channel error +4 points);
#: opt-in with CLI ``--fs-b-only``
FS_B_ONLY = False
#: encoder columns never clipped (ModelConfig.clip_exempt; every network of
#: the pipeline): the final-product stoichiometric mass ratio and its
#: MW-ratio twin are exact structural quantities with a tiny training spread;
#: with clip_inputs >= 0 a held-out product's ratio is clipped onto the
#: value of other products (round 5: HP's 0.80 to 1.00, the value of the
#: products without a mass-changing conversion, while its upgrading features
#: say dehydration; HP LOPO 13 -> 82 %); CLI ``--no-clip-exempt``
CLIP_EXEMPT = ('final_product_stoich_mass_ratio', 'phys_pp_fp_MW_ratio')
#: channel loss of the pipeline PI-NN: the channel / intensity term compared
#: as surfaces at the sample's points relative to the total numerator
#: (LossConfig.ch_space = 'surface') with a FIXED weight chosen on the
#: out-of-fold channel / price-intensity fidelity (not by the MPSP HPO
#: criterion, which cannot see the decomposition; w_ch is then no HPO
#: dimension); CLI ``--v2-channel-loss``: coefficient space, w_ch from the HPO
CHANNEL_LOSS = dict(ch_space='surface', w_ch=3.0)
#: channel targets masked in the PI-NN loss (train.mask_channel_targets),
#: by biorefinery: known defective decompositions (DESIGN.md section 1,
#: Known data defects, and section 6, round 5, 5f): TAL's feedstock channel
#: (D8: y is the yield on all substrates, but the sugar demand per kg TAL is
#: 1/y - 1.5 a/t with the spiked acetate a, a titer term no input describes
#: and the product-blind feedstock head cannot represent; c' < 0 in every
#: TAL sample). Lactic's feedstock targets (a recovery shape, not a defect)
#: are NOT masked: on the corrected data the extra mask helped only the corn
#: leave-one-feedstock-out fold, in 2 of 3 seeds, and cost
#: leave-one-product-out 2 points (5f); CLI ``--no-target-mask``
DEFECTIVE_TARGETS = {'feedstock': ('TAL',)}
#: report text (changelog v2 -> v3): the re-check of the round-5 settings on
#: the corrected development rows (DESIGN.md section 6, round 5, 5f)
ROUND5_RECHECK = (
    '**Re-check on the corrected development rows** (data hash c345efda41c4fdca, 1353 samples / 37 '
    'configurations; v2 model settings vs. changes 1-4 of v3, v2-selected hyperparameters, one member per '
    'fit, paired over 3 seeds on the same folds; mean over seeds (difference, seed t); `*` = without TAL and '
    'lactic; `wrong-sign el.` = share of the held-out samples of configurations without power export predicted '
    'to export; %; DESIGN.md 5f):\n\n'
    '| scheme | MPSP fold mean | MPSP geo-mean (configurations) | folds better / worse (median diff) | '
    'channel error ch_rms | feedstock sensitivity error* | electricity response | wrong-sign el. |\n'
    '|---|---|---|---|---|---|---|---|\n'
    '| leave-one-configuration-out (37) | 9.7 -> 9.3 (-0.4, -0.9) | 5.6 -> 5.7 (+0.1, +1.9) | 14 / 23 (+0.12) | '
    '25.2 -> 15.7 (-9.5, -92) | 15.1 -> 4.2 | 9.8 -> 7.7 | 85 -> 0 |\n'
    '| leave-one-product-out (7) | 57.2 -> 46.5 (-10.7, -4.4) | 52.0 -> 28.1 | 6 / 1 (-5.2) | '
    '38.6 -> 31.5 (-7.1, -11.7) | 19.4 -> 18.1 | 14.8 -> 14.0 | 100 -> 0 |\n'
    '| leave-one-feedstock-out (4) | 30.3 -> 34.1 (+3.8, +2.0) | 17.6 -> 16.6 | 2 / 2 (+0.9) | '
    '38.9 -> 35.0 (-3.9, -12.6) | 16.6 -> 12.7 | 16.1 -> 14.4 | 100 -> 0 |\n\n'
    'Per configuration (LOCO): the corn feedstock-sensitivity bias is -32 ... -16 % under v2 for every corn '
    'configuration and -0.6 ... +1.0 % now (isobutanol +8.7 %); isobutanol 70.5 -> 39.7 %, oxalic 52.1 -> '
    '47.3 %, HP in LOPO 54.9 -> 13.2 %. Worse: the three ethanol configurations (LOCO corn 18.0 -> 27.9, corn '
    'stover 27.4 -> 38.0, sugarcane 17.8 -> 20.2 %; LOFO ethanol corn 48 -> 111 %, ethanol corn stover 22 -> '
    '35 %; LOPO ethanol 60 -> 77 %), which carry the LOFO deterioration. Attribution (the same settings without the power '
    'gate): the ethanol corn stover loss comes from the surface channel loss (LOCO 27.4 -> 51.9 % without the '
    'gate, 38.0 % with it), the ethanol corn loss from the gate, which removes its compensating wrong-sign '
    'electricity credit (LOCO 19.3 -> 27.9 %; LOFO corn fold 48.8 -> 63.3 %, all of the LOFO change); the gate '
    'improves the LOCO fold mean (10.2 -> 9.3 %) and LOPO (52.0 -> 46.5 %; isobutanol 72 -> 40 %).\n')
#: channel-wise input masking preset of the pipeline PI-NN ('physics' = v1;
#: CLI ``--v1-g-inputs``): 'physics_v2' keeps the reagent prices away from the
#: recovery g
MASKING_PRESET = 'physics_v2'
#: development configurations of the families with a held-out test
#: configuration (the frozen test holds out one configuration of each):
#: LOCO restricted to them is the matched CV estimate of the OOD test
OOD_FAMILIES = ('HP', 'TAL', 'succinic', 'lactic')
#: one-line captions of the report figures (section 9)
FIGURE_CAPTIONS = {
    'cv_mape_by_scheme.png': 'MPSP MAPE of every method per CV scheme (mean over folds ± s.e.; dashed: '
                             'irreducible fit MAPE).',
    'ablations.png': 'Paired fold differences (variant - full PI-NN fold MAPE, pp; dots = folds, bar = mean ± '
                     '1.96 s.e.) on the folds where the ablations ran (LOCO: the stratified subset of folds).',
    'mpsp_parity.png': 'Predicted vs. simulated MPSP at the simulated TY points: final test (models/final_dev) '
                       'and leave-one-product-out out-of-fold.',
    'parity_test.png': 'Predicted vs. fitted coefficients a′, b′, c′, d′, g on the final test set.',
    'parity_lopo.png': 'Predicted vs. fitted coefficients, leave-one-product-out out-of-fold (colour = held-out '
                       'product family).',
    'mpsp_contours_test.png': 'MPSP over yield and titer for held-out test configurations: simulation vs. '
                              'prediction.',
    'price_sensitivity.png': 'Predicted price derivatives of the coefficients (central secants) vs. the simulated '
                             'price intensities, final dev-trained models on the test set.',
    'channels_test.png': 'Per-channel coefficient breakdown on the held-out test configurations: prediction vs. '
                         'decomposed simulation data.',
    'publication_comparison.png': 'Baseline coefficients of this work (current BioSTEAM, harness v2-v4, '
                                  'D3-patched; D16 / D17 corrected) vs. the publication.',
    'hpo_search.png': 'HPO candidates ranked by the criterion (configuration-balanced held-out MPSP MAPE, fold '
                      'mean); the selected one in blue.',
}
#: samples per configuration planned by data_generation/scripts/run_all.sh
#: (1 baseline + N_PROD productivity + N_RANDOM random samples)
PLANNED_SAMPLES = dict(oxalic=211, isobutanol=109, ethanol=109, HP=29, TAL=29, lactic=36, succinic=18)
#: modules whose edits change the processed tables (their hash is part of the
#: processed-data fingerprint)
PROCESSING_SOURCES = ('dataset.py', 'physics.py', 'data_patches.py', 'downstream.py')
#: nominal miscoverage of the conformal intervals
ALPHA = 0.1

HP_DEFAULT = dict(width=128, depth=3, dropout=0.05, weight_decay=1e-4, lr=2e-3,
                  w_mpsp=1.0, w_ch=0.5, w_phys=0.1, colloc_ratio=1.0, early_stop='data', clip_inputs=-1.0)
#: HPO selection criterion: fold mean of the configuration-balanced held-out
#: MPSP MAPE (every held-out configuration counts once). The point-weighted
#: MAPE of a fold is dominated by its largest configuration (e.g. oxalic: 189
#: of 911 development samples), so selecting on it selects on one
#: extrapolation.
HPO_METRIC = 'mpsp_mape_cfg'
#: weight of leave-one-product-out in the HPO criterion (round 5, DESIGN.md
#: section 6, 5e-5f; CLI ``--hpo-lopo-weight``): criterion = (1 - w) * (LOCO
#: fold mean of HPO_METRIC) + w * (mean over the leave-one-product-out folds
#: of HPO_METRIC). 0 = the v2 criterion, the pipeline default: proposed in
#: round 5 because the LOCO-only criterion of v2 selected clip_inputs = 0,
#: which clipped the held-out HP's final-product stoichiometric ratio 0.80
#: to 1.00 (HP LOPO 13 -> 82 %); CLIP_EXEMPT removes that mechanism for
#: every candidate, and with it w = 0.5 still switches the selection to a
#: candidate whose single-seed LOPO advantage (46.1 vs 49.7 %) reverses over
#: three seeds (54.5 +- 8.0 vs 49.5 %) at +1.3 points of the LOCO criterion,
#: +21 % run time and a post-selection LOPO estimate, so it is not adopted
HPO_LOPO_WEIGHT = 0.0

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

def processing_source_hash():
    """Hash of the modules that build the processed tables
    (:data:`PROCESSING_SOURCES`): editing one rebuilds the tables on the
    next run. (The fit cache is keyed by ``CODE_VERSION`` instead, which must
    be bumped by hand when model.py, train.py, features.py or evaluate.py
    change what a fit computes.)"""
    h = hashlib.sha256()
    for f in PROCESSING_SOURCES:
        with open(os.path.join(PKG_DIR, f), 'rb') as fh:
            h.update(f.encode() + b'\0' + fh.read())
    return h.hexdigest()[:16]


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
    src = processing_source_hash()
    have = all(os.path.exists(os.path.join(proc, f'{n}.csv.gz')) for n in ('samples', 'points', 'configs'))
    old = {}
    if os.path.exists(fp_path):
        with open(fp_path) as fh:
            old = json.load(fh)
    src_changed = old.get('processing_sources') not in (None, src)
    if have and not files and not a.rebuild:
        # a clone without the raw simulation records (not in git): the shipped
        # processed tables are the data (README section 8, "Data availability")
        log(f'data: no raw *.jsonl files under {a.raw}: using the shipped processed tables in {proc}'
            + (' (WARNING: the data-processing sources changed since they were built; they cannot be rebuilt '
               'without the raw records)' if src_changed else ''))
    elif have and (a.no_rebuild or (old.get('raw_fingerprint') == fp and not a.rebuild and not src_changed)):
        log(f'data: cached processed tables in {proc} '
            f'({"--no-rebuild" if a.no_rebuild else "raw files and processing sources unchanged"})')
        if old.get('processing_sources') is None and not a.no_rebuild:
            # tables built before the source hash was recorded: adopted as built
            # by the current sources (run with --rebuild if they are not)
            log('data: NOTE the cached tables predate the processing-source hash '
                f'({", ".join(PROCESSING_SOURCES)}); recorded now -- use --rebuild if these modules changed '
                'since the tables were built')
            _dump_json(dict(old, processing_sources=src), fp_path)
    else:
        if src_changed:
            log('data: data-processing sources changed since the tables were built: rebuilding')
        if not files:
            raise SystemExit(f'no raw *.jsonl files under {a.raw}')
        t0 = time.time()
        log(f'data: building processed tables from {len(files)} raw files ...')
        dataset.build(a.raw, out_dir=proc)
        _dump_json(dict(raw_fingerprint=fp, raw_dirs=a.raw, files=files, processing_sources=src,
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
                recent_raw_files=recent, n_raw_files=len(files) or len(old.get('files') or []))


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
    import pandas as pd
    from fermentation_insights.coeff_predictor.evaluate import fit_method, compute_metrics
    from fermentation_insights.coeff_predictor.train import calibrate_conformal_oof, conformal_scores
    t0 = time.time()
    out = dict(key=job['key'], stage=job['stage'], tag=job['tag'], method=job['method'], error=None,
               evals={}, info={}, epochs=[], best_epochs=[], refit_epochs=None, n_train=len(job['train_ids']))
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter('always')
        try:
            s_tr, p_tr = _subset(job['train_ids'])
            if job.get('target_mask'):
                from fermentation_insights.coeff_predictor.train import mask_channel_targets
                s_tr = mask_channel_targets(s_tr, job['target_mask'])
            model = fit_method(job['method'], s_tr, p_tr, job['members'], job['mk'], job['lk'], job['tk'],
                               group_col=job['group_col'], calib_frac=job['calib_frac'], seed=job['seed'],
                               spec_kw=job.get('spec_kw') or None)
            out['fit_seconds'] = time.time() - t0
            out['info'] = {k: v for k, v in (getattr(model, 'fit_info', {}) or {}).items()}
            ens = getattr(model, 'ens', None)
            if ens is not None:
                out['epochs'] = [int(e) for e in ens.meta.get('epochs', [])]
                out['best_epochs'] = [int(e) for e in ens.meta.get('best_epochs', [])]
                out['refit_epochs'] = ens.meta.get('refit_epochs')
                out['n_fit_rows'] = int(ens.meta.get('n_fit', 0))
                # intervals from out-of-fold scores of the same procedure (before the
                # evaluations and the save); reporting-only calibrations go to ens.meta
                own = dict(ens.conformal or {})
                for name, (rec, info) in (job.get('calib_extra') or {}).items():
                    ens.meta[f'conformal_{name}'] = calibrate_conformal_oof(
                        ens, pd.DataFrame(rec), alpha=ALPHA, group_col='config_id', fold_col='fold', **info)
                ens.conformal = own
                if job.get('calib_scores') is not None:
                    calibrate_conformal_oof(ens, pd.DataFrame(job['calib_scores']), alpha=ALPHA,
                                            group_col='config_id', fold_col='fold', **(job.get('calib_info') or {}))
                out['conformal'] = dict(ens.conformal)
                out['conformal_extra'] = {k[len('conformal_'):]: v for k, v in ens.meta.items()
                                          if k.startswith('conformal_') and isinstance(v, dict)}
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
                if 'scores' in job['extras'] and ens is not None:
                    ev['scores'] = conformal_scores(ens, s_te, p_te).to_dict('list')
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
                                protocol=job.get('protocol'), calib_info=job.get('calib_info'),
                                best_epochs=out['best_epochs'],
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
             extras=(), save_model=None, meta=None, spec_kw=None, calib_scores=None, calib_info=None,
             calib_extra=None):
    """One fit + evaluation. ``spec_kw``: FeatureSpec options of every method
    (e.g. ``{'downstream': True}``; default ``ctx['spec_kw']``, off with CLI
    ``--no-downstream``; part of the cache key only when non-empty).
    ``ctx['mk_extra']`` (water gate, off with CLI ``--no-water-gate``; channel-wise
    input masking, off with ``--no-masking``; EXTRAPOLATION_SAFEGUARDS, off with
    ``--no-extrapolation-safeguards``) is merged into ``mk``. A masked
    PI-NN without explicit sub-network sizes gets ``masked_width = width // 2``
    and ``masked_depth = depth - 1`` (64 / 2 for the defaults, as in the
    development CV), so that the HPO dimensions width and depth still apply.

    ``calib_scores`` (records of a ``train.conformal_scores`` table with a
    ``fold`` column, e.g. the dev LOCO out-of-fold scores) calibrates the
    fitted ensemble's intervals with ``train.calibrate_conformal_oof``
    (``calib_info``: its metadata, e.g. source / scheme / weighting) right
    after the fit; ``calib_extra`` = {name: (records, info)} adds
    reporting-only calibrations as ``ens.meta['conformal_<name>']``. Both
    are part of the cache key. ``extras`` 'scores' stores the per-sample
    conformal scores of every evaluation set (``ev['scores']``)."""
    if spec_kw is None:
        spec_kw = ctx.get('spec_kw')
    if ctx.get('mk_extra'):
        mk = {**mk, **ctx['mk_extra']}
    pinn = method.split('-')[0] == 'pinn'
    if ctx.get('lk_extra') and pinn:
        # pipeline channel loss (CHANNEL_LOSS); ablations override it (evaluate._configs)
        lk = {**lk, **ctx['lk_extra']}
    target_mask = {k: list(v) for k, v in (ctx.get('target_mask') or {}).items()} if pinn else {}
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
    if target_mask:
        spec['target_mask'] = target_mask
    if calib_scores is not None or calib_extra:
        spec['calib'] = _sha([calib_scores, calib_info, calib_extra], 24)
    return dict(key=_sha(spec, 24), stage=stage, tag=tag, method=method, train_ids=list(train_ids),
                evals={k: list(v) for k, v in evals.items()}, members=members, mk=mk, lk=lk, tk=tk,
                spec_kw=dict(spec_kw or {}), target_mask=target_mask, calib_scores=calib_scores,
                calib_info=dict(calib_info or {}),
                calib_extra=calib_extra or None, protocol=(ctx.get('protocol') or PROTOCOL)['name'],
                group_col=group_col, calib_frac=calib_frac, seed=seed, extras=tuple(extras),
                save_model=save_model, meta=meta or {}, data_hash=ctx['data_hash'],
                cost=(members if method.split('-')[0] in NN_BASES else 0.3) * len(train_ids))


class Runner:
    """Runs jobs in-process (``workers == 1``) or in a pool of spawned
    single-threaded worker processes; caches each finished job."""

    def __init__(self, workers, processed_dir, dhash, cache_dir, force_pool=False):
        self.workers, self.processed, self.dhash, self.cache = int(workers), processed_dir, dhash, cache_dir
        self.force_pool = bool(force_pool)
        self.n_run = 0          # fits computed (not cached) in this run
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
        self.n_run += len(todo)
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
    """Plan of a run. The ensembles of the CV, final_dev and release fits
    have the same number of members: the conformal scores divide by the
    member spread, so out-of-fold scores of the CV ensembles calibrate the
    final ensembles only at equal M (``check_plan``). LOCO covers every
    development configuration for the PI-NN and the plain MLP; the ablations
    run on a stratified subset of ``ablation_loco_folds`` configurations
    (None = all)."""
    # --ablation-loco-folds: None = the mode's default (full 12, quick 4), 0 = every LOCO fold
    abl = getattr(a, 'ablation_loco_folds', None)
    if a.quick:
        return dict(mode='quick', hpo_k=3, hpo_candidates=3, hpo_members=1, cv_members=2, final_members=2,
                    release_members=2, random_k=3, loco_folds=None,
                    ablation_loco_folds=(4 if abl is None else abl or None), epochs=12, patience=6,
                    min_epochs=4, batch_size=128)
    abl = ABLATION_LOCO_FOLDS if abl is None else abl
    # hpo_members: single-member fits differ between seeds by up to ~50 % (relative)
    # in held-out MAPE on the development folds, so one member per candidate
    # and fold mostly ranks seeds, not hyperparameters
    return dict(mode='full', hpo_k=5, hpo_candidates=None, hpo_members=a.hpo_members, cv_members=a.members,
                hpo_lopo_weight=hpo_lopo_weight(a),
                final_members=a.members, release_members=a.members, random_k=5, loco_folds=None,
                ablation_loco_folds=(abl or None), epochs=300, patience=40, min_epochs=20,
                batch_size=128, request=plan_request(a))


#: default number of stratified LOCO folds of the PI-NN ablations (full mode)
ABLATION_LOCO_FOLDS = 12


def plan_request(a):
    """The CLI settings a full-mode plan is sized from; plan.json is reused
    only when they are unchanged (``stage_hpo``)."""
    abl = getattr(a, 'ablation_loco_folds', None)
    return dict(budget_hours=float(a.budget_hours), hpo_candidates=int(a.hpo_candidates or 0),
                hpo_members=int(a.hpo_members), hpo_lopo_weight=hpo_lopo_weight(a),
                ablation_loco_folds=int(ABLATION_LOCO_FOLDS if abl is None else abl))


def check_plan(plan):
    """Raise unless the CV, final_dev and release ensembles have the same
    number of members (out-of-fold calibration, ``base_plan``)."""
    m = (plan['cv_members'], plan['final_members'], plan['release_members'])
    if len(set(m)) != 1:
        raise SystemExit(f'plan: cv_members, final_members, release_members = {m}; they must be equal '
                         '(the final intervals are calibrated from the CV ensembles\' out-of-fold scores)')
    return plan


def protocol_holdouts(protocol, stage, scheme=None):
    """``(group_col, calib_frac)`` of a network fit of ``stage`` ('hpo',
    'cv', 'final', 'release') under ``protocol`` (PROTOCOL or PROTOCOL_V1).
    The sklearn baselines get the same values, so they fit exactly the rows
    the network members train on."""
    if protocol['name'] == 'v1':
        return dict(hpo=('config_id', 0.0), cv=(INNER_GROUP_V1.get(scheme), 0.15), final=('config_id', 0.15),
                    release=('config_id', 0.0))[stage]
    return protocol['es_group'], float(protocol['calib_frac'])


def hpo_candidates(n, seed, fixed=None):
    """Fixed seeded sequence (candidate 0 = package defaults); the first n.
    ``fixed`` = {dimension: value} overrides a dimension the pipeline sets
    itself (e.g. ``w_ch`` under CHANNEL_LOSS); the random draws are the same,
    so the other dimensions of the candidates are unchanged, and candidates
    that become identical are dropped."""
    import numpy as np
    rng = np.random.RandomState(seed + 2027)
    fixed = dict(fixed or {})
    out = [{**HP_DEFAULT, **fixed}]
    sig = lambda v: float(f'{v:.2g}')
    while len(out) < max(n, 1):
        c = dict(width=int(rng.choice([64, 128, 256])), depth=int(rng.choice([2, 3, 4])),
                 dropout=float(rng.choice([0., 0.05, 0.1, 0.2])),
                 weight_decay=sig(10 ** rng.uniform(-5, -3)), lr=sig(10 ** rng.uniform(np.log10(5e-4), np.log10(5e-3))),
                 w_mpsp=float(rng.choice([0.3, 1., 3.])), w_ch=float(rng.choice([0.1, 0.3, 0.5, 1.])),
                 w_phys=float(rng.choice([0.03, 0.1, 0.3, 1.])), colloc_ratio=float(rng.choice([0.5, 1., 2.])),
                 early_stop=str(rng.choice(['data', 'mpsp'])), clip_inputs=float(rng.choice([-1., 0.])))
        c.update(fixed)
        if c not in out:
            out.append(c)
    return out[:n]


def hpo_fixed(ctx):
    """HPO dimensions fixed by the pipeline settings of ``ctx``: ``w_ch``
    under the channel loss CHANNEL_LOSS (``ctx['lk_extra']``)."""
    lk = (ctx or {}).get('lk_extra') or {}
    return {'w_ch': float(lk['w_ch'])} if 'w_ch' in lk else {}


def hp_kwargs(hp, plan, protocol=None):
    """(ModelConfig, LossConfig, TrainConfig) kwargs of hyperparameters
    ``hp``; ``TrainConfig.refit`` from ``protocol`` (default PROTOCOL)."""
    mk = dict(width=int(hp['width']), depth=int(hp['depth']), dropout=float(hp['dropout']),
              clip_inputs=float(hp.get('clip_inputs', -1.0)))
    lk = dict(w_mpsp=float(hp['w_mpsp']), w_ch=float(hp['w_ch']), w_phys=float(hp['w_phys']),
              colloc_ratio=float(hp['colloc_ratio']))
    tk = dict(lr=float(hp['lr']), weight_decay=float(hp['weight_decay']), epochs=plan['epochs'],
              patience=plan['patience'], min_epochs=plan['min_epochs'], batch_size=plan['batch_size'], threads=1,
              early_stop=str(hp.get('early_stop', 'data')), refit=str((protocol or PROTOCOL)['refit']))
    return mk, lk, tk


def ablation_fold_labels(dev, plan, seed):
    """Labels (configurations) of the LOCO folds the ablations run on: a
    stratified subset of ``plan['ablation_loco_folds']`` configurations, or
    None = every LOCO fold."""
    n = plan.get('ablation_loco_folds')
    if not n or n >= dev['config_id'].nunique():
        return None
    return set(stratified_subset(dev, int(n), seed))


def cv_methods(scheme, label, abl_labels):
    """Methods fit on one CV fold: the base methods everywhere, the PI-NN
    ablations on the ABLATION_SCHEMES (LOCO: only the folds in
    ``abl_labels``, None = all)."""
    meths = list(BASE_METHODS)
    if scheme in ABLATION_SCHEMES and (scheme != 'config' or abl_labels is None or label in abl_labels):
        meths += list(ABLATIONS)
    return meths


def plan_budget(plan, a, pilot, dev, test_n, n_all, seed):
    """Full mode: size the HPO search (and, if needed, the number of LOCO
    folds of the ABLATIONS; the PI-NN and the plain MLP always run on every
    LOCO fold) so that the estimated wall time fits ``--budget-hours`` with
    ``--threads`` workers. Cost model: seconds per member fit proportional to
    the number of training samples, measured on the HPO pilot (candidate 0,
    trained with the pipeline protocol, so the longer sample-level
    early-stopping runs are included)."""
    per = pilot['seconds_per_member_sample']
    hgb = pilot.get('hgb_seconds_per_sample', 0.004)

    def fit(n, members, method='pinn'):
        return members * per * n if method.split('-')[0] in NN_BASES else hgb * n

    def cv_cost(abl_folds):
        pl = dict(plan, ablation_loco_folds=abl_folds)
        abl = ablation_fold_labels(dev, pl, seed)
        tot = 0.
        for sch in SCHEMES:
            for _, tr, _, lab in scheme_splits(dev, sch, plan, seed):
                tot += sum(fit(len(tr), plan['cv_members'], m) for m in cv_methods(sch, lab, abl)
                           if m not in ('ridge', 'knn'))
        return tot

    n_dev = len(dev)
    # final_dev PI-NN + plain MLP, release PI-NN, three baselines
    final = (2 * fit(n_dev, plan['final_members']) + fit(n_all, plan['release_members'])
             + 3 * fit(n_dev, 1, 'hgb'))
    hpo_one = plan['hpo_k'] * fit(n_dev * (plan['hpo_k'] - 1) / plan['hpo_k'], plan.get('hpo_members', 1))
    if plan.get('hpo_lopo_weight', 0.) > 0:      # the leave-one-product-out part of the criterion
        hpo_one += sum(fit(len(tr), plan.get('hpo_members', 1)) for _, tr, _, _ in
                       scheme_splits(dev, 'product', plan, seed))
    budget = a.budget_hours * 3600 * max(1, a.threads)
    n_cand = a.hpo_candidates or 12
    abl = plan.get('ablation_loco_folds')
    total = lambda nc, af: nc * hpo_one + cv_cost(af) + final
    if not a.hpo_candidates:
        # shrink the ablations' LOCO folds (min 6), then candidates (min 8); grow candidates (max 16) with slack
        n_loco = dev['config_id'].nunique()
        while total(n_cand, abl) > budget and (abl or n_loco) > 6:
            abl = max(6, (abl or n_loco) - 2)
        while total(n_cand, abl) > budget and n_cand > 8:
            n_cand -= 1
        while n_cand < 16 and total(n_cand + 1, abl) <= budget:
            n_cand += 1
    est = total(n_cand, abl)
    plan.update(hpo_candidates=int(n_cand), ablation_loco_folds=abl,
                estimate=dict(seconds_per_member_sample=per, hpo_hours=n_cand * hpo_one / 3600 / max(1, a.threads),
                              cv_hours=cv_cost(abl) / 3600 / max(1, a.threads),
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

def hpo_lopo_weight(a):
    """Weight of the leave-one-product-out folds in the HPO criterion
    (``--hpo-lopo-weight``, default HPO_LOPO_WEIGHT; quick mode: 0)."""
    w = getattr(a, 'hpo_lopo_weight', None)
    w = (0. if getattr(a, 'quick', False) else HPO_LOPO_WEIGHT) if w is None else float(w)
    if not 0. <= w <= 1.:
        raise SystemExit(f'--hpo-lopo-weight must be in [0, 1], got {w}')
    return w


def stage_hpo(a, plan, runner, dev, ctx, n_test, n_all):
    folds = stratified_group_kfold(dev, plan['hpo_k'], a.seed)
    w_lopo = float(plan.get('hpo_lopo_weight', 0.))
    lopo = [(tr, te, lab) for _, tr, te, lab in scheme_splits(dev, 'product', plan, a.seed)] if w_lopo > 0 else []
    ids = dev['sample_id'].astype(str).to_numpy()
    n_max = 16 if plan['hpo_candidates'] is None else plan['hpo_candidates']
    cands = hpo_candidates(max(n_max, 1), a.seed, fixed=hpo_fixed(ctx))

    protocol = ctx.get('protocol') or PROTOCOL
    gcol, cfrac = protocol_holdouts(protocol, 'hpo')

    def jobs_for(ci):
        # the pipeline protocol (as in CV and the final models): the HPO selects for the procedure used
        mk, lk, tk = hp_kwargs(cands[ci], plan, protocol)
        js = [make_job('hpo', f'cand{ci} {lab}', 'pinn', ids[tr], {'test': ids[te]}, plan.get('hpo_members', 1),
                       mk, lk, tk, gcol, cfrac, a.seed + 101 * ci + f, ctx, meta=dict(cand=ci, fold=f))
              for f, (tr, te, lab) in enumerate(folds)]
        # leave-one-product-out part of the criterion (HPO_LOPO_WEIGHT)
        js += [make_job('hpo', f'cand{ci} lopo {lab}', 'pinn', ids[tr], {'test': ids[te]},
                        plan.get('hpo_members', 1), mk, lk, tk, gcol, cfrac, a.seed + 101 * ci + 50 + f, ctx,
                        meta=dict(cand=ci, fold=f, part='lopo', label=lab))
               for f, (tr, te, lab) in enumerate(lopo)]
        return js

    results = runner.run(jobs_for(0), 'hpo pilot (candidate 0 = defaults)')
    if plan['hpo_candidates'] is None:      # full mode: size the search from the pilot
        pl_path = os.path.join(a.out, 'plan.json')
        old = {}
        if os.path.exists(pl_path) and not a.replan:
            with open(pl_path) as fh:
                old = json.load(fh)
        same = dict(data_hash=old.get('data_hash') == ctx['data_hash'],
                    split_hash=old.get('split_hash') == ctx['split_hash'], mode=old.get('mode') == 'full',
                    threads=old.get('threads') == a.threads, code_version=old.get('code_version') == CODE_VERSION,
                    protocol=old.get('protocol') == protocol['name'],
                    cv_members=old.get('cv_members') == plan['cv_members'],
                    request=old.get('request') == plan.get('request'))
        if old and all(same.values()):
            plan.update({k: old[k] for k in ('hpo_candidates', 'ablation_loco_folds', 'estimate') if k in old})
            log(f'plan: reusing {pl_path}')
        else:
            if old:
                log(f'plan: not reusing {pl_path} (changed: {", ".join(k for k, v in same.items() if not v)}); '
                    'recomputing it from the pilot')
            ok = [r for r in results.values() if not r['error'] and r.get('n_fit_rows')]
            per = (sum(r['fit_seconds'] for r in ok) / sum(r['n_fit_rows'] * len(r['epochs'] or [1]) for r in ok)
                   if ok else 0.05)
            plan_budget(plan, a, dict(seconds_per_member_sample=per), dev, n_test, n_all, a.seed)
        e = plan['estimate']
        log(f"plan: {plan['hpo_candidates']} HPO candidates, LOCO folds = {plan['loco_folds'] or 'all'} "
            f"(ablations: {plan.get('ablation_loco_folds') or 'all'}); "
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
        js = [j for j in jobs if j['meta']['cand'] == ci and j['meta'].get('part', 'loco') == 'loco']
        jl = [j for j in jobs if j['meta']['cand'] == ci and j['meta'].get('part') == 'lopo']
        rs = [results[j['key']] for j in js]
        ml = []
        for j in jl:
            r = results[j['key']]
            m = {} if r['error'] else {(x, t): v for x, t, v in r['evals']['test']['rows']}
            ml.append(m)
            fold_rows.append(dict(candidate=ci, fold=f"lopo {j['meta']['label']}", error=r['error'],
                                  mpsp_mape=m.get(('mpsp_mape', 'MPSP'), np.nan),
                                  mpsp_mape_cfg=m.get(('mpsp_mape_cfg', 'MPSP'), np.nan),
                                  mpsp_mape_cfg_gmean=m.get(('mpsp_mape_cfg_gmean', 'MPSP'), np.nan),
                                  fit_mape=m.get(('fit_mape', 'MPSP'), np.nan), seconds=r['seconds'],
                                  epochs=np.mean(r['epochs']) if r['epochs'] else np.nan))
        crit_l = np.array([m.get((HPO_METRIC, 'MPSP'), np.nan) for m in ml])
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
        c_loco = float(np.nanmean(crit)) if n_ok == len(crit) else np.inf
        c_lopo = (float(np.mean(crit_l)) if np.isfinite(crit_l).all() else np.inf) if jl else np.nan
        rows.append(dict(candidate=ci, **hp, n_folds_ok=n_ok,
                         criterion=((1 - w_lopo) * c_loco + w_lopo * c_lopo) if jl else c_loco,
                         criterion_loco=c_loco, criterion_lopo=c_lopo,
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
    log(f'hpo: selected candidate {best}: {sel} (criterion {tab.loc[best, "criterion"]:.4f} = '
        f'{1 - w_lopo:g} x LOCO {tab.loc[best, "criterion_loco"]:.4f} + {w_lopo:g} x LOPO '
        f'{tab.loc[best, "criterion_lopo"]:.4f}; point-weighted LOCO {tab.loc[best, "mpsp_mape"]:.4f})')
    return dict(table=tab, folds=pd.DataFrame(fold_rows), selected=sel, best=best, fold_info=fold_info,
                jobs=jobs, raw=results, lopo_weight=w_lopo)


def stage_cv(a, plan, runner, dev, ctx, hp, D):
    """Generalization CV of every method with the pipeline protocol (the
    same procedure as the final models). The PI-NN and plain-MLP fits also
    return the per-sample conformal scores of their held-out samples:
    ``cv['oof_scores']`` (scheme, method, fold = fold label, fold_index and
    the ``train.conformal_scores`` columns; every development sample once
    per scheme and method, except in the productivity scheme, which holds
    out and scores only its two tails), from which the final models are
    calibrated (LOCO) and
    ``cv['crossfit']`` estimates the coverage."""
    import pandas as pd
    protocol = ctx.get('protocol') or PROTOCOL
    mk, lk, tk = hp_kwargs(hp, plan, protocol)
    ids = dev['sample_id'].astype(str).to_numpy()
    abl = ablation_fold_labels(dev, plan, a.seed)
    jobs, meta = [], {}
    for sch in SCHEMES:
        gcol, cfrac = protocol_holdouts(protocol, 'cv', sch)
        for f, tr, te, lab in scheme_splits(dev, sch, plan, a.seed):
            for m in cv_methods(sch, lab, abl):
                ex = ('scores', 'channels') if m == 'pinn' else ('scores',) if m in NN_BASES else ()
                j = make_job('cv', f'{sch} {lab}', m, ids[tr], {'test': ids[te]}, plan['cv_members'], mk, lk, tk,
                             gcol, cfrac, a.seed + f, ctx, extras=ex, meta=dict(scheme=sch, fold=f, label=lab))
                jobs.append(j)
    results = runner.run(jobs, 'cv')
    S = D['samples'].assign(sample_id=D['samples']['sample_id'].astype(str)).set_index('sample_id', drop=False)
    P = {k: v for k, v in D['points'].assign(sample_id=D['points']['sample_id'].astype(str)).groupby('sample_id')}
    rows, oofs, scores = [], [], []
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
            if not r['error'] and r['evals']['test'].get('scores') is not None:
                scores.append(pd.DataFrame(r['evals']['test']['scores']).assign(
                    scheme=sch, method=j['method'], fold=str(j['meta']['label']), fold_index=j['meta']['fold']))
        if failed:
            log(f'cv {sch}: pooled metrics exclude fold(s) {sorted(failed)} (a method failed there)')
        if len({j['meta']['fold'] for j in js}) > 1:
            pr, oof = pooled_rows(sch, parts, S, P)
            rows += pr
            oofs += oof
    res = pd.DataFrame(rows)
    oof = pd.concat(oofs, ignore_index=True) if oofs else pd.DataFrame()
    sc = pd.concat(scores, ignore_index=True) if scores else pd.DataFrame()
    if len(sc):
        lead = ['scheme', 'method', 'fold', 'fold_index']
        sc = sc[lead + [c for c in sc.columns if c not in lead]]
    return dict(results=res, oof=oof, jobs=jobs, raw=results, oof_scores=sc, crossfit=crossfit_table(sc),
                ablation_folds=sorted(abl) if abl is not None else None,
                fidelity=oof_fidelity(jobs, results, S, P))


def oof_fidelity(jobs, results, samples_by_id, points_by_id, method='pinn'):
    """Out-of-fold channel / price-intensity fidelity of ``method``
    (``evaluate.channel_fidelity`` of each held-out sample, from the fold
    ensemble that did not see it; member-mean channels): one row per scheme
    x held-out sample, with ``fold`` (index) and ``label``."""
    import pandas as pd
    from fermentation_insights.coeff_predictor.evaluate import channel_fidelity
    out = []
    for j in jobs:
        if j['method'] != method:
            continue
        r = results.get(j['key']) or {}
        ev = (r.get('evals') or {}).get('test')
        if r.get('error') or not ev or ev.get('channels') is None or ev.get('intensities') is None:
            continue
        ids = [str(i) for i in ev['ids']]
        s = samples_by_id.loc[ids].reset_index(drop=True)
        parts = [points_by_id[i] for i in ids if i in points_by_id]
        if not parts:
            continue
        f = channel_fidelity(s, pd.concat(parts, ignore_index=True), ev['channels'], ev['channel_names'],
                             ev['intensities'], ev['price_keys'], ev['g'])
        out.append(f.assign(scheme=j['meta']['scheme'], fold=j['meta']['fold'], label=str(j['meta']['label'])))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


#: columns of the report's fidelity summary
FIDELITY_REPORT_COLS = ('ch_rms', 'ch_feedstock', 'ch_electricity', 'ch_natural_gas', 'ch_materials_other',
                        'ch_capex_downstream', 'sens_fs', 'bias_fs', 'resp_fs', 'resp_el', 'resp_ng')


#: per-configuration relative price-sensitivity errors (``sens_*``, ``bias_*``)
#: are reported only where the simulated price share of MPSP is at least this
FIDELITY_MIN_SHARE = 0.02


def fidelity_tables(fid):
    """(summary per scheme, per scheme x configuration) of the out-of-fold
    fidelity rows: configuration-balanced means (``bias_*`` absolute in the
    summary, signed medians per configuration); feedstock metrics also
    without the configurations whose feedstock targets are masked as
    defective (``*_valid``, DEFECTIVE_TARGETS)."""
    import numpy as np
    import pandas as pd
    from fermentation_insights.coeff_predictor.evaluate import fidelity_summary
    if fid is None or not len(fid):
        return pd.DataFrame(), pd.DataFrame()
    d = fid.copy()
    bad = d['biorefinery'].isin(DEFECTIVE_TARGETS.get('feedstock', ())) | \
        d['config_id'].isin(DEFECTIVE_TARGETS.get('feedstock', ()))
    for c in ('ch_feedstock', 'sens_fs', 'bias_fs', 'resp_fs'):
        if c in d:
            d[c + '_valid'] = np.where(bad, np.nan, d[c])
    summ = fidelity_summary(d, by='scheme')
    cols = [c for c in d.columns if c.startswith(('ch_', 'sens_', 'bias_', 'resp_', 'share_'))
            and not c.endswith('_valid')]
    per = d.groupby(['scheme', 'config_id'])[cols].median().reset_index()
    per.insert(2, 'n_samples', d.groupby(['scheme', 'config_id']).size().to_numpy())
    # relative errors of a price key are undefined where its simulated share
    # of MPSP is about 0 (DESIGN.md 5a): NaN below FIDELITY_MIN_SHARE (the
    # absolute resp_* columns remain)
    for k in ('fs', 'el', 'ng'):
        if f'share_{k}' in per:
            low = per[f'share_{k}'].abs() < FIDELITY_MIN_SHARE
            for c in (f'sens_{k}', f'bias_{k}'):
                if c in per:
                    per.loc[low, c] = np.nan
    return summ, per


def crossfit_table(scores, alpha=ALPHA):
    """Cross-fitted coverage estimate (``train.crossfit_coverage``, leave
    one CV fold out) per scheme x method x interval key of the out-of-fold
    score table, with both weightings."""
    import pandas as pd
    from fermentation_insights.coeff_predictor.train import crossfit_coverage, CONF_KEYS
    rows = []
    if not len(scores):
        return pd.DataFrame(columns=['scheme', 'method', 'weighting', 'key', 'mean', 'pooled', 'min', 'q10',
                                     'frac_ok', 'n_folds'])
    for (sch, meth), g in scores.groupby(['scheme', 'method'], sort=False):
        if g['fold'].nunique() < 2:
            continue
        for w in ('sample', 'group'):
            cf = crossfit_coverage(g, alpha, 'fold', CONF_KEYS, weighting=w, group_col='config_id')
            rows += [dict(scheme=sch, method=meth, weighting=w, key=k, **v) for k, v in cf.items()]
    return pd.DataFrame(rows)


def oof_calibration_scores(cv, scheme, method, dev, sp, strict=True, drop=False):
    """Records of the out-of-fold conformal scores of ``method`` in CV
    ``scheme``, checked for use as a calibration set: every development
    sample scored exactly once, none of them a test sample, and >=
    train.MIN_CALIB_GROUPS configurations. Otherwise: raises (``strict``;
    a test sample in the table raises always), or logs a warning and
    returns None (``drop``: the model gets no intervals) or the incomplete
    records (quick mode)."""
    from fermentation_insights.coeff_predictor.train import MIN_CALIB_GROUPS, EPS_KEYS
    sc = cv.get('oof_scores')
    if sc is None or not len(sc) or not ((sc['scheme'] == scheme) & (sc['method'] == method)).any():
        msg = f'no out-of-fold conformal scores of {method} in the {scheme} CV'
        if strict:
            raise SystemExit('final: ' + msg)
        log('final: WARNING ' + msg + ('; no intervals for it' if drop else ''))
        return None
    sc = sc[(sc['scheme'] == scheme) & (sc['method'] == method)].reset_index(drop=True)
    ids = sc['sample_id'].astype(str)
    dev_ids = set(dev['sample_id'].astype(str))
    test = set(map(str, sp['ood_test_sample_ids'])) | set(map(str, sp['id_test_sample_ids']))
    problems = []
    if ids.duplicated().any():
        problems.append(f'{int(ids.duplicated().sum())} samples scored more than once')
    if set(ids) != dev_ids:
        problems.append(f'{len(dev_ids - set(ids))} development samples not scored, '
                        f'{len(set(ids) - dev_ids)} scored samples outside the development set')
    if set(ids) & test:
        problems.append(f'{len(set(ids) & test)} scored samples are test samples')
    n_cfg = int(sc['config_id'].astype(str).nunique()) if len(sc) else 0
    if n_cfg < MIN_CALIB_GROUPS:
        problems.append(f'{n_cfg} configurations (< {MIN_CALIB_GROUPS})')
    if problems:
        msg = f'out-of-fold scores of {method} ({scheme} CV) unusable for calibration: ' + '; '.join(problems)
        if set(ids) & test or strict:
            raise SystemExit('final: ' + msg)
        if drop:
            log(f'final: WARNING {msg}; {method} gets no intervals')
            return None
        log('final: WARNING (quick mode) ' + msg)
    keep = ['sample_id', 'config_id', 'fold', 'invalid', 'a', 'b', 'c', 'd', 'g', 'mpsp']
    # raw |error| and spread: the final ensemble re-scores with its own eps (train.rescore_conformal)
    keep += [c for k in EPS_KEYS for c in (f'{k}_abserr', f'{k}_sd') if c in sc]
    return sc[keep].to_dict('list')


def stage_final(a, plan, runner, dev, ctx, hp, sp, D, models_dir, cv):
    """``models/final_dev`` (PI-NN, plain MLP and the baselines trained on the
    development set, evaluated once on the frozen test split) and
    ``models/release`` (PI-NN on ALL samples), all with the pipeline protocol.
    Under PROTOCOL the PI-NN / plain intervals are calibrated from their own
    out-of-fold scores of the development LOCO CV (same procedure;
    ``train.calibrate_conformal_oof``); the release PI-NN uses the SAME dev
    LOCO PI-NN scores. Test samples are never used to calibrate or select."""
    if sp['test_configs_missing']:
        msg = (f"final: held-out test configuration(s) {sp['test_configs_missing']} have no samples; the frozen "
               f"split ({a.split_file}) is never changed -- generate their data before the final evaluation")
        if not a.quick:
            raise SystemExit(msg)
        log('WARNING (quick mode only) ' + msg)
    check_plan(plan)
    protocol = ctx.get('protocol') or PROTOCOL
    mk, lk, tk = hp_kwargs(hp, plan, protocol)
    s = D['samples']
    dev_ids = dev['sample_id'].astype(str).tolist()
    evals = dict(test_ood=list(sp['ood_test_sample_ids']), test_id=list(sp['id_test_sample_ids']),
                 test_all=sorted(set(sp['ood_test_sample_ids']) | set(sp['id_test_sample_ids'])))
    all_ids = s['sample_id'].astype(str).tolist()
    # a missing or incomplete PI-NN score table stops the run (full mode); the
    # plain MLP's only costs the plain final_dev its intervals
    cal = {m: oof_calibration_scores(cv, 'config', m, dev, sp, strict=not a.quick and m == 'pinn',
                                     drop=not a.quick and m != 'pinn') for m in NN_BASES}
    rnd = {m: oof_calibration_scores(cv, 'random', m, dev, sp, False, drop=not a.quick) for m in NN_BASES}
    src = ('out-of-fold scores of the development leave-one-configuration-out CV '
           f'({plan["cv_members"]}-member ensembles, same training procedure)')

    def info(m, release=False):
        return dict(source=src + ('; the release ensemble trains on development + test samples' if release else ''),
                    scheme='config', cv_members=plan['cv_members'], weighting='sample', method_scores=m)

    def extra(m):
        if rnd[m] is None:
            return None
        return {'random_cv': (rnd[m], dict(source='out-of-fold scores of the development random CV (reporting '
                                                  'only: in-distribution diagnostics)', scheme='random',
                                           cv_members=plan['cv_members'], weighting='sample'))}

    gcol, cfrac = protocol_holdouts(protocol, 'final')
    # v1 protocol: final_dev keeps its own calibration hold-out (as in v1)
    own_hold_out = cfrac > 0
    jobs = []
    for m in BASE_METHODS:
        nn = m in NN_BASES
        jobs.append(make_job('final', 'final(dev->test)', m, dev_ids, evals, plan['final_members'], mk, lk, tk,
                             gcol, cfrac, a.seed, ctx,
                             extras=(('channels', 'sens', 'scores') if m == 'pinn' else
                                     ('sens', 'scores') if nn else ('sens',)),
                             save_model=os.path.join(models_dir, 'final_dev') if m == 'pinn' else None,
                             calib_scores=cal[m] if nn and not own_hold_out else None,
                             calib_info=info(m) if nn and not own_hold_out and cal[m] is not None else None,
                             calib_extra=extra(m) if nn else None,
                             meta=dict(kind='final', nn=nn)))
    # release: ALL samples are fitted, intervals from the dev LOCO PI-NN scores
    gcol_r, cfrac_r = protocol_holdouts(protocol, 'release')
    jobs.append(make_job('final', 'release(all data)', 'pinn', all_ids, dict(train_all=all_ids),
                         plan['release_members'], mk, lk, tk, gcol_r, cfrac_r, a.seed, ctx,
                         save_model=os.path.join(models_dir, 'release'), calib_scores=cal['pinn'],
                         calib_info=info('pinn', release=True), calib_extra=extra('pinn'),
                         meta=dict(kind='release')))
    results = runner.run(jobs, 'final + release')
    for j in jobs:
        # the checkpoint of the current models directory, not the absolute
        # path stored in a cached fit result (a moved or re-cloned checkout)
        if j.get('save_model') and not results[j['key']].get('error'):
            results[j['key']]['saved_model'] = j['save_model']
    refresh_registries(jobs, results, s)
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
    out['release_conformal'] = out['release'].get('conformal') if not out['release']['error'] else None
    c = out['release_conformal'] or {}
    if c:
        log(f"release: intervals from {c.get('source')}: {c.get('n')} samples of {c.get('n_calib_groups')} "
            f"configurations, cross-fitted MPSP coverage {c.get('mpsp_coverage_cv', float('nan')):.3f} -> {rel}")
    # release check through the user-facing API: it must encode every
    # training row, and the by-name query of every training configuration,
    # exactly as the validated path (Ensemble.predict_samples)
    out['release_check'] = None
    if not out['release']['error']:
        try:
            from fermentation_insights.coeff_predictor.predict import CoefficientPredictor
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                P = CoefficientPredictor.load(rel)
                base = s[s['sample_kind'].astype(str) == 'baseline'].reset_index(drop=True)
                api = api_consistency(P, s, base)
                pr = P.predict(base.head(1), warn=False)
                q = P.predict(base[API_QUERY_FIELDS].head(1), warn=False)
                fdp = os.path.join(models_dir, 'final_dev')
                ood = s[s['sample_id'].astype(str).isin(set(map(str, sp.get('ood_test_sample_ids', []))))]
                if os.path.isdir(fdp) and len(ood):
                    api['new_configs'] = api_consistency_new_rows(CoefficientPredictor.load(fdp),
                                                                  ood.reset_index(drop=True))
            out['release_check'] = dict(ok=True, n=len(pr), columns=int(pr.shape[1]),
                                        example={k: float(pr[k].iloc[0]) for k in TARGETS},
                                        example_by_name={k: float(q[k].iloc[0]) for k in TARGETS},
                                        intervals=all(f'{k}_lo' in pr for k in TARGETS),
                                        example_config=str(base['config_id'].iloc[0]), **api)
            log(f'release: {rel} loads with CoefficientPredictor; prediction for {base["config_id"].iloc[0]} '
                f'baseline: ' + ', '.join(f'{k}={pr[k].iloc[0]:.4g}' for k in TARGETS)
                + f'; intervals: {out["release_check"]["intervals"]}; API vs. validated path: max rel. diff '
                f'{api["rows_max_rel"]:.2g} over {api["n_rows"]} training rows, by-name max rel. MPSP diff '
                f'{api["by_name_max_rel_mpsp"]:.2g} over {api["n_configs"]} configurations')
        except Exception as e:
            out['release_check'] = dict(ok=False, error=f'{type(e).__name__}: {e}')
            log(f'release: WARNING CoefficientPredictor check failed: {type(e).__name__}: {e}')
    return out


#: fields of a by-name query of a training configuration (release check)
API_QUERY_FIELDS = ['config_id', 'feedstock', 'fermentation_product', 'final_product', 'r', 'feedstock_price',
                    'electricity_price', 'natural_gas_price']


def api_consistency(P, samples, base):
    """Max relative difference between the user API
    (``CoefficientPredictor.predict``) and the validated path
    (``Ensemble.predict_samples``, behind every CV / test number): over all
    ``samples`` passed as full rows (a', b', c', d', g), and for the by-name
    query (:data:`API_QUERY_FIELDS`) of each configuration's baseline sample
    ``base`` (MPSP at its baseline yield and titer)."""
    import numpy as np

    def ref(rows):
        _, o = P.ensemble.predict_samples(rows, strict=False)
        return np.column_stack([o['coef'].mean(0), o['g'].mean(0)])

    def rel(a, b):
        return np.abs(a - b) / np.maximum(np.abs(b), 1e-9)

    full = P.predict(samples, warn=False)[list(TARGETS)].to_numpy(float)
    qd = P.predict(base[API_QUERY_FIELDS], warn=False)
    q = qd[list(TARGETS)].to_numpy(float)
    y, t = base['baseline_yield'].to_numpy(float), base['baseline_titer'].to_numpy(float)
    mp = lambda C: (C[:, 0] + C[:, 1] / y + C[:, 2] / t + C[:, 3] / (y * t)) / C[:, 4]
    out = dict(n_rows=int(len(samples)), rows_max_rel=float(np.nanmax(rel(full, ref(samples)))),
               n_configs=int(base['config_id'].nunique()),
               by_name_max_rel_mpsp=float(np.nanmax(rel(mp(q), mp(ref(base))))))
    # the applicability domain of the by-name query equals that of the
    # baseline row, and a query without prices uses the baseline prices
    if 'ad_ratio' in qd:
        fb = P.predict(base, warn=False)
        out['by_name_ad_max_abs'] = float(np.nanmax(np.abs(qd['ad_ratio'].to_numpy(float)
                                                           - fb['ad_ratio'].to_numpy(float))))
        out['by_name_in_domain_mismatch'] = int((qd['in_domain'].to_numpy() != fb['in_domain'].to_numpy()).sum())
    noprice = [c for c in API_QUERY_FIELDS if c not in ('electricity_price', 'natural_gas_price')]
    qn = P.predict(base[noprice], warn=False)[list(TARGETS)].to_numpy(float)
    out['by_name_no_prices_max_rel_mpsp'] = float(np.nanmax(rel(mp(qn), mp(ref(base)))))
    return out


def api_consistency_new_rows(P, rows):
    """The user API on rows of configurations NOT in the checkpoint's
    training data, passed as stored (e.g. ``final_dev`` on the held-out test
    configurations): max relative difference of a', b', c', d', g from the
    validated path and the number of rows whose ``in_domain`` differs from
    the applicability domain of the validated encoding."""
    import numpy as np
    b, o = P.ensemble.predict_samples(rows, strict=False)
    ref = np.column_stack([o['coef'].mean(0), o['g'].mean(0)])
    pr = P.predict(rows, warn=False)
    got = pr[list(TARGETS)].to_numpy(float)
    out = dict(n_rows=int(len(rows)), n_configs=int(rows['config_id'].nunique()),
               max_rel=float(np.nanmax(np.abs(got - ref) / np.maximum(np.abs(ref), 1e-9))))
    if 'in_domain' in pr and P.ensemble.ad is not None:
        ratio, oor = P._ad(b)
        dom = np.isfinite(ratio) & (ratio <= 1.) & np.array([not x for x in oor])
        out['in_domain_mismatch'] = int((pr['in_domain'].to_numpy() != dom).sum())
    return out


def refresh_registries(jobs, results, samples):
    """Spec registries (by-name defaults of the user API) of the saved
    checkpoints in the current layout (``predict.refresh_spec_registry``;
    their networks, encoding and intervals are unchanged): a checkpoint
    saved by an earlier run with a cached fit gets the per-configuration
    defaults. A no-op for checkpoints trained by this code."""
    from fermentation_insights.coeff_predictor.predict import refresh_spec_registry
    for j in jobs:
        path = j.get('save_model')
        if not path or results[j['key']].get('error') or not os.path.isdir(path):
            continue
        tr = samples[samples['sample_id'].astype(str).isin(set(map(str, j['train_ids'])))]
        try:
            if refresh_spec_registry(path, tr):
                log(f'final: spec registry of {os.path.relpath(path, PKG_DIR)} refreshed (per-configuration '
                    'defaults; model unchanged)')
        except ValueError as e:
            log(f'final: WARNING spec registry of {path} not refreshed: {e}')


def annotate_unseen_coverage(shc, paths):
    """Store the coverage of the shipped interval quantiles for a new
    feedstock / product (leave-one-feedstock-out / -product-out out-of-fold
    scores, pooled; ``report_checks.shipped_quantile_coverage``) in the
    checkpoints' ``meta.json`` (``meta['coverage_unseen']``), so that
    ``CoefficientPredictor`` warns with the measured numbers. Report-only
    metadata: the networks, quantiles and domain are unchanged."""
    cov = {sch: {r.key: round(float(r.pooled), 4) for r in shc[shc['scheme'] == sch].itertuples()}
           for sch in ('feedstock', 'product') if (shc['scheme'] == sch).any()}
    if not cov:
        return
    cov['source'] = ('coverage of the release quantiles on the leave-one-feedstock-out / leave-one-product-out '
                     'out-of-fold scores, pooled over samples (REPORT.md section 6, '
                     'tables/interval_coverage_by_scheme.csv)')
    for path in paths:
        mp = os.path.join(path, 'meta.json')
        if not os.path.exists(mp):
            continue
        with open(mp) as fh:
            m = json.load(fh)
        if (m.get('meta') or {}).get('coverage_unseen') == cov:
            continue
        m.setdefault('meta', {})['coverage_unseen'] = cov
        with open(mp, 'w') as fh:
            json.dump(m, fh, indent=1)
        log(f'tables: {os.path.relpath(path, PKG_DIR)}/meta.json: coverage for a new feedstock / product recorded')


#: interval schemes of ``meta['coverage_by_family']``: the CV schemes whose
#: out-of-fold scores measure the coverage for a known configuration (new
#: parameter draws; productivities beyond the sampled range) and for a new
#: configuration of a known product (leave-one-configuration-out)
FAMILY_COVERAGE_SCHEMES = {'known_configuration': ('random', 'productivity'), 'new_configuration': ('config',)}


def annotate_family_coverage(sc, samples, paths):
    """Store, per fermentation product (family) and interval key, the
    coverage of each checkpoint's OWN shipped quantiles on the development
    out-of-fold scores (``report_checks.family_coverage``; random CV and
    productivity tails = a known configuration, the minimum of the two;
    LOCO = a new configuration of a known product) in its ``meta.json``
    (``meta['coverage_by_family']``), so that ``CoefficientPredictor`` flags
    the keys that under-cover the row's family (``intervals_undercovered``).
    Report-only metadata: the networks, quantiles and domain are unchanged,
    and no test sample is read (the scores are development CV scores)."""
    import pandas as pd
    from fermentation_insights.coeff_predictor import report_checks as RC
    from fermentation_insights.coeff_predictor.train import CONF_KEYS
    fam = dict(zip(samples['config_id'].astype(str), samples['fermentation_product'].astype(str)))
    for path in paths:
        mp = os.path.join(path, 'meta.json')
        if not os.path.exists(mp):
            continue
        with open(mp) as fh:
            m = json.load(fh)
        conf = m.get('conformal') or {}
        if not conf:
            continue
        cov = dict(nominal=round(1 - float(conf.get('alpha', ALPHA)), 4))
        for name, schemes in FAMILY_COVERAGE_SCHEMES.items():
            per = {}
            for sch in schemes:
                f_, _ = RC.family_coverage(sc, conf, scheme=sch, keys=CONF_KEYS, families=fam)
                if not len(f_):
                    continue
                cov.setdefault('by_scheme', {})[sch] = {
                    str(fp): {x.key: round(float(x.coverage), 4) for x in g.itertuples()}
                    for fp, g in f_.groupby('family')}
                cov.setdefault('n_configs', {})[sch] = {
                    str(fp): int(g['n_configs'].max()) for fp, g in f_.groupby('family')}
                for x in f_.itertuples():
                    d = per.setdefault(str(x.family), {})
                    d[x.key] = round(min(float(x.coverage), d.get(x.key, float('inf'))), 4)
            if per:
                cov[name] = per
        cov['source'] = ('coverage of this checkpoint\'s interval quantiles on the development out-of-fold scores '
                         '(tables/cv_oof_scores.csv.gz), pooled over the held-out samples of each fermentation '
                         'product: known_configuration = min(random CV, productivity tails), new_configuration = '
                         'leave-one-configuration-out (REPORT.md section 6, '
                         'tables/interval_coverage_in_distribution_by_family.csv, '
                         'tables/interval_coverage_by_family.csv)')
        if (m.get('meta') or {}).get('coverage_by_family') == cov:
            continue
        m.setdefault('meta', {})['coverage_by_family'] = cov
        with open(mp, 'w') as fh:
            json.dump(m, fh, indent=1)
        log(f'tables: {os.path.relpath(path, PKG_DIR)}/meta.json: interval coverage per fermentation product recorded')


def annotate_channel_sum(samples, train_ids, masked=None):
    """Note on the ``channel_sum`` entry of the checkpoints' constraint
    report: the fit checks the identity (channel coefficients sum to a', b',
    c', d') on its training targets AFTER the defective channel targets were
    masked for the loss (``DEFECTIVE_TARGETS``, ``train.mask_channel_targets``),
    so the masked samples count as failing. Adds the identity on the
    unmasked targets (``frac_satisfied_unmasked``, max error relative to the
    coefficient spread) and a note. ``train_ids``: {checkpoint path: its
    training sample ids}. Report-only metadata."""
    import numpy as np
    import pandas as pd
    from fermentation_insights.coeff_predictor.train import COEFF_KEYS
    masked = DEFECTIVE_TARGETS if masked is None else masked
    chs = sorted({c.split('__')[1] for c in samples.columns if c.startswith('ch__') and c.count('__') == 2})
    if not chs:
        return
    col = lambda c: pd.to_numeric(samples[c], errors='coerce').to_numpy(float) if c in samples else \
        np.full(len(samples), np.nan)
    coef = np.column_stack([col(k) for k in COEFF_KEYS])
    S_ = np.column_stack([np.nansum(np.column_stack([col(f'ch__{c}__{k}') for c in chs]), 1) for k in COEFF_KEYS])
    fp = samples['fermentation_product'].astype(str).to_numpy() if 'fermentation_product' in samples else None
    for path, ids in train_ids.items():
        mp = os.path.join(path, 'meta.json')
        if not os.path.exists(mp):
            continue
        with open(mp) as fh:
            m = json.load(fh)
        rep = m.get('constraint_report') or []
        ent = [r for r in rep if r.get('constraint') == 'channel_sum']
        rows = np.flatnonzero(samples['sample_id'].astype(str).isin(set(map(str, ids))).to_numpy())
        if not ent or len(rows) != int(ent[0].get('n_cases') or -1):
            continue
        s_tot = np.nanstd(coef[rows], 0) + 1e-12
        err = (np.abs(S_[rows] - coef[rows]) / s_tot).max(1)
        ok = np.isfinite(err)
        nm = int(sum((fp[rows] == p).sum() for ps in masked.values() for p in ps)) if fp is not None else 0
        what = '; '.join(f'{k}: ' + ', '.join(v) for k, v in masked.items())
        upd = dict(frac_satisfied_unmasked=float((err[ok] <= 1e-4).mean()) if ok.any() else float('nan'),
                   max_rel_err_unmasked=float(err[ok].max()) if ok.any() else float('nan'),
                   note=(f'frac_satisfied is computed by the fit on its training targets after the defective '
                         f'channel targets were masked for the loss ({what}'
                         f'; {nm} samples, whose masked channel counts as 0 in the sum): it measures the mask, not '
                         f'the data. On the unmasked targets the identity holds for '
                         f'{int((err[ok] <= 1e-4).sum())} of {int(ok.sum())} samples (max error '
                         f'{float(err[ok].max()):.1e} of the coefficient spread)'))
        if all(ent[0].get(k) == v for k, v in upd.items()):
            continue
        ent[0].update(upd)
        with open(mp, 'w') as fh:
            json.dump(m, fh, indent=1)
        log(f'tables: {os.path.relpath(path, PKG_DIR)}/meta.json: channel_sum note recorded')


def report_only_checks(a, D, cv, fin, checks, S, Pts, tdir, fp):
    """Report-only statistics from the saved outputs (``report_checks``;
    none is used to fit, tune, select or calibrate anything): paired
    fold-level ablation statistics, final-test comparison with
    configuration-cluster bootstrap CIs, coverage of the shipped interval
    quantiles on each CV scheme's out-of-fold scores, applicability-domain
    validation on the test sets, per-configuration sample coverage, and the
    random-CV error on random samples per method."""
    import warnings
    from fermentation_insights.coeff_predictor import report_checks as RC, dataset
    checks['ablation_stats'], checks['ablation_per_fold'] = RC.ablation_stats(cv['results'], ABLATION_SCHEMES,
                                                                              ABLATIONS)
    checks['ablation_stats'].to_csv(os.path.join(tdir, 'cv_ablation_paired.csv'), index=False)
    checks['ablation_per_fold'].to_csv(os.path.join(tdir, 'cv_ablation_per_fold.csv'), index=False)
    checks['test_paired'] = RC.paired_test_comparison(fin.get('final'), S, Pts, seed=a.seed)
    checks['test_paired'].to_csv(os.path.join(tdir, 'final_test_paired.csv'), index=False)
    conf = fin.get('release_conformal') or {}
    sc = cv.get('oof_scores')
    checks['shipped_coverage'] = (RC.shipped_quantile_coverage(sc, conf) if conf and sc is not None and len(sc)
                                  else None)
    if checks['shipped_coverage'] is not None:
        checks['shipped_coverage'].to_csv(os.path.join(tdir, 'interval_coverage_by_scheme.csv'), index=False)
        annotate_unseen_coverage(checks['shipped_coverage'], [os.path.join(a.models, m)
                                                              for m in ('release', 'final_dev')])
        annotate_family_coverage(sc, D['samples'], [os.path.join(a.models, m) for m in ('release', 'final_dev')])
    annotate_channel_sum(D['samples'], {j['save_model']: j['train_ids'] for j in (fin.get('jobs') or ())
                                        if j.get('save_model')})
    checks['baseline_paired'] = RC.baseline_fold_stats(cv['results'], methods=[m for m in BASE_METHODS
                                                                               if m != 'pinn'])
    checks['baseline_paired'].to_csv(os.path.join(tdir, 'cv_baselines_paired.csv'), index=False)
    checks['family_coverage'] = checks['config_coverage'] = None
    checks['id_family_coverage'] = checks['id_config_coverage'] = None
    if conf and sc is not None and len(sc):
        fam = dict(zip(D['samples']['config_id'].astype(str), D['samples']['biorefinery'].astype(str)))
        from fermentation_insights.coeff_predictor.train import CONF_KEYS as _CK
        checks['family_coverage'], checks['config_coverage'] = RC.family_coverage(sc, conf, keys=_CK, families=fam)
        checks['family_coverage'].to_csv(os.path.join(tdir, 'interval_coverage_by_family.csv'), index=False)
        checks['config_coverage'].to_csv(os.path.join(tdir, 'interval_coverage_loco_per_config.csv'), index=False)
        # a KNOWN configuration: random-CV (new parameter draws) and
        # productivity-tail out-of-fold scores, per family and configuration
        import pandas as pd
        from fermentation_insights.coeff_predictor.train import CONF_KEYS
        fr, cr = [], []
        for sch in ('random', 'productivity'):
            f_, c_ = RC.family_coverage(sc, conf, scheme=sch, keys=CONF_KEYS, families=fam)
            if len(f_):
                fr.append(f_.assign(scheme=sch))
                cr.append(c_.assign(scheme=sch, family=c_['config_id'].map(fam)))
        if fr:
            checks['id_family_coverage'] = pd.concat(fr, ignore_index=True)[['scheme', 'key', 'family', 'n',
                                                                               'n_configs', 'coverage']]
            checks['id_config_coverage'] = pd.concat(cr, ignore_index=True)[['scheme', 'key', 'family',
                                                                               'config_id', 'coverage']]
            checks['id_family_coverage'].to_csv(os.path.join(tdir, 'interval_coverage_in_distribution_by_family.csv'),
                                                index=False)
            checks['id_config_coverage'].to_csv(os.path.join(tdir, 'interval_coverage_in_distribution_per_config.csv'),
                                                index=False)
    checks['ad'] = checks['ad_per_config'] = None
    path = (fp or {}).get('saved_model') or os.path.join(a.models, 'final_dev')
    if fp and not fp.get('error') and os.path.isdir(path):
        try:
            from fermentation_insights.coeff_predictor.predict import CoefficientPredictor
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                Pfd = CoefficientPredictor.load(path)
                checks['ad'], checks['ad_per_config'] = RC.ad_validation(
                    Pfd, {st: fp['evals'][st] for st in ('test_ood', 'test_id') if st in fp['evals']}, S, Pts)
            checks['ad'].to_csv(os.path.join(tdir, 'applicability_domain_test.csv'), index=False)
            checks['ad_per_config'].to_csv(os.path.join(tdir, 'applicability_domain_test_per_config.csv'),
                                           index=False)
        except Exception as e:
            log(f'tables: WARNING applicability-domain validation failed: {type(e).__name__}: {e}')
    try:
        raw_files = dataset._find_files(a.raw)
        cov_path = os.path.join(tdir, 'data_coverage.csv')
        if not raw_files and os.path.exists(cov_path):
            # no raw records (a clone of the repository): keep the table of the
            # run that had them rather than counting zero started samples
            import pandas as pd
            checks['data_coverage'] = pd.read_csv(cov_path)
            log(f'tables: no raw files under {a.raw}: sample coverage kept from {cov_path}')
        else:
            checks['data_coverage'] = RC.raw_sample_coverage(raw_files, D['samples'], PLANNED_SAMPLES)
            checks['data_coverage'].to_csv(cov_path, index=False)
    except Exception as e:
        checks['data_coverage'] = None
        log(f'tables: WARNING raw sample coverage failed: {type(e).__name__}: {e}')
    checks['price_response'] = {}
    for name, path_, rows in (('models/release, in-sample (all training rows)', os.path.join(a.models, 'release'),
                               D['samples']),
                              ('models/final_dev, test set', path,
                               S.loc[(fp or {}).get('evals', {}).get('test_all', {}).get('ids', [])])):
        if not os.path.isdir(path_) or not len(rows):
            continue
        try:
            from fermentation_insights.coeff_predictor.predict import CoefficientPredictor
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                t = RC.price_response_error(CoefficientPredictor.load(path_), rows)
            checks['price_response'][name] = t
            t.to_csv(os.path.join(tdir, f'price_response_{os.path.basename(path_)}.csv'), index=False)
        except Exception as e:
            log(f'tables: WARNING price response of {path_} failed: {type(e).__name__}: {e}')
    checks['capacity_response'] = None
    cap_csv = os.path.join(PKG_DIR, RC.CAPACITY_RESPONSE_CSV)
    rel_path = os.path.join(a.models, 'release')
    if os.path.exists(cap_csv) and os.path.isdir(rel_path):
        try:
            import pandas as pd
            from fermentation_insights.coeff_predictor.predict import CoefficientPredictor
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                t = RC.capacity_response(CoefficientPredictor.load(rel_path), pd.read_csv(cap_csv), D['samples'])
            checks['capacity_response'] = t
            t.to_csv(os.path.join(tdir, 'capacity_response_release.csv'), index=False)
        except Exception as e:
            log(f'tables: WARNING capacity response failed: {type(e).__name__}: {e}')
    r = cv['results']
    r = r[(r['scheme'] == 'random') & (r['fold'] >= 0) & (r['metric'] == 'mpsp_mape_cfg_gmean_random')]
    checks['random_cv_random_samples'] = r.groupby('method')['value'].mean().to_dict()
    clustered_and_twin_checks(a, D, cv, fin, checks, S, Pts, tdir, fp, path)


def clustered_and_twin_checks(a, D, cv, fin, checks, S, Pts, tdir, fp, saved_model):
    """Report-only checks of the CV and test designs (``report_checks``):
    the LOCO comparison with the baselines at the family level (the 37
    configurations cluster in 7 families), LOCO without the folds that hold
    out the only configuration of their product (new-product folds) and per
    family, and configurations whose model inputs are identical (input
    twins) -- across the development / test boundary they make a held-out
    configuration not new in input space, so the OOD test is also reported
    without them, and the LOCO folds whose twin is in training are marked."""
    import warnings
    import pandas as pd
    from fermentation_insights.coeff_predictor import report_checks as RC
    from fermentation_insights.coeff_predictor.train import CONF_KEYS
    res = cv['results']
    checks['yt_check'] = RC.yt_definition_check(D['samples'], D['points'])
    checks['yt_check'].to_csv(os.path.join(tdir, 'yt_definition_check.csv'), index=False)
    checks['sugar_check'] = RC.sugar_balance_check(D['samples'], D['points'])
    checks['sugar_check'].to_csv(os.path.join(tdir, 'sugar_balance_check.csv'), index=False)
    fam = dict(zip(D['samples']['config_id'].astype(str), D['samples']['biorefinery'].astype(str)))
    checks['baseline_family'] = RC.baseline_family_stats(res, fam, methods=[m for m in BASE_METHODS if m != 'pinn'])
    checks['baseline_family'].to_csv(os.path.join(tdir, 'cv_baselines_paired_family.csv'), index=False)
    loco = set(res.loc[res['scheme'] == 'config', 'label'].astype(str))
    newp = RC.single_config_products(D['samples'][D['samples']['config_id'].astype(str).isin(loco)])
    # the same over the known-product families only (the new-product folds
    # are single-configuration families of their own)
    res_k = res[(res['scheme'] != 'config') | ~res['label'].astype(str).isin(newp)]
    checks['baseline_family_known'] = RC.baseline_family_stats(res_k, fam, methods=[m for m in BASE_METHODS
                                                                                    if m != 'pinn'])
    checks['baseline_family_known'].insert(1, 'excluded_folds', ', '.join(newp))
    checks['baseline_family_known'].to_csv(os.path.join(tdir, 'cv_baselines_paired_family_known_product.csv'),
                                           index=False)
    checks['loco_new_product_folds'] = newp
    checks['loco_known_product'] = RC.loco_subset(res, newp, methods=BASE_METHODS)
    checks['loco_known_product'].insert(0, 'excluded_folds', ', '.join(newp))
    checks['loco_known_product'].to_csv(os.path.join(tdir, 'cv_loco_known_product.csv'), index=False)
    checks['loco_per_family'] = RC.loco_per_family(res, fam)
    fm = RC._fold_mape(res, 'config')
    checks['loco_per_fold_pinn'] = fm['pinn'].dropna().to_dict() if 'pinn' in fm else {}
    checks['loco_per_family'].to_csv(os.path.join(tdir, 'cv_loco_per_family.csv'), index=False)
    conf = fin.get('release_conformal') or {}
    sc = cv.get('oof_scores')
    checks['known_product_coverage'] = None
    if conf and sc is not None and len(sc):
        grp = {c: ('new-product folds' if c in newp else 'known-product folds') for c in fam}
        checks['known_product_coverage'], _ = RC.family_coverage(sc, conf, keys=CONF_KEYS, families=grp)
        checks['known_product_coverage'].rename(columns={'family': 'folds'}).to_csv(
            os.path.join(tdir, 'interval_coverage_known_product.csv'), index=False)
    # input twins (features of models/final_dev)
    for k in ('twins', 'twin_nearest', 'test_no_twin', 'test_no_twin_cov', 'ad_no_twin'):
        checks[k] = None
    checks['ood_twins'], checks['loco_twin_folds'] = [], []
    if not (fp and not fp.get('error') and os.path.isdir(saved_model) and 'test_ood' in fp.get('evals', {})):
        return
    try:
        from fermentation_insights.coeff_predictor.predict import CoefficientPredictor
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            Pfd = CoefficientPredictor.load(saved_model)
            ood = set(S.loc[list(fp['evals']['test_ood']['ids']), 'config_id'].astype(str))
            sets = {c: ('test_ood' if c in ood else 'development') for c in fam}
            tw, near = RC.input_twins(Pfd, D['samples'], sets)
            checks['twins'], checks['twin_nearest'] = tw, near
            tw.to_csv(os.path.join(tdir, 'input_twins.csv'), index=False)
            cross = tw[tw['crosses_split']]
            oodt = sorted({c for x in cross.itertuples() for c, st in ((x.config_id, x.set), (x.twin, x.twin_set))
                           if st == 'test_ood'})
            devt = tw[(tw['set'] == 'development') & (tw['twin_set'] == 'development')]
            checks['ood_twins'] = oodt
            checks['loco_twin_folds'] = sorted(set(devt['config_id']) | set(devt['twin']))
            if oodt:
                log(f'tables: WARNING held-out test configuration(s) {oodt} encode identically to development '
                    'configuration(s) (input twins; tables/input_twins.csv): the OOD test is also reported without '
                    'them')
                checks['test_no_twin'], checks['test_no_twin_cov'] = RC.final_test_without_configs(
                    fin.get('final'), S, Pts, oodt, seed=a.seed)
                checks['test_no_twin'].to_csv(os.path.join(tdir, 'final_test_without_twins.csv'), index=False)
                checks['test_no_twin_cov'].to_csv(os.path.join(tdir, 'final_intervals_without_twins.csv'),
                                                  index=False)
                ev = RC.drop_configs(fp['evals']['test_ood'], S, oodt)
                checks['ad_no_twin'], _ = RC.ad_validation(Pfd, {'test_ood without the input twins': ev}, S, Pts)
                checks['ad_no_twin'].to_csv(os.path.join(tdir, 'applicability_domain_test_without_twins.csv'),
                                            index=False)
    except Exception as e:
        log(f'tables: WARNING input-twin check failed: {type(e).__name__}: {e}')
    checks['twin_bias'] = None
    rel = os.path.join(a.models, 'release')
    if checks.get('twins') is not None and len(checks['twins']) and os.path.isdir(rel):
        try:
            from fermentation_insights.coeff_predictor.predict import CoefficientPredictor
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                tw = checks['twins']
                checks['twin_bias'] = RC.baseline_bias(CoefficientPredictor.load(rel), D['samples'],
                                                       sorted(set(tw['config_id']) | set(tw['twin'])))
            checks['twin_bias'].to_csv(os.path.join(tdir, 'release_twin_bias.csv'), index=False)
        except Exception as e:
            log(f'tables: WARNING twin bias of the release model failed: {type(e).__name__}: {e}')
    if checks['loco_twin_folds']:
        checks['loco_matched_no_twin'] = loco_matched(res, S, exclude=checks['loco_twin_folds'])
        checks['loco_matched_no_twin'].to_csv(os.path.join(tdir, 'cv_loco_ood_families_without_twin_folds.csv'),
                                              index=False)


def data_violation_table(dev, plan, seed, fin, S, Pts):
    """Constraint-violation rates of the data themselves (fitted
    coefficients; ``report_checks.data_violation_rates``) on the held-out
    samples of every CV scheme (mean over its folds, as the methods' rates)
    and of every final test set: what a method that reproduces the data
    exactly would show."""
    import numpy as np
    import pandas as pd
    from fermentation_insights.coeff_predictor import report_checks as RC
    ids = dev['sample_id'].astype(str).to_numpy()
    rows = []
    for sch in SCHEMES:
        per = [RC.data_violation_rates(S, Pts, list(ids[te])) for _, _, te, _ in scheme_splits(dev, sch, plan, seed)]
        keys = sorted({k for p in per for k in p})
        rows.append(dict(kind='cv', set=sch, n_folds=len(per),
                         **{k: float(np.nanmean([p.get(k, np.nan) for p in per])) for k in keys}))
    fp = (fin.get('final') or {}).get('pinn') or {}
    for st, ev in (fp.get('evals') or {}).items():
        if st in ('test_ood', 'test_id', 'test_all'):
            rows.append(dict(kind='test', set=st, n_folds=1, **RC.data_violation_rates(S, Pts, list(ev['ids']))))
    return pd.DataFrame(rows)


def final_intervals(fin, samples_by_id, alpha=ALPHA):
    """Realized test coverage of the final models' intervals (report only;
    never used for any choice): per model x test set x key, the
    sample-simultaneous coverage (a sample counts as covered only if its
    score <= q, i.e. for MPSP every TY point is inside its interval) pooled,
    averaged over configurations and in the worst configuration, the
    pointwise MPSP coverage, the mean / median relative half-width, next to
    the cross-fitted estimate of the calibration (fold-balanced mean,
    sample-pooled, min, q10 over the LOCO folds: compare ``coverage`` with
    ``crossfit_pooled`` and ``coverage_cfg_mean`` with ``crossfit_mean``).
    ``calibration``: 'loco' (out-of-fold LOCO scores), 'holdout' (own
    calibration hold-out, ``--v1-protocol``) or a reporting-only calibration
    (e.g. 'random_cv'). Returns (summary, per-configuration table)."""
    import numpy as np
    import pandas as pd
    from fermentation_insights.coeff_predictor.train import CONF_KEYS
    rows, per = [], []
    cals = []
    for meth, r in (fin.get('final') or {}).items():
        if r is None or r.get('error'):
            continue
        if r.get('conformal'):
            # 'loco': out-of-fold LOCO scores (v2); 'holdout': own calibration hold-out (--v1-protocol)
            lab = 'holdout' if r['conformal'].get('source') == 'calibration hold-out' else 'loco'
            cals.append((meth, r, lab, r['conformal']))
        for name, c in (r.get('conformal_extra') or {}).items():
            cals.append((meth, r, name, c))
    for meth, r, cal, conf in cals:
        for st, ev in r['evals'].items():
            if ev.get('scores') is None:
                continue
            sc = pd.DataFrame(ev['scores'])
            mm = {(m, t): v for m, t, v in ev['rows']}
            cfg = sc['config_id'].astype(str).to_numpy()
            for j, key in enumerate(CONF_KEYS):
                q = conf.get(key)
                v = sc[key].to_numpy(float)
                ok = ~np.isnan(v)
                if q is None or not ok.any():
                    continue
                cov = v[ok] <= q
                byc = pd.Series(cov).groupby(cfg[ok]).mean()
                cf = (conf.get('crossfit') or {}).get(key, {})
                if key == 'mpsp':
                    hw_mean = mm.get(('interval_halfwidth_rel', 'MPSP'), np.nan)
                    hw_med = np.nan
                else:
                    if key == 'g':
                        mu, sd = np.asarray(ev['g'], float), np.asarray(ev.get('g_std', np.nan), float)
                    else:
                        mu = np.asarray(ev['coef'], float)[:, j]
                        sd = np.asarray(ev.get('coef_std'), float)[:, j] if ev.get('coef_std') is not None \
                            else np.full(len(mu), np.nan)
                    with np.errstate(divide='ignore', invalid='ignore'):
                        rel = q * (sd + conf.get(f'{key}_eps', 0.)) / np.abs(mu)
                    rel = rel[np.isfinite(rel)]
                    hw_mean = float(np.mean(rel)) if len(rel) else np.nan
                    hw_med = float(np.median(rel)) if len(rel) else np.nan
                rows.append(dict(method=meth, calibration=cal, set=st, key=key, n=int(ok.sum()),
                                 n_configs=int(byc.size),
                                 q=float(q), level=1 - alpha, coverage=float(cov.mean()),
                                 coverage_cfg_mean=float(byc.mean()), coverage_cfg_min=float(byc.min()),
                                 coverage_pointwise=(mm.get(('coverage', 'MPSP'), np.nan) if key == 'mpsp'
                                                     else float(cov.mean())),
                                 halfwidth_rel_mean=hw_mean, halfwidth_rel_median=hw_med,
                                 crossfit_mean=cf.get('mean', np.nan), crossfit_pooled=cf.get('pooled', np.nan),
                                 crossfit_min=cf.get('min', np.nan),
                                 crossfit_q10=cf.get('q10', np.nan), source=conf.get('source'),
                                 weighting=conf.get('weighting'), n_calib=conf.get('n'),
                                 n_calib_groups=conf.get('n_calib_groups'),
                                 coverage_guaranteed=conf.get('coverage_guaranteed')))
                if st != 'test_all':
                    per += [dict(method=meth, calibration=cal, set=st, key=key, config_id=c,
                                 n=int((cfg[ok] == c).sum()), coverage=float(v_), crossfit_mean=cf.get('mean'))
                            for c, v_ in byc.items()]
    return pd.DataFrame(rows), pd.DataFrame(per)


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


def loco_matched(res, samples_by_id, families=OOD_FAMILIES, exclude=()):
    """LOCO CV restricted to the configurations of the families that have a
    held-out test configuration: the matched CV estimate of the OOD test
    (which holds out ONE configuration of each family). Per method: the
    family-balanced mean (mean over families of the mean over their LOCO
    folds; comparable to the configuration-balanced OOD test MAPE), the
    configuration mean, and one column per family. ``exclude``: fold labels
    left out (e.g. the folds whose held-out configuration has an
    input-identical twin in training)."""
    import numpy as np
    import pandas as pd
    r = res[(res['scheme'] == 'config') & (res['fold'] >= 0) & (res['metric'] == 'mpsp_mape')
            & (res['target'] == 'MPSP') & ~res['label'].astype(str).isin([str(x) for x in exclude])]
    if not len(r):
        return pd.DataFrame()
    bio = samples_by_id.groupby('config_id')['biorefinery'].first()
    r = r.assign(family=r['label'].map(bio))
    r = r[r['family'].isin(families)]
    rows = []
    for meth, g in r.groupby('method', sort=False):
        fam = g.groupby('family')['value'].mean()
        rows.append(dict(method=meth, n_configs=int(g['label'].nunique()), family_mean=float(fam.mean()),
                         config_mean=float(g['value'].mean()),
                         **{f: float(fam.get(f, np.nan)) for f in families}))
    out = pd.DataFrame(rows)
    out['_m'] = out['method'].map({m: i for i, m in enumerate(METHOD_ORDER)})
    return out.sort_values('_m').drop(columns='_m').reset_index(drop=True)


def protocol_table(cv, fin, hpo, plan):
    """Training protocol in practice, per stage (HPO, CV per scheme,
    final_dev, release) and network method: member fits, best-epoch
    distribution (median, min, max; members whose best epoch is <= 20 and
    members that ran into the epoch limit), epochs run and CPU time."""
    import numpy as np
    import pandas as pd
    groups = {}

    def add(stage, method, r):
        if r is None or r.get('error') or not r.get('epochs'):
            return
        g = groups.setdefault((stage, method), dict(best=[], run=[], sec=[], fits=0))
        g['fits'] += 1
        g['run'] += list(r['epochs'])
        g['best'] += list(r.get('best_epochs') or r['epochs'])
        g['sec'].append(float(r.get('seconds', np.nan)))

    for j in (hpo or {}).get('jobs', []):
        add('hpo', j['method'], hpo['raw'].get(j['key']))
    for j in cv.get('jobs', []):
        if j['method'].split('-')[0] in NN_BASES:
            add(f"cv {j['meta']['scheme']}", j['method'], cv['raw'].get(j['key']))
    for j in fin.get('jobs', []):
        if j['method'] in NN_BASES:
            add('final_dev' if j['meta']['kind'] == 'final' else 'release', j['method'], fin['raw'].get(j['key']))
    rows = []
    for (st, m), g in groups.items():
        b, e = np.asarray(g['best'], float), np.asarray(g['run'], float)
        rows.append(dict(stage=st, method=m, fits=g['fits'], members=len(b), best_epoch_median=float(np.median(b)),
                         best_epoch_min=int(b.min()), best_epoch_max=int(b.max()),
                         members_best_le_20=int((b <= 20).sum()),
                         members_at_max_epochs=int((e >= plan['epochs']).sum()),
                         epochs_run_median=float(np.median(e)), cpu_hours=float(np.nansum(g['sec']) / 3600),
                         seconds_per_fit=float(np.nanmean(g['sec']))))
    return pd.DataFrame(rows)


#: earlier full runs of the pipeline that evaluated on the frozen test split
#: (report comparison; their report and headline tables are kept in
#: ``results/<run>/``): (label, provenance)
PRIOR_RUNS = (('v1', 'commit fd3f386'),
              ('v2', 'commit 725ae11 (report regenerated from cached fits in b7b0f53)'),
              ('v3', 'commit 2060257 (report regenerated from cached fits in 056c0fc)'))
#: label of the run this code produces (the fourth use of the frozen test split)
THIS_RUN = 'v4'
#: every run that evaluated on the frozen test split, this one last
ALL_RUNS = tuple(r for r, _ in PRIOR_RUNS) + (THIS_RUN,)
#: title of the report's changelog section
CHANGES_TITLE = 'Changes ' + ' -> '.join(ALL_RUNS)
#: how many times this run uses the frozen test split ('fourth')
TEST_USE = ('first', 'second', 'third', 'fourth', 'fifth', 'sixth', 'seventh', 'eighth')[len(ALL_RUNS) - 1]
#: the data and code of the v3 run (results/v3; DESIGN.md section 1,
#: "Previous data set"): the v3 -> v4 changelog compares this run with them
V3_RUN = dict(data_hash='c345efda41c4fdca', code_version='5', n_samples=1600, n_points=35583, n_dev=1353,
              n_dev_configs=37, fit_mape_median=0.0132)
#: plan keys that define the training protocol and the search (v3 -> v4 check)
PLAN_KEYS = ('hpo_k', 'hpo_candidates', 'hpo_members', 'hpo_lopo_weight', 'cv_members', 'final_members',
             'release_members', 'random_k', 'loco_folds', 'ablation_loco_folds', 'epochs', 'patience', 'min_epochs',
             'batch_size')


def _and(items):
    """'a', 'a and b', 'a, b and c'."""
    items = [str(x) for x in items]
    return ' and '.join(items) if len(items) < 3 else ', '.join(items[:-1]) + ' and ' + items[-1]


def load_run_tables(out_dir, run, fallback=True):
    """Headline tables of an earlier run ``run`` (a label of PRIOR_RUNS;
    ``results/<run>/tables``) for the side-by-side comparison of the report:
    {name: DataFrame, 'dir': path}. Looks in ``<out_dir>/<run>/tables``, then
    (``fallback``; not in quick mode, whose numbers are not comparable with
    the full runs) in the package's ``results/<run>/tables``; {} if absent."""
    import pandas as pd
    dirs = [os.path.join(out_dir, run, 'tables')] + ([os.path.join(DEFAULT_RESULTS, run, 'tables')] if fallback
                                                      else [])
    for d in dirs:
        if os.path.isdir(d):
            res = {}
            for name in ('final_test_headline', 'cv_headline', 'final_test_per_config', 'cv_lopo_per_product',
                         'timings', 'timings_fits', 'hpo_candidates'):
                f = os.path.join(d, f'{name}.csv')
                if os.path.exists(f):
                    try:
                        res[name] = pd.read_csv(f)
                    except Exception:
                        pass
            # selected hyperparameters and plan of that run ({'selected', 'plan', ...})
            f = os.path.join(d, 'selected_hyperparameters.json')
            if os.path.exists(f):
                try:
                    with open(f) as fh:
                        res['selected_hyperparameters'] = json.load(fh)
                except Exception:
                    pass
            res['dir'] = d
            return res
    return {}


def load_v1_tables(out_dir, fallback=True):
    """Headline tables of the v1 run (commit fd3f386; ``results/v1/tables``);
    see :func:`load_run_tables`."""
    return load_run_tables(out_dir, 'v1', fallback)


# %% Stage 6: checks

#: publication comparison: a coefficient ratio ours / publication is shown
#: only where both values have the same sign and the smaller |value| is at
#: least this fraction of the median |publication value| of the coefficient
PUB_RATIO_MIN = 0.2


def _nanmedian(v):
    """Median of the finite values (NaN, without a warning, if none)."""
    import numpy as np
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(np.median(v)) if len(v) else float('nan')


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
    out = pd.DataFrame(rows)
    # a ratio of two near-zero (or opposite-sign) coefficients is undefined
    # (e.g. the TAL a', |a'| <= 0.16): NaN where the smaller |value| is below
    # PUB_RATIO_MIN x the median |publication value| of that coefficient
    for k in TARGETS:
        if len(out):
            o, q = out[f'{k}_ours'].to_numpy(float), out[f'{k}_pub'].to_numpy(float)
            bad = (np.minimum(np.abs(o), np.abs(q)) < PUB_RATIO_MIN * np.nanmedian(np.abs(q))) | (np.sign(o) != np.sign(q))
            out.loc[bad, f'{k}_ratio'] = np.nan
    return out


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


def fig_ablation(path, stats, per, made):
    """Paired fold differences (variant - full PI-NN fold MPSP MAPE, pp) of
    the ablations and the plain MLP, on the folds where the ablations ran."""
    import numpy as np
    plt = _style()
    schemes = [x for x in ABLATION_SCHEMES if (stats['scheme'] == x).any()]
    fig, axs = plt.subplots(1, max(len(schemes), 1), figsize=(5.6 * max(len(schemes), 1), 3.6), squeeze=False)
    for ax, sch in zip(axs[0], schemes):
        st = stats[stats['scheme'] == sch]
        abl = st[st['variant'].isin(ABLATIONS)]
        n = int(abl['n_folds'].max()) if len(abl) else 0
        sub = 'ablation folds' if (st['subset'] == 'ablation folds').any() else 'all folds of the variant'
        rows = [(m, 'all folds of the variant') for m in ABLATIONS if (abl['variant'] == m).any()] + \
            [('plain', sub)] * bool(((st['variant'] == 'plain') & (st['subset'] == sub)).any())
        y = np.arange(len(rows))[::-1]
        for yy, (m, ss) in zip(y, rows):
            d = 100 * per[(per['scheme'] == sch) & (per['variant'] == m) & (per['subset'] == ss)]['delta'].to_numpy()
            x = st[(st['variant'] == m) & (st['subset'] == ss)].iloc[0]
            jit = np.linspace(-0.18, 0.18, len(d)) if len(d) > 1 else np.zeros(len(d))
            ax.scatter(d, yy + jit, s=12, color=METHOD_COLOR[m], alpha=0.75, edgecolors='white', linewidths=0.3,
                       zorder=3)
            ax.errorbar([100 * x['mean']], [yy], xerr=[1.96 * 100 * np.nan_to_num(x['se'])], fmt='D', ms=5,
                        color=INK, ecolor=INK, elinewidth=1.2, capsize=3, zorder=4)
            ax.text(1.0, yy, f"  {100 * x['mean']:+.1f} ± {100 * np.nan_to_num(x['se']):.1f} pp; better on "
                    f"{int(x['n_variant_better'])}/{int(x['n_folds'])}", transform=ax.get_yaxis_transform(),
                    va='center', fontsize=7.5, color=INK2, clip_on=False)
        ax.axvline(0, color=INK, lw=1)
        ax.set_yticks(y)
        ax.set_yticklabels([METHOD_LABEL[m].replace('PI-NN ', 'PI-NN\n').replace('Plain MLP', 'Plain\nMLP')
                            for m, _ in rows])
        ax.set_xlabel('variant - full PI-NN fold MPSP MAPE [pp]\n(dots = folds; diamond = mean ± 1.96 s.e.)')
        lab = _slabel(sch).replace('\n', '')
        ax.set_title(f'{lab}: {n} ' + ('stratified LOCO folds' if sch == 'config' else 'folds'), loc='left')
        ax.set_xscale('symlog', linthresh=5)
        ax.grid(axis='y', visible=False)
    fig.tight_layout()
    fig.subplots_adjust(right=0.80 if len(schemes) == 1 else 0.86, wspace=0.75)
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
            med = _nanmedian(pub[f'{k}_ratio'])
            nr = int(np.isfinite(pub[f'{k}_ratio'].to_numpy(float)).sum())
            ax.text(0.04, 0.96, f'median ratio = {med:.3f}\nn = {f.sum()}'
                    + (f' ({nr} defined ratios)' if nr != f.sum() else ''), transform=ax.transAxes, va='top',
                    fontsize=8, color=INK2)
        ax.set_xlabel(f'publication {TARGET_LABEL[k]}')
        ax.set_ylabel(f'this work {TARGET_LABEL[k]}')
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=len(l), bbox_to_anchor=(0.5, -0.06))
    fig.suptitle('Baseline coefficients: this work (current BioSTEAM, harness v2-v4, D3-patched, D16 / D17 '
                 'corrected) vs. publication',
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


def ablation_results(res, ablation_folds):
    """Per-fold CV results of the ablation schemes restricted to the folds
    the ablations ran on (LOCO: ``ablation_folds``; None = all), for a
    like-for-like ablation figure."""
    r = res[res['scheme'].isin(ABLATION_SCHEMES) & (res['fold'] >= 0)]
    if ablation_folds is not None:
        keep = set(map(str, ablation_folds))
        r = r[(r['scheme'] != 'config') | r['label'].astype(str).isin(keep)]
    return r.reset_index(drop=True)


def _round_up(x, sig=3):
    """``x`` > 0 rounded UP to ``sig`` significant digits (a printed upper
    bound that holds: 1.4032 -> '1.41')."""
    import math
    if not math.isfinite(x) or x <= 0:
        return f'{x:g}'
    e = math.floor(math.log10(x)) - sig + 1
    return f'{math.ceil(x / 10 ** e - 1e-9) * 10 ** e:.{max(-e, 0)}f}'


def _round_down(x, sig=2):
    """``x`` > 0 rounded DOWN to ``sig`` significant digits (a printed lower
    bound that holds: 0.0771 -> '0.077')."""
    import math
    if not math.isfinite(x) or x <= 0:
        return f'{x:g}'
    e = math.floor(math.log10(x)) - sig + 1
    return f'{math.floor(x / 10 ** e + 1e-9) * 10 ** e:.{max(-e, 0)}f}'


def _abl_label(m):
    return METHOD_LABEL[m].replace('PI-NN ', '')


def ablation_summary(st):
    """Summary bullets of the paired fold-level ablation statistics
    (``report_checks.ablation_stats``)."""
    import numpy as np
    if st is None or not len(st):
        return ''
    L = []
    own = st[st['subset'] == 'all folds of the variant']
    for sch in ABLATION_SCHEMES:
        r = own[(own['scheme'] == sch) & own['variant'].isin(ABLATIONS)]
        if not len(r):
            continue
        n = int(r['n_folds'].max())
        L.append(f'* Ablations ({_slabel(sch)}; paired over the {n} folds where they ran; ablation - full PI-NN '
                 'fold MAPE: mean ± s.e., median, folds where the ablation is better, the fold with the largest '
                 'difference): ' + '; '.join(
                     f'{_abl_label(x.variant)} {100 * x.mean:+.1f} ± {100 * x.se:.1f} pp (median '
                     f'{100 * x.median:+.1f}; better on {x.n_variant_better}/{x.n_folds}; {x.dominant_fold} '
                     f'{100 * x.dominant_delta:+.1f})' for x in r.itertuples()) + '.')
    a = own[own['variant'].isin(ABLATIONS)]
    sig = a[(a['t'].abs() >= 2) & (a['p_wilcoxon'] < 0.05)]
    if len(a) and not len(sig):
        L.append(f'* No ablation effect is established at the fold level (|t| <= {_round_up(np.nanmax(a["t"].abs()))}, '
                 f'Wilcoxon p >= {_round_down(np.nanmin(a["p_wilcoxon"]))}; the means are driven by single folds). The '
                 'channel (intermediate-data) supervision, the price-explicit structure and the physics residuals / '
                 'collocation are kept for what they guarantee by construction -- exact linearity in the prices, '
                 'dMPSP/dr <= 0, coefficient bounds, an additive channel breakdown -- not for a measured accuracy '
                 'gain.')
    elif len(sig):
        L.append('* Ablation effects established at the fold level (|t| >= 2, Wilcoxon p < 0.05): '
                 + ', '.join(f'{_abl_label(x.variant)} ({_slabel(x.scheme)})' for x in sig.itertuples()) + '.')
    pl = own[(own['variant'] == 'plain') & (own['scheme'] == 'config')]
    if len(pl):
        x = pl.iloc[0]
        L.append('* Fold-level evidence for the physics-informed model as a whole: PI-NN vs. the plain MLP (same '
                 'inputs and training protocol; it reuses the hyperparameters selected for the PI-NN) on '
                 f'leave-one-configuration-out, plain - PI-NN {100 * x["mean"]:+.1f} ± {100 * x["se"]:.1f} pp, '
                 f'PI-NN better on {int(x["n_variant_worse"])}/{int(x["n_folds"])} folds (Wilcoxon p '
                 f'{x["p_wilcoxon"]:.2g}).')
    return '\n'.join(L)


def baseline_paired_clause(bp, sch, best):
    """'; paired over the folds ...' clause of a CV summary bullet: the
    fold-level comparison of the PI-NN with the best baseline and the plain
    MLP (``report_checks.baseline_fold_stats``) and whether an advantage over
    the best baseline is established (mean > 0 and Wilcoxon p < 0.05)."""
    import numpy as np
    if bp is None or not len(bp) or sch not in ('config', 'feedstock', 'product'):
        return ''
    r = bp[bp['scheme'] == sch].set_index('method')
    ms = [m for m in dict.fromkeys([best, 'plain']) if m in r.index]
    if not ms:
        return ''

    def one(m):
        x = r.loc[m]
        p = (f'Wilcoxon p {x["p_wilcoxon"]:.2g}, ' if np.isfinite(x['p_wilcoxon']) else '') \
            + f'sign test p {x["p_sign"]:.2g}'
        return (f'{METHOD_LABEL[m]} {100 * x["mean"]:+.1f} ± {100 * x["se"]:.1f} pp (median '
                f'{100 * x["median"]:+.1f}; PI-NN better on {int(x["n_pinn_better"])}/{int(x["n_folds"])}; {p}; '
                f'without {x["dominant_fold"]} {100 * x["mean_without_dominant"]:+.1f})')

    x = r.loc[ms[0]]
    ok = x['mean'] > 0 and np.isfinite(x['p_wilcoxon']) and x['p_wilcoxon'] < 0.05
    n = int(x['n_folds'])
    concl = (f'the advantage over {METHOD_LABEL[ms[0]]} holds at the fold level' if ok else
             f'no advantage of the PI-NN over {METHOD_LABEL[ms[0]]} is established at the fold level'
             + (f' ({n} folds, too few for a signed-rank test)' if n < 5 else ''))
    return (f'; paired over the {n} folds (baseline - PI-NN fold MAPE): ' + ', '.join(one(m) for m in ms)
            + f' -- {concl}')


def _test_vs_baseline(tp, st):
    """' vs. best baseline ...' clause of a final-test summary bullet
    (``report_checks.paired_test_comparison``)."""
    if tp is None or not len(tp):
        return ''
    r = tp[(tp['set'] == st) & (tp['method'] != 'pinn')]
    if not len(r):
        return ''
    b = r.loc[r['mpsp_mape'].idxmin()]
    n, lab = int(b['n_configs']), METHOD_LABEL.get(b['method'], b['method'])
    txt = (f' vs. best baseline {lab} {100 * b["mpsp_mape"]:.1f} % (PI-NN - {lab} '
           f'{100 * b["diff_reference_minus_method"]:+.1f} pp, 95 % configuration-cluster bootstrap CI '
           f'[{100 * b["ci_lo"]:+.1f}, {100 * b["ci_hi"]:+.1f}] over {n} configurations)')
    if n < 10:
        txt += f'; with {n} held-out configurations this test cannot rank the methods'
    return txt


def ablation_deltas(res):
    """{scheme: dict(n_folds, delta={ablation: mean over the common folds of
    (ablation - PI-NN) fold MAPE})}: the ablations may run on a subset of the
    LOCO folds, so they are compared with the PI-NN on the same folds."""
    out = {}
    r = res[(res['fold'] >= 0) & (res['metric'] == 'mpsp_mape') & (res['target'] == 'MPSP')]
    for sch in ABLATION_SCHEMES:
        g = r[r['scheme'] == sch].pivot_table(index='fold', columns='method', values='value', aggfunc='mean')
        if 'pinn' not in g:
            continue
        d, n = {}, 0
        for m in ABLATIONS:
            if m in g:
                ok = g[[m, 'pinn']].dropna()
                if len(ok):
                    d[m] = float((ok[m] - ok['pinn']).mean())
                    n = max(n, len(ok))
        if d:
            out[sch] = dict(n_folds=n, delta=d)
    return out


def _regen_test_clause(regen):
    """'k of the n out-of-distribution test samples, those of <configs>, and k
    of the n in-distribution ones' (the v3 -> v4 corrected test rows)."""
    if not regen:
        return 'the test samples of the regenerated configurations'
    oc = f', those of {_and(regen["ood_configs"])},' if regen['ood_configs'] else ''
    return (f'{regen["test_ood"]} of the {regen["n_test_ood"]} out-of-distribution test samples{oc} and '
            f'{regen["test_id"]} of the {regen["n_test_id"]} in-distribution ones')


def runs_comparison(fh, head, prior, regen=None, dev3=()):
    """Markdown: the final-test and CV headline numbers of the earlier runs
    (``prior`` = {run: tables of :func:`load_run_tables`}, PRIOR_RUNS) and of
    this run (THIS_RUN) side by side; ``regen``: ``report_checks.regenerated_rows``;
    ``dev3``: :func:`settings_vs_v3`."""
    import pandas as pd
    runs = [r for r, _ in PRIOR_RUNS]
    allr = runs + [THIS_RUN]
    L = [f'## {" / ".join(allr)} headline numbers\n']
    have = [r for r in runs if (prior.get(r) or {}).get('final_test_headline') is not None]
    if not have or not len(fh):
        L.append('Tables of the earlier runs not found (`results/<run>/tables/`); no side-by-side comparison.\n')
        return '\n'.join(L)
    tabs = {r: prior[r]['final_test_headline'] for r in have}
    tabs[THIS_RUN] = fh
    g = lambda d, c: float(d[c].iloc[0]) if len(d) and c in d else float('nan')
    rows = []
    for st in ('test_ood', 'test_id', 'test_all'):
        for m in BASE_METHODS:
            sel = {r: t[(t['set'] == st) & (t['method'] == m)] for r, t in tabs.items()}
            if not any(len(d) for d in sel.values()):
                continue
            row = dict(set=st, method=m)
            for c, lab in (('mpsp_mape', 'mape'), ('mpsp_mape_median', 'median'), ('R2_g', 'R2_g'),
                           ('n_test', 'n_test')):
                for r in allr:
                    if r in sel:
                        row[f'{r}_{lab}'] = g(sel[r], c)
            rows.append(row)
    T = pd.DataFrame(rows)
    for c in [c for c in T if c.endswith('_n_test')]:
        T[c] = T[c].map(lambda v: '' if v != v else f'{int(v)}')
    prov = '; '.join(f'{r} = {p}' for r, p in PRIOR_RUNS)
    L.append(f'**Final test** (the same frozen split of configurations and sample ids, used for the {TEST_USE} '
             f'time: the {THIS_RUN} numbers are likely optimistic, see the notice above; {prov}; MAPE in %):\n')
    L.append(md_table(T, pct_cols=tuple(c for c in T if c.endswith(('_mape', '_median')))) + '\n')
    L.append('v1 fit every method, the baselines included, without its 6 calibration configurations (15 % of the '
             'development configurations, configuration-grouped calibration hold-out of the v1 protocol); v2 and '
             'later runs fit every method on all development rows. So the v1 -> v2 change of the baselines '
             '(ridge, HGB, kNN, plain MLP) comes from the protocol, not from the methods. Between v2 and v3 the data '
             'were corrected (D15, H1) and the PI-NN changed (round 5); between v3 and ' + THIS_RUN
             + (' only the data were corrected' if not dev3 else ' the data were corrected and the settings differ ('
                + ', '.join(f'`{f}`' for f in dev3) + ')')
             + f' (D16, D17; "{CHANGES_TITLE}"). Each correction replaced the training rows AND '
             'the test rows of the regenerated configurations: v2 -> v3 also their number (`n_test`), v3 -> '
             f'{THIS_RUN} the simulated values of the same samples ({_regen_test_clause(regen)}). So the v2 -> v3 '
             f'and v3 -> {THIS_RUN} test differences of every method, the baselines included, mix any change of the '
             'methods with the corrected training data and changed test targets'
             + (f'; v3 -> {THIS_RUN} measures no model change.\n' if not dev3 else '.\n'))
    cvs = {r: prior[r].get('cv_headline') for r in have if prior[r].get('cv_headline') is not None}
    if cvs and len(head):
        cvs[THIS_RUN] = head
        rows = []
        for sch in SCHEMES:
            for m in BASE_METHODS:
                sel = {r: t[(t['scheme'] == sch) & (t['method'] == m)] for r, t in cvs.items()}
                if any(len(d) for d in sel.values()):
                    rows.append(dict(scheme=sch, method=m, **{f'{r}_mpsp_mape': g(sel[r], 'mpsp_mape')
                                                             for r in allr if r in sel}))
        T = pd.DataFrame(rows)
        later = [r for r in allr if r != 'v1']
        L.append('**Development CV** (fold-mean MPSP MAPE, %; v1 used the v1 protocol and settings, and the '
                 f'development rows of {_and(later)} differ by the corrected data (section 1; v3 -> {THIS_RUN}: '
                 'the same samples, the values of the regenerated configurations corrected), so the CV estimates '
                 'changed as well; each version\'s changes were designed and chosen on these development folds, '
                 'so its numbers are post-selection estimates'
                 + (f' ({THIS_RUN} changed only the data and is as post-selection as v3)' if not dev3 else '')
                 + ' and the differences, above all LOFO and LOPO with only 4 and 7 held-out '
                 f'groups, are likely optimistic -- see "{CHANGES_TITLE}"):\n')
        L.append(md_table(T, pct_cols=tuple(c for c in T if c.endswith('_mpsp_mape'))) + '\n')
    return '\n'.join(L)


def settings_vs_v3(ctx, plan=None):
    """CLI flags by which this run's model, protocol and HPO settings differ
    from those of the v3 run (the pipeline defaults since round 5); [] if
    none."""
    mk, lk = ctx.get('mk_extra') or {}, ctx.get('lk_extra') or {}
    down = bool((ctx.get('spec_kw') or {}).get('downstream'))
    return [flag for off, flag in (
        (not down, '--no-downstream'),
        (down and not mk.get('water_gate'), '--no-water-gate'),
        (mk.get('channel_inputs') is None, '--no-masking'),
        (mk.get('channel_inputs') == 'physics', '--v1-g-inputs'),
        (not mk.get('r_structure'), '--no-r-structure'),
        (down and not mk.get('fs_prior'), '--no-extrapolation-safeguards'),
        (down and not mk.get('power_gate'), '--no-power-gate'),
        (not mk.get('clip_exempt'), '--no-clip-exempt'),
        (bool(mk.get('fs_b_only')), '--fs-b-only'),
        (lk.get('ch_space') != 'surface', '--v2-channel-loss'),
        (not ctx.get('target_mask'), '--no-target-mask'),
        ((ctx.get('protocol') or PROTOCOL)['name'] != 'v2', '--v1-protocol'),
        (float((plan or {}).get('hpo_lopo_weight', 0.) or 0.) > 0, '--hpo-lopo-weight')) if off]


def changelog(ctx, plan=None, info=None):
    """Markdown section: what changed between the runs that used the frozen
    test split (ALL_RUNS) and why (development evidence only; details in
    DESIGN.md section 1 and section 6, rounds 4 and 5). ``info``: see
    :func:`changes_v3_v4`."""
    L = [f'## {CHANGES_TITLE} (development-set evidence only)\n',
         'Each change was chosen on the development set (paired over seeds or draws, test split never used); the '
         'data corrections (v2 -> v3, v3 -> ' + THIS_RUN + ') fix simulation defects and were verified on '
         'development rows only. DESIGN.md section 1 (data corrections, known data defects) and section 6 '
         '("Development-set diagnostics", rounds 4 and 5) have the full evidence. Each model change has a CLI '
         'flag that restores the behaviour of the previous version. "Seed t" below is a t over training seeds on '
         'fixed folds: it measures seed noise only, not whether a new configuration, feedstock or product '
         'benefits; the fold-level statistics (seeds averaged per fold, paired over folds) answer that.\n']
    return ('\n'.join(L) + '\n' + changes_since_v1(ctx) + '\n' + changes_v2_v3(ctx, plan) + '\n'
            + changes_v3_v4(ctx, plan, info))


#: report text (changelog v3 -> v4): the effect of the D16 / D17
#: corrections on the data (development samples only; DESIGN.md section 1,
#: "Known data defects" D16 / D17 and "Data set")
D16_D17_EFFECTS = (
    '1. **D17, the anaerobic IC reactor of the wastewater treatment counted its biogas twice** (BioSTEAM\'s '
    '`InternalCirculationRx`, default lumped method: it reacted a helper copy of the influent and degassed it '
    'into the biogas, then degassed the reacted waste and effluent into it again; every HP and TAL '
    'configuration, ethanol_sugarcane, '
    'IBO_EtOH_corn and oxalic_neutral_sugarcane). The boiler burned methane that was never produced, so the '
    'steam / power credit was too high and the MPSP too low, most where the wastewater load is high (corn, '
    'corn stover, low yields). Harness v4 applies the runtime correction `ic_biogas` '
    '(`data_generation/biosteam_fixes.py`; BioSTEAM itself is unchanged) in every campaign process and records '
    'it in every sample (`biosteam_fixes`); older samples of these configurations are rejected '
    '(`unpatchable_D17`). Verified: the reactor\'s mass balance closes (HP_neutral_cornstover baseline: -45 '
    'kg/h on an 829 t/h inlet, before +18.5 t/h) and the biogas halves. Regenerated vs. v3 rows at identical '
    'TY points (development samples; median per configuration): HP glucose +0.7-0.9 %, HP sugarcane +1.0-1.2 %, '
    'HP corn +5.4-6.4 %, HP corn stover +7.5-10.0 % (5-95 %: 2.1-15.7), TAL_SA_glucose +0.9 %, TAL_glucose '
    '+2.1 %, TAL_SA_corn +3.6 %, TAL_corn +7.9 %, TAL_SA_cornstover +6.7 %, TAL_cornstover +12.2 %, '
    'ethanol_sugarcane +3.4 % (up to 18 % at the lowest yields), IBO_EtOH_corn +0.57 %, oxalic < 0.01 %; b′ '
    '(median over samples) +9.5-22 % for HP and TAL on corn and corn stover and +11 % for ethanol_sugarcane '
    '(the phantom gas grew as the yield fell), within ±2 % on glucose and sugarcane.',
    '2. **D16, TAL_SA_sugarcane fed the juice recycle to the cane mill twice** (a duplicate `ins.append` in '
    'the Bioindustrial-Park system file; the six other TAL systems comment it out): its fermentation received '
    '1.088 x the sugar the cane holds, now 0.946 as the other sugarcane configurations '
    '(`tables/sugar_balance_check.csv`). Fixed on the fi-TAL branch (`data_generation/bip_patches/fi-TAL.patch`); '
    'TAL adapter v3 for this configuration, its maximum-titer line recalibrated on the corrected system (679 y '
    '- 1.3 g/L, before 614 y - 0.6); older rows rejected (`unpatchable_D16`). With D17, the regenerated '
    'samples\' fitted MPSP at the baseline point is +4.0 % (median; 2.7-5.2 %) and b′ +11.9 % above the v3 '
    'rows.',
    '3. **Unchanged**: the 14 configurations without an IC reactor (lactic, succinic, ethanol corn and corn '
    'stover; identical rows), the seeds, sample plans and TY designs, the quality-filter rejections (the same '
    '23 lactic samples) and the frozen split (the same sample ids). Isobutanol: 14 of its development samples '
    'changed their feedstock composition slightly (up to 0.027 mass fraction; the history-dependent starch '
    'setter, open defect D10).')


def changes_v3_v4(ctx, plan=None, info=None):
    """Markdown subsection: what changed between v3 and THIS_RUN -- the
    training data only (the D16 / D17 simulation corrections of the final
    review of v3); the model, protocol and search settings are compared with
    v3's. ``info``: dict(regen=``report_checks.regenerated_rows``, n_samples, n_points,
    n_dev, n_dev_configs, n_configs, fit_mape_median, data_hash, selected (this run),
    prior (v3 tables of :func:`load_run_tables`))."""
    import numpy as np
    info = info or {}
    regen, v3 = info.get('regen'), (info.get('prior') or {}).get('selected_hyperparameters') or {}
    dev3 = settings_vs_v3(ctx, plan)
    fmt = lambda v: ('' if v is None else (f'{v:,}'.replace(',', ' ') if v >= 10000 else str(v))
                     if isinstance(v, (int, np.integer)) else str(v))
    L = [f'### v3 -> {THIS_RUN} (data corrections D16 / D17 only)\n',
         '**Data** (simulation defects found by the final review of v3; DESIGN.md section 1, "Data corrections", '
         '"Known data defects" D16 / D17 and "Data set"). Both were fixed at the source and the affected '
         'configurations regenerated with the seeds, sample plans and TY designs of `run_all.sh`; the corrected '
         'rows replace the defective ones in training and in the test split. The simulation checks that found '
         'them and define the affected configurations (`data_generation/sim_checks.py`, '
         '`data_generation/checks/`) re-simulate every configuration, the test configurations included, at '
         'baseline parameters (not their recorded rows), and the correction applies by configuration to every '
         'configuration with the defect; every number below is from development samples.\n']
    L += list(D16_D17_EFFECTS)
    if regen:
        L.append(f'4. **Rows**: {fmt(regen["n"])} samples of {regen["n_configs"]} configurations regenerated '
                 f'(harness v4): {fmt(regen["dev"])} of the {fmt(regen["n_dev"])} development samples, '
                 f'{_regen_test_clause(regen)}. Data hash {V3_RUN["data_hash"]} -> '
                 f'{info.get("data_hash", "?")}; {fmt(V3_RUN["n_samples"])} -> {fmt(info.get("n_samples"))} '
                 f'samples, {fmt(V3_RUN["n_points"])} -> {fmt(info.get("n_points"))} feasible points; development '
                 f'set {fmt(V3_RUN["n_dev"])} -> {fmt(info.get("n_dev"))} samples ({V3_RUN["n_dev_configs"]} -> '
                 f'{info.get("n_dev_configs")} configurations); median fit MAPE of the samples '
                 f'{100 * V3_RUN["fit_mape_median"]:.2f} -> '
                 + (f'{100 * info["fit_mape_median"]:.2f} %.' if info.get('fit_mape_median') is not None
                    else '? %.'))
    L.append('')
    # model / protocol / search
    p3 = v3.get('plan') or {}
    pdiff = [f'`{k}` {p3.get(k)} -> {(plan or {}).get(k)}' for k in PLAN_KEYS
             if p3 and k in p3 and (plan or {}).get(k) != p3.get(k)]
    sel, sel3 = info.get('selected') or {}, v3.get('selected') or {}
    sdiff = [f'{k} {sel3.get(k)} -> {sel.get(k)}' for k in sorted(set(sel) | set(sel3)) if sel.get(k) != sel3.get(k)]
    differs = bool(dev3) or CODE_VERSION != V3_RUN['code_version'] or bool(pdiff)
    same = not differs and bool(p3)
    L.append('**Model, training protocol and hyperparameter search: '
             + ('NOT identical to v3**' if differs else 'unchanged**' if same else
                'settings unchanged, plan not compared**') + ' (round-5 settings, code version '
             f'{CODE_VERSION}' + ('' if CODE_VERSION == V3_RUN['code_version'] else
                                   f' (v3: {V3_RUN["code_version"]})')
             + (f'; this run differs from v3 by {", ".join(f"`{f}`" for f in dev3)}' if dev3 else '')
             + ('; the same seed, deterministic folds, search candidates and criterion' if same else '')
             + (f'; plan differences from v3: {"; ".join(pdiff)}' if pdiff else
                ('' if p3 else '; v3 plan not found (`results/v3/tables/selected_hyperparameters.json`)'))
             + '). The HPO was re-run on the corrected development rows and '
             + ((('selected the same hyperparameters as v3.' if not sdiff else
                  'selected different hyperparameters from v3: ' + '; '.join(sdiff) + '.')) if sel3 else
                'its selection is in section 4 (v3 selection not found).') + '\n')
    L.append(f'**Reading the v3 -> {THIS_RUN} differences.** The targets themselves changed for '
             + (f'{regen["n_configs"]} of the {info.get("n_configs") or "?"}' if regen else 'the regenerated')
             + ' configurations (a higher '
             'MPSP, a larger b′ on corn and corn stover), in training and in every held-out fold and test set, so '
             f'the v3 -> {THIS_RUN} changes of the CV and test errors of every method, the baselines included, '
             'measure how well the methods learn the corrected data, not a model improvement (the MAPE '
             f'denominators changed as well). The v3 checkpoints (`results/v3`) learned the D16 / D17 bias (MPSP '
             f'about 1-13 % low at the baseline point for HP and TAL); {THIS_RUN} replaces them whatever the '
             'comparison shows.\n')
    return '\n'.join(L)


def changes_since_v1(ctx):
    """Markdown subsection: what changed between v1 and v2 and why
    (development evidence only; details in DESIGN.md section 6, round 4)."""
    mk = ctx.get('mk_extra') or {}
    v2p = (ctx.get('protocol') or PROTOCOL)['name'] == 'v2'
    on = lambda b, flag: '' if b else f' **Not used in this run (`{flag}`).**'
    L = ['### v1 -> v2 (round 4: training protocol, recovery inputs, structures)\n',
         '1. **Training protocol** (`PROTOCOL`): every network fit -- HPO, all CV schemes, `final_dev` and '
         '`release` -- trains on all its samples with sample-level early stopping (each member on its own '
         'random 15 % hold-out), without a calibration hold-out and without refit, so the CV estimates the '
         'training procedure of the final models with the selected hyperparameters (the HPO is not nested in '
         'the CV). v1 held out 15 % of the configurations for calibration and '
         'early-stopped every member on about five held-out configurations: noisy and often far too short '
         '(dev proxy of the frozen test: best epochs 6-226, in 4 of 6 draws a member stopped at <= 20 '
         'epochs). On that proxy the v1 protocol gave OOD 12.0 ± 6.9 % / ID 13.9 ± 5.5 %, the new one '
         '7.2 ± 1.2 % / 3.9 ± 0.2 % (6 draws, paired t = -1.8 / -4.6; the draws share configurations: 24 '
         'picks, 16 distinct). Refitting after early stopping did not '
         'help (+0.2 pp at twice the cost) and is not used.' + on(v2p, '--v1-protocol'),
         '2. **Intervals**: the final and release ensembles have no calibration hold-out; their split-conformal '
         'quantiles come from the out-of-fold scores of the development leave-one-configuration-out CV of the '
         'same procedure and ensemble size (`train.calibrate_conformal_oof`, re-scored with the final '
         'ensemble\'s spread floors, `coverage_guaranteed = False`), with a cross-fitted coverage estimate. v1 '
         'calibrated `final_dev` on its own 6 calibration configurations and transferred the release intervals '
         'from the dev model scored on its calibration configurations AND the test configurations; v2 never uses '
         'test samples for calibration.' + on(v2p, '--v1-protocol'),
         "3. **Hard productivity structure** (`ModelConfig.r_structure='physics'`): the coefficients of the "
         "fermentation capital and of the electricity intensity are theta0(x) + theta1(x) r_ref/r with "
         'theta1 >= 0, so dMPSP/dr <= 0 holds by construction (v1 violated it for 0.8-8 % of held-out samples). '
         'In the data only the fermentation capital depends materially on r, as 1/r (beta = 1.0, R² = 0.9999 '
         'within configurations). Dev productivity CV with sample-level hold-outs: 5.16 -> 4.29 % (seed t = '
         '-4.5; low-r side 6.34 -> 4.39 %, seed t = -12.5; per configuration of the single productivity fold 28 '
         'of 36 better, t = -2.9); LOCO and random CV unchanged.'
         + on(bool(mk.get('r_structure')), '--no-r-structure'),
         "4. **Recovery inputs** (masking preset 'physics_v2'): the recovery g no longer sees the reagent "
         '(base / acid) prices, which set what the reagents cost, not how much product is recovered '
         '(dg/dprice = 0 exactly; within-configuration correlation -0.05 / -0.02). Adopted as a hard physical '
         'constraint, not for accuracy: dev CV in-family g unchanged, LOPO -10.5 points (seed t = -1.2; g error '
         'seed t = -2.9; per held-out product 3 better / 4 worse, median +0.1), LOFO configuration-balanced -2.5 '
         'points (seed t = -2.4; per held-out feedstock 2 better / 2 worse, median -1.0, t over folds -0.7), '
         'LOCO neutral. Feedstock, fermentation, scale or titer-domain inputs for g, and '
         'a larger g loss weight, did not help and are not used.' + on(mk.get('channel_inputs') == 'physics_v2', '--v1-g-inputs'),
         '5. **Extrapolation safeguards** (`EXTRAPOLATION_SAFEGUARDS`): a sugar-availability prior of the '
         'feedstock intensity (free sugars and starch fully available, structural polysaccharides by the '
         'saccharification yield; the network predicts only the residual), clipping of the feedstock and '
         'pretreatment inputs to their training range, and a neutralization gate (no neutralization cost '
         'without an acid group or neutralizing agent; verified on the data). v1 extrapolated a new '
         "feedstock's minor composition fractions (corn: 9.5 training sd) into a saturated feedstock efficiency "
         '(corn 1/eta 1.50 vs 0.93-1.01). Dev LOFO fold mean 63.4 -> 34.4 % (seed t = -5.3, 5 seeds; corn 128 -> '
         '43 %); over the 4 held-out feedstocks 3 better / 1 worse, median -15.8, t over folds -1.5, corn 72 % of '
         'the gain: not established at the fold level, but the corrected feedstock channel is (DESIGN.md 4d: '
         'feedstock-channel error 23 -> 3.5 % of the numerator). LOCO and random not worse.' + on(bool(mk.get('fs_prior')), '--no-extrapolation-safeguards'),
         '6. **Plan**: the CV, final_dev and release ensembles have the same number of members (the scores '
         'divide by the member spread); the PI-NN and the plain MLP run on every LOCO fold, the PI-NN '
         'ablations on a stratified subset of configurations (`--ablation-loco-folds`). The hyperparameter '
         'search was re-run under the new protocol.',
         '',
         '**All changes together** (development set only; DESIGN.md round 4e): PI-NN with the v1 settings vs v2 on '
         "the pipeline's CV folds, 3 members, paired over 3 seeds, the v1-selected hyperparameters for both "
         '(fold-mean MPSP MAPE, %). "Folds" columns: seeds averaged per fold, then paired over the folds '
         '(held-out configurations for the single productivity fold):\n',
         '| scheme | v1 -> v2 | seed t | folds better / worse | median diff | t over folds |',
         '|---|---|---|---|---|---|',
         '| random (5 folds) | 4.10 -> 3.91 | -9.1 | 5 / 0 | -0.19 | -5.7 |',
         '| leave-one-configuration-out (37) | 16.1 -> 9.7 | -10.8 | 35 / 2 | -3.7 | -4.4 |',
         '| ... HP / TAL / succinic / lactic configurations (32; matched OOD estimate) | 13.3 -> 6.2 | -9.5 | 31 / 1 | '
         '-3.3 | -4.9 |',
         '| leave-one-feedstock-out (4) | 64.1 -> 34.9 | -24 | 3 / 1 (cornstover +8) | -14.9 | -1.3 |',
         '| leave-one-product-out (7) | 86.1 -> 50.1 | -7.5 | 6 / 1 (succinic +23) | -5.6 | -1.5 |',
         '| productivity (1 fold; 36 configurations) | 18.3 -> 4.3 | -4.2 | 36 / 0 | -8.0 | -9.2 |',
         '',
         'The LOCO and productivity gains hold at the fold (configuration) level. The leave-one-feedstock-out and '
         'leave-one-product-out gains are **not statistically established at the fold level**: one fold '
         'dominates each (corn 72 %, isobutanol 48 % of the summed absolute change). On the development proxy of '
         'the frozen test (6 draws that share configurations: 24 picks, 16 distinct): OOD 11.8 -> 7.1 %, ID '
         '14.3 -> 3.8 % (6 of 6 draws better), better than HGB (9.2 / 4.8 %), the plain MLP (10.3 / 6.1 %) and '
         'ridge (11.4 / 7.7 %) on both.\n',
         '**Post-selection caveat.** These v2 figures, and the CV of the v2 and later reports (same folds: LOCO, '
         'LOFO, LOPO and productivity are deterministic, and the random CV with the default `--seed 0` too), come '
         'from the folds and proxy draws on which round 4 '
         'was diagnosed and its variants were selected, so they are optimistic -- most of all leave-one-feedstock-out, '
         'whose safeguards were designed from the corn failure. Confirmation on 6 fresh proxy draws and a fresh '
         'random-CV fold seed that played no role in any decision (same configurations, new combinations, ID '
         'samples and seeds; DESIGN.md 4e): OOD 13.3 -> 6.9 %, ID 13.9 -> 4.0 % (6 of 6 draws better), v2 better '
         'than HGB (8.5 / 4.9 %), the plain MLP (8.9 / 5.7 %) and ridge (9.8 / 8.0 %) on both; random CV 4.12 -> '
         '3.95 % (4 of 5 folds better). No fresh LOFO or LOPO folds exist.\n',
         '']
    return '\n'.join(L)


def changes_v2_v3(ctx, plan=None):
    """Markdown subsection: what changed between v2 and v3 and why -- the
    data corrections of the final review of v2 and the round-5 model changes
    (development evidence only; DESIGN.md section 1 and section 6, round 5).
    The "Not used in this run" marks refer to the run that writes the report
    (THIS_RUN keeps the round-5 settings)."""
    mk = ctx.get('mk_extra') or {}
    lk = ctx.get('lk_extra') or {}
    tm = ctx.get('target_mask') or {}
    w_lopo = float((plan or {}).get('hpo_lopo_weight', 0.) or 0.)
    on = lambda b, flag: '' if b else f' **Not used in this run (`{flag}`).**'
    masked = '; '.join(f"{k}: {', '.join(v)}" for k, v in tm.items()) or 'none'
    L = ['### v2 -> v3 (data corrections and round 5)\n',
         '**Data** (defects found by the final review of v2; DESIGN.md section 1, "Data corrections" and "Known '
         'data defects"). The corrected rows replace the defective ones in training and in the test split: '
         '1550 -> 1600 samples (TAL 185 -> 232, succinic 141 -> 144), development set 1306 -> 1353.\n',
         '1. **D15, succinic acid at neutral pH** (4 configurations, 72 samples): succinic adapter v2 set only the '
         'neutralization flag, and the reactor kept the purchased lime out of the broth, so the salt formation, '
         'the acidulation (R401) and the gypsum filtration (S405) were never simulated: the v2 rows equalled the '
         'low-pH process except the base cost, while the downstream specification (the model inputs) described '
         'acidulation and gypsum filtration. Regenerated with adapter v3 (72 / 72 samples, no failed point): MPSP '
         '+0.22 USD/kg vs. the low-pH process at identical TY points and parameters (+0.10 vs. the v2 rows): lime '
         '1 mol/mol succinic acid, 0.89-0.90 kg H2SO4 (93 wt%) and 1.44-1.45 kg wet gypsum per kg; recovery '
         'unchanged. Rows of adapter < 3 are rejected (`unpatchable_D15`).',
         '2. **H1, TAL_SA_glucose / TAL_SA_sugarcane**: after a convergence failure whose simulator reset also '
         'failed, the old harness wrote the later samples as complete without a simulated point (46 samples, and '
         'one cut after 3 points); the quality filter rejected them, so v2 trained on 1 and 10 of 29 samples. '
         'Cause: the IPA evaporator F407 (a root search that evaporates all the IPA in one effect); the TAL adapter '
         'now returns the exact overall evaporation fraction in that region (re-simulated recorded points change '
         'by < 1e-4 in MPSP), the harness reloads after a failed reset, and the 47 samples were re-run (29 / 29 '
         'each).',
         '3. **D8, TAL yield basis**: analysed, y kept (the module\'s yield on all substrates, glucose-eq. by '
         'theoretical yield). Its consequence, a feedstock-channel decomposition with a negative titer term from '
         'the spiked acetate\'s sugar substitution (c\' < 0 in every TAL sample; no input describes the acetate), '
         'is handled by model change 3.',
         '',
         '**Model** (round 5; DESIGN.md section 6, round 5: designed on a frozen copy of the v2 tables, '
         'development rows only, one-member fits paired over 3 seeds; re-checked on the corrected development '
         'rows, table below). Each has a CLI flag that restores v2 (all together `--v2-model`):\n',
         "1. **Channel loss in surface space with a fixed weight** (`LossConfig.ch_space='surface'`, `w_ch = 3`, "
         'no longer an HPO dimension): v2 compared each channel\'s four coefficients in units of their '
         'cross-configuration spread, at the HPO\'s weight 0.1, which left the decomposition free once MPSP fit '
         '(in-sample: every corn configuration\'s feedstock sensitivity 19-27 % low, the isobutanol electricity '
         'response off by 94 % of MPSP). The surface term compares every channel numerator (price channels as '
         'price x intensity) at the sample\'s simulated TY points relative to the total numerator, which is the '
         'out-of-fold fidelity metric of section 5. `w_ch` was chosen on out-of-fold fidelity (1 < 3 ~ 10 on LOCO; '
         '3 best on LOFO and LOPO), not by the MPSP criterion, which cannot see the decomposition.'
         + on(lk.get('ch_space') == 'surface', '--v2-channel-loss'),
         '2. **Power gate** (`ModelConfig.power_gate`): where electricity is not among the coproducts (no '
         'turbogenerator: corn ethanol, the lactic configurations), every coefficient of the net electricity '
         'intensity is mapped to s softplus(q / s) >= 0, so such a plant cannot be predicted to export power '
         '(price linearity and dMPSP/dr <= 0 kept). A hard physical constraint, checked on the data (the net '
         'electricity is never negative at the 3667 points of those configurations; the gate switches itself off '
         'if the data violate it). It removes a qualitative failure (out of fold, v2 predicted export for 85-100 % '
         'of the held-out samples of those configurations), at a cost for corn ethanol: its v2 prediction offset a '
         'missed DDGS / corn-oil credit with a wrong-sign electricity credit, which the gate removes (open issue).'
         + on(bool(mk.get('power_gate')), '--no-power-gate'),
         f'3. **Known defective channel targets masked in the loss** (`experiments.DEFECTIVE_TARGETS`; this run: '
         f'{masked}): the TAL feedstock channel (D8) carries the acetate\'s sugar substitution, which no input '
         'describes and the product-blind feedstock head cannot represent; its totals, g and other channels keep '
         'their targets. A correction of defective targets, neutral on the development CV. The study also '
         'proposed masking lactic acid\'s feedstock targets (a recovery shape, not a defect); on the corrected '
         'data that helped only the corn leave-one-feedstock-out fold in 2 of 3 seeds and cost leave-one-product-out '
         'about 2 points, so it is not used.' + on(bool(tm), '--no-target-mask'),
         '4. **The exact final-product stoichiometric ratio is never clipped** (`ModelConfig.clip_exempt`): the '
         'v2 HPO chose `clip_inputs = 0`, which clipped the held-out HP\'s ratio (0.80, acrylic acid) to the '
         'training minimum 1.00, a combination with its dehydration upgrade that no training product has: the '
         'cause of the v2 leave-one-product-out HP failure (development LOPO, v2 hyperparameters, 2 members, 3 '
         'seeds: HP 82.0 -> 14.8 %; every other LOPO fold bit-identical; LOCO unaffected).'
         + on(bool(mk.get('clip_exempt')), '--no-clip-exempt'),
         '5. **HPO criterion unchanged** (leave-configurations-out only'
         + (f'; **this run: leave-one-product-out weight {w_lopo:g}**, `--hpo-lopo-weight`' if w_lopo > 0 else '')
         + '): the study proposed adding the leave-one-product-out folds with weight 0.5; with the clip exemption '
         'that criterion still switches the selection, to a candidate whose single-seed LOPO advantage (46.1 vs '
         '49.7 %) reverses over three seeds (54.5 ± 8.0 vs 49.5 %), at +1.3 points of the LOCO criterion, +21 % '
         'run time and a post-selection LOPO estimate.',
         '6. **Reporting**: out-of-fold channel and price-intensity fidelity per scheme and held-out '
         'configuration (section 5).',
         '',
         "Not used: a pure-b' feedstock intensity (`--fs-b-only`: LOPO channel error +4 points), `w_ch` 1 or 10, "
         'wider input-clipping margins (they do not resolve the HP / oxalic trade-off).\n',
         ROUND5_RECHECK,
         '**Post-selection caveat.** The round-5 options were designed and selected on the same deterministic '
         'development folds as this report\'s CV (first on the v2 rows, then re-checked on the corrected rows), '
         'so the CV estimates of v3 and later runs are post-selection like v2\'s; the v2 -> v3 differences of '
         'the CV tables also include the data corrections (the development rows changed). Fresh folds do not '
         'exist for LOFO and LOPO.\n',
         '']
    return '\n'.join(L)


def intervals_section(fin, protocol):
    """Markdown: how the final intervals were calibrated and their realized
    test coverage (report only)."""
    import numpy as np
    L = []
    r = (fin.get('final') or {}).get('pinn') or {}
    conf = r.get('conformal') or {}
    if not conf:
        return 'PI-NN intervals: none (no usable calibration scores; see log).\n'
    L.append(f'**Intervals.** Source: {conf.get("source")}. {conf.get("n")} out-of-fold samples of '
             f'{conf.get("n_calib_groups")} configurations, weighting `{conf.get("weighting")}`, nominal level '
             f'{1 - conf.get("alpha", ALPHA):.2f}, `coverage_guaranteed = {conf.get("coverage_guaranteed")}`: the '
             'calibration scores come from the CV fold ensembles (same procedure, same number of members), not '
             'from the final ensemble itself, so the split-conformal guarantee does not apply exactly; the '
             'cross-fitted estimate (leave one LOCO fold out) is reported instead. Quantiles: '
             + ', '.join(f'{k} {_fmt(conf.get(k, float("nan")))}' for k in ('a', 'b', 'c', 'd', 'g', 'mpsp'))
             + '. A sample counts as covered for MPSP only if every TY point of it is inside its interval '
             '(sample-simultaneous); `coverage_pointwise` counts points.\n'
             if protocol['name'] == 'v2' or conf.get('source') != 'calibration hold-out' else
             f'**Intervals** (v1 protocol): calibrated on the final_dev calibration hold-out ({conf.get("n")} samples '
             f'of {conf.get("n_calib_groups")} configurations).\n')
    iv = fin.get('intervals')
    if iv is not None and len(iv):
        own = iv['calibration'].isin(['loco', 'holdout'])
        sub = iv[(iv['set'] != 'test_all') & (iv['method'] == 'pinn') & (own | iv['key'].isin(['mpsp', 'g']))]
        cols = ['method', 'calibration', 'set', 'key', 'n', 'q', 'coverage', 'crossfit_pooled', 'coverage_cfg_mean',
                'crossfit_mean', 'coverage_cfg_min', 'crossfit_min', 'crossfit_q10', 'coverage_pointwise',
                'halfwidth_rel_mean', 'halfwidth_rel_median']
        L.append('Realized test coverage (report only, never used for any choice) next to the cross-fitted '
                 'estimate of the calibration (`tables/final_intervals.csv`). Compare like with like: `coverage` '
                 '(pooled over the test samples) with `crossfit_pooled`, `coverage_cfg_mean` (mean over the test '
                 'configurations) with `crossfit_mean` (mean over the LOCO folds), the worst configuration with '
                 '`crossfit_min`. `calibration` = loco: the intervals of the saved models (out-of-fold LOCO scores); '
                 'holdout: own calibration hold-out (`--v1-protocol` only); random_cv: reporting-only calibration on '
                 'the random-CV out-of-fold scores, stored in the checkpoint meta as `conformal_random_cv`. The '
                 'plain MLP\'s rows are in the CSV:\n')
        L.append(md_table(sub[[c for c in cols if c in sub]]) + '\n')
    pc = fin.get('intervals_per_config')
    if pc is not None and len(pc):
        sub = pc[(pc['method'] == 'pinn') & pc['calibration'].isin(['loco', 'holdout']) & (pc['set'] == 'test_ood')
                 & pc['key'].isin(['mpsp', 'g'])]
        if len(sub):
            t = sub.pivot_table(index='config_id', columns='key', values='coverage', aggfunc='mean').reset_index()
            t.columns = [str(c) if c == 'config_id' else f'coverage_{c}' for c in t.columns]
            est = sub.groupby('key')['crossfit_mean'].first()
            L.append('PI-NN realized coverage per held-out test configuration (estimate = cross-fitted mean: '
                     + ', '.join(f'{k} {_fmt(v)}' for k, v in est.items()) + '):\n')
            L.append(md_table(t) + '\n')
    L.append('Coverage is marginal over the samples of new configurations, with the configurations weighted by '
             'their numbers of samples (weighting `sample`), not guaranteed per configuration: the '
             'per-configuration coverage of the out-of-fold scores spreads widely (cross-fitted `min` / `q10` in '
             'section 5), and a new configuration of a family whose scores are larger than the calibration set\'s '
             'is under-covered (development proxy of the frozen test with LOCO-calibrated intervals, 6 draws: '
             'MPSP coverage 0.86-0.99 per draw, mean 0.93, but TAL_glucose / TAL_corn 0.70 in 2 draws and '
             'lactic_lowpH_SHF_cornstover 0.58 in one; DESIGN.md round 4f). For a known configuration (in '
             'distribution) the intervals are conservative pooled over samples but not for every family: the '
             'paragraph "Coverage for a KNOWN configuration" below (in the v1 and v2 reports this sentence read '
             '"conservative (proxy 0.95-0.99)" without the per-family breakdown). No test sample was used for '
             'selection or calibration.\n')
    return '\n'.join(L)


#: tables of the ``unit_bounds`` simulation check (known data defect D18),
#: relative to the package directory
UNIT_BOUNDS_CSV = os.path.join('data_generation', 'checks', 'unit_bounds.csv')
COLUMN_TRAYS_CSV = os.path.join('data_generation', 'checks', 'column_trays_bias.csv')
COLUMN_TRAYS_CAPACITY_CSV = os.path.join('data_generation', 'checks', 'column_trays_capacity.csv')


def _family_undercoverage(checks, key, tol=0.1):
    """Text: the families whose coverage of the shipped ``key`` quantile on
    the LOCO out-of-fold scores is more than ``tol`` below nominal."""
    fc = checks.get('family_coverage')
    if fc is None or not len(fc):
        return ''
    f = fc[fc['key'] == key]
    low = f[f['coverage'] < 1 - ALPHA - tol]
    if not len(low):
        return ''
    cc = checks.get('config_coverage')
    worst = ''
    if cc is not None and len(cc):
        c = cc[cc['key'] == key].sort_values('coverage').head(2)
        worst = '; worst configurations ' + ', '.join(f'{x.config_id} {x.coverage:.2f}' for x in c.itertuples())
    return ('; but per family the shipped quantile UNDER-covers the new configurations of '
            + ', '.join(f'{x.family} ({x.coverage:.2f})' for x in low.itertuples())
            + f' (LOCO out-of-fold scores{worst}; `predict()` / `mpsp()` flag such rows in '
            '`intervals_undercovered`, section 6)')


def _column_trays_note(rise_cfgs):
    """Text: the mechanism of the super-linear capital of the configurations
    ``rise_cfgs`` from the ``unit_bounds`` simulation check (D18), or the
    plain statement when its tables are missing."""
    import pandas as pd
    try:
        ub = pd.read_csv(os.path.join(PKG_DIR, UNIT_BOUNDS_CSV))
        cb = pd.read_csv(os.path.join(PKG_DIR, COLUMN_TRAYS_CSV))
        cc = pd.read_csv(os.path.join(PKG_DIR, COLUMN_TRAYS_CAPACITY_CSV))
    except (OSError, ValueError):
        return ': cost correlations used beyond their size ranges, not investigated'
    wide = set(cb.loc[cb['affected'].astype(str) == 'True', 'config_id'].astype(str))
    if not wide:
        return ''
    d = ub[ub['key'].astype(str).str.lower().str.endswith('diameter') & (ub['side'] == 'above')]
    dmax = d.groupby('config_id')['value'].max()
    cov = [c for c in rise_cfgs if c in wide]
    out = (f'; known data defect D18, DESIGN.md section 1: the purification columns of all {len(wide)} '
           f'HP configurations are wider than BioSTEAM\'s own 24 ft diameter bound already at the baseline '
           f'({dmax.min():.1f}-{dmax.max():.1f} ft, `data_generation/checks/unit_bounds.csv`), and BioSTEAM prices '
           'their trays with an exponential correlation in the diameter (fitted to about 16 ft), so their '
           f'capital grows super-linearly with the flow; {len(cov)} of the {len(rise_cfgs)} configurations above '
           'have such a column')
    a = cc[(cc['mode'] == 'as_is')].set_index(['config_id', 'mult'])
    f = cc[(cc['mode'] == 'fixed')].set_index(['config_id', 'mult'])
    if (len(a) and (1.15 in a.index.get_level_values('mult')) and (2. in a.index.get_level_values('mult'))):
        fci = a.xs(1.15, level='mult')['FCI_rel_to_1x']
        worst = str(fci.idxmax())
        m15, m2 = a.loc[(worst, 1.15), 'MPSP_rel_to_1x'], a.loc[(worst, 2.), 'MPSP_rel_to_1x']
        out += (f'. The effect starts INSIDE the sampled capacity range (0.81-1.19x): {worst} FCI '
                f'x{fci.max():.2f} and MPSP {100 * (m15 - 1):+.1f} % at 1.15x, MPSP {100 * (m2 - 1):+.0f} % at 2x')
        if (worst, 2.) in f.index and (worst, 1.15) in f.index:
            out += (f'; with the trays of a column wider than 24 ft priced as parallel stacks of the same area (a '
                    f'sensitivity, `data_generation/checks/column_trays_capacity.csv`) {100 * (f.loc[(worst, 1.15), "MPSP_rel_to_1x"] - 1):+.1f} '
                    f'% at 1.15x and {100 * (f.loc[(worst, 2.), "MPSP_rel_to_1x"] - 1):+.0f} % at 2x, a normal '
                    'economy of scale')
    return out


def coverage_reading(c, fam_cov=None, cfg_cov=None, nominal=None, tol=0.02):
    """Data-driven reading of the cross-fitted LOCO coverage per key
    (pooled, mean and 10 % quantile over configurations vs. nominal) and,
    with the per-family / per-configuration coverage of the shipped
    quantiles (``report_checks.family_coverage``), where it fails."""
    nominal = 1 - ALPHA if nominal is None else nominal
    L = []
    for x in c.itertuples():
        st = ('below nominal' if x.pooled < nominal - tol else 'above nominal' if x.pooled > nominal + tol
              else 'at nominal')
        if x.q10 < nominal - 0.2:
            L.append(f'{x.key}: pooled {x.pooled:.2f} ({st}), but the 10 % quantile over configurations is '
                     f'{x.q10:.2f}')
        elif st != 'at nominal':
            L.append(f'{x.key}: pooled {x.pooled:.2f} ({st})')
    out = ('Pooled over samples every key is within ±0.02 of nominal' if not L else
           'Departures from nominal: ' + '; '.join(L)) + '.'
    if fam_cov is not None and len(fam_cov):
        parts = []
        for k in ('g', 'mpsp'):
            f = fam_cov[fam_cov['key'] == k]
            if not len(f):
                continue
            seg = f'{k}: ' + ', '.join(f'{x.family} {x.coverage:.2f}' for x in f.itertuples())
            cc = cfg_cov[cfg_cov['key'] == k] if cfg_cov is not None and len(cfg_cov) else None
            if cc is not None and len(cc):
                seg += (f' ({int((cc["coverage"] >= 0.999).sum())} of {len(cc)} configurations covered at 1.0, '
                        f'{int((cc["coverage"] < 0.7).sum())} below 0.7)')
            parts.append(seg)
        low = fam_cov[fam_cov['key'].isin(('a', 'b', 'c', 'd')) & (fam_cov['coverage'] < nominal - 0.1)]
        if len(low):
            parts.append('the coefficient intervals more than 0.1 below nominal: ' + ', '.join(
                f"{k}' " + ', '.join(f'{x.family} {x.coverage:.2f}' for x in g.itertuples())
                for k, g in low.groupby('key')))
        if parts:
            out += (' Per family, coverage of the shipped quantiles on the LOCO out-of-fold scores (pooled over '
                    'the family\'s held-out samples; `tables/interval_coverage_by_family.csv`, per configuration '
                    '`tables/interval_coverage_loco_per_config.csv`): ' + '; '.join(parts)
                    + '. A pooled coverage near nominal can average strongly over- and under-covered families: '
                    'a new configuration of an over-covered family gets intervals that are too wide, one of an '
                    'under-covered family intervals that are too narrow.')
    return out


def in_distribution_coverage(shc, idf, idc=None, nominal=None, tol=0.1):
    """(summary, markdown) of the coverage of the shipped (release)
    quantiles for a KNOWN configuration: pooled on the random-CV and
    productivity-tail out-of-fold scores (``shc``,
    ``report_checks.shipped_quantile_coverage``) and the families whose
    coverage is more than ``tol`` below nominal (``idf``, per family;
    ``idc``, per configuration)."""
    import numpy as np
    nominal = 1 - ALPHA if nominal is None else nominal
    if shc is None or not len(shc) or idf is None or not len(idf):
        return '', ''
    lab = {'random': 'random CV (new parameter draws)', 'productivity': 'productivity tails'}
    pooled = []
    for sch in ('random', 'productivity'):
        t = shc[shc['scheme'] == sch].set_index('key')
        if len(t):
            pooled.append(f'{lab[sch]} ' + ', '.join(f'{k} {t.loc[k, "pooled"]:.2f}' for k in ('mpsp', 'g', 'a', 'b',
                                                                                           'c', 'd') if k in t.index))
    low = idf[idf['coverage'] < nominal - tol]
    parts = []
    for x in low.sort_values(['scheme', 'coverage']).itertuples():
        seg = f'{x.key} {x.family} {x.coverage:.2f}'
        if idc is not None and len(idc):
            c = idc[(idc['scheme'] == x.scheme) & (idc['key'] == x.key) & (idc['family'] == x.family)]
            if x.n_configs > 1 and len(c):
                seg += f' (per configuration {c["coverage"].min():.2f}-{c["coverage"].max():.2f})'
        parts.append((x.scheme, seg))
    other = (idf[idf['coverage'] >= nominal - tol]['coverage'].min() if (idf['coverage'] >= nominal - tol).any()
             else float('nan'))
    by = {sch: [seg for s_, seg in parts if s_ == sch] for sch in ('random', 'productivity')}
    low_txt = '; '.join(f'{lab[sch]}: ' + ', '.join(v) for sch, v in by.items() if v)
    tal_g = bool(((low['family'] == 'TAL') & (low['key'] == 'g')).any())
    md = ('**Coverage for a KNOWN configuration** (the shipped `models/release` quantiles, re-scored with its '
          'spread floors, on the random-CV out-of-fold scores -- a training configuration with new parameter '
          'draws -- and on the productivity-tail scores; `report_checks.shipped_quantile_coverage`, '
          '`tables/interval_coverage_by_scheme.csv`; per family `tables/interval_coverage_in_distribution_by_family'
          '.csv`, per configuration `tables/interval_coverage_in_distribution_per_config.csv`; nominal '
          f'{nominal:.2f}). Pooled over the samples the intervals are conservative: ' + '; '.join(pooled) + '. '
          + ('They are NOT conservative for every family: ' + low_txt
             + (f'; every other family and key {other:.2f} or more' if np.isfinite(other) else '') + '.'
             if low_txt else f'Every family and key is at {other:.2f} or more.')
          + (' The TAL recovery under-coverage holds for every TAL configuration, known or new: the '
             'within-configuration spread of g comes from a parameter that is not a model input (the TAL '
             'decarboxylation conversion of the separation, `bp_separation_tal_decarboxylation_conversion`; '
             'sd(ln g) 0.11, DESIGN.md section 6, round 4), i.e. noise for the model, which the score |error| / '
             '(member spread + floor) does not represent; a per-family calibration or a score with a noise term '
             'for g is an open item.' if tal_g else '')
          + f' `predict()` / `mpsp()` list, per row, the keys whose coverage for the row\'s fermentation product is '
          f'more than {tol:.1f} below nominal in `intervals_undercovered`, with an ApplicabilityWarning (the '
          'coverage of each checkpoint\'s own quantiles on these development out-of-fold scores, '
          '`meta[\'coverage_by_family\']`: for a training `config_id` the lower of the random-CV and '
          'productivity-tail coverage, for any other row of a known product the LOCO coverage; report-only '
          'metadata, the quantiles are unchanged).\n')
    summ = ('* For a KNOWN configuration the shipped intervals are conservative pooled ('
            + '; '.join(pooled) + ')' + (', but not for every family: ' + low_txt if low_txt else '')
            + ' (section 6).') if pooled else ''
    return summ, md


def interval_extrapolation_section(shc, cfx, fam_cov=None, cfg_cov=None):
    """Markdown: per-key cross-fitted coverage (new configurations) and the
    coverage of the shipped quantiles on the leave-one-feedstock-out /
    leave-one-product-out out-of-fold scores (new feedstocks / products)."""
    L = []
    if cfx is not None and len(cfx):
        c = cfx[(cfx['scheme'] == 'config') & (cfx['method'] == 'pinn') & (cfx['weighting'] == 'sample')]
        if len(c):
            L.append('**Coverage per interval key for a new configuration** (cross-fitted LOCO estimate, PI-NN, '
                     'nominal 0.90; pooled over samples = the quantity the intervals target, mean / 10 % quantile '
                     'over configurations): ' + '; '.join(
                         f'{x.key} {x.pooled:.2f} / {x.mean:.2f} / {x.q10:.2f}' for x in c.itertuples())
                     + '. ' + coverage_reading(c, fam_cov, cfg_cov) + '\n')
    if shc is not None and len(shc):
        t = shc[shc['scheme'].isin(['feedstock', 'product'])]
        if len(t):
            L.append('**The intervals are not calibrated for a new feedstock or a new product.** Coverage of the '
                     'shipped quantiles (`models/release`, calibrated on the LOCO scores) on the out-of-fold scores '
                     'of leave-one-feedstock-out (`feedstock`) and leave-one-product-out (`product`), re-scored '
                     'with the release spread floors (`report_checks.shipped_quantile_coverage`, '
                     '`tables/interval_coverage_by_scheme.csv`; nominal 0.90; `group_*` over the held-out '
                     'feedstocks / products):\n')
            L.append(md_table(t[['scheme', 'key', 'q', 'n', 'n_groups', 'pooled', 'group_mean', 'group_min',
                                 'worst_group']]) + '\n')
            L.append('`CoefficientPredictor.predict` / `mpsp` flag rows whose feedstock, fermentation product or '
                     'product pair is not in training (`intervals_calibrated = False`, with an '
                     'ApplicabilityWarning).\n')
    return '\n'.join(L)


def with_data_rows(t, dv, kind, col):
    """Table ``t`` (one row per ``col`` x method) with a ``data (fit)`` row
    after the rows of each ``col`` value: the violation rates of the data
    (``data_violation_table``); its other columns are empty."""
    import pandas as pd
    if dv is None or not len(dv) or col not in t:
        return t
    d = dv[dv['kind'] == kind].set_index('set')
    out = []
    for v, g in t.groupby(col, sort=False):
        out.append(g)
        if v in d.index:
            row = {c: d.loc[v, c] for c in t.columns if c.startswith('viol_') and c in d.columns}
            out.append(pd.DataFrame([dict({col: v, 'method': 'data (fit)'}, **row)], columns=t.columns))
    res = pd.concat(out, ignore_index=True)
    for c in t.columns:      # keep integer columns integer (the data rows leave them empty)
        if pd.api.types.is_integer_dtype(t[c]):
            res[c] = pd.Series([int(x) if pd.notna(x) else float('nan') for x in res[c]], dtype=object)
    return res


def data_violation_note(byc):
    """Which configurations make the data violate the monotonicity priors
    (``report_checks.data_violations_by_config`` over all samples), why,
    and the consequence: the residuals are enabled only where the targets
    satisfy them in >= 99 % of the samples, so they are off for every
    product."""
    if byc is None or not len(byc):
        return ''
    why = {('mono_y', 'isobutanol'): 'the MPSP rises with the yield because the yield deficit is sold as ethanol at '
                                     'a fixed price (b′ < 0 and d′ < 0 in most samples)',
           ('mono_t', 'oxalic'): 'at high titer the MPSP rises with the titer, because evaporating the juice to '
                                 'reach it cuts the exported electricity',
           ('g_le_gmax', 'oxalic'): 'the recovery slightly exceeds the stoichiometric ratio 1'}
    parts = []
    for (k, g) in byc.groupby('constraint', sort=False):
        lab = {'mono_y': 'dMPSP/dy <= 0', 'mono_t': 'dMPSP/dt <= 0', 'positive': 'MPSP > 0',
               'g_le_gmax': 'g <= g_max'}.get(k, k)
        tot = int(g['n_samples_all_configs'].iloc[0])
        cfgs = []
        for x in g.itertuples():
            fam = x.config_id.split('_')[0].lower()
            fam = {'ibo': 'isobutanol'}.get(fam, fam)
            extra = f' (max g / g_max {x.max_g_over_gmax:.4f})' if k == 'g_le_gmax' else ''
            cfgs.append(f'`{x.config_id}` {x.n_violating} of {x.n_samples} samples{extra}'
                        + (f': {why[(k, fam)]}' if (k, fam) in why else ''))
        parts.append(f'{lab}: violated by {int(g["n_violating"].sum())} of the {tot} samples, all of '
                     + '; '.join(cfgs))
    return ('**Violations in the data** (all samples; `tables/data_constraint_violations_by_config.csv`): '
            + '. '.join(parts) + '. Every other configuration satisfies these constraints. A residual is enabled '
            'only where the targets satisfy it in >= 99 % of the samples (`train.check_constraints`), so the '
            'yield- and titer-monotonicity residuals are disabled for EVERY product because of these two '
            'configurations (DESIGN.md section 1, "Product-specific cost structure"), and a method\'s '
            '`viol_mono_y` near the data\'s is not a model error. The predicted g is bounded by g_max.\n')


def capacity_response_text(t, samples=None):
    """Section 8 paragraph: the release model's MPSP response to the
    feedstock capacity vs. BioSTEAM (``report_checks.capacity_response``)."""
    import numpy as np
    import pandas as pd
    if t is None or not len(t):
        return ''
    span = ''
    if samples is not None and 'p_feedstock_capacity' in samples:
        c = samples.groupby('config_id')['p_feedstock_capacity'].agg(['min', 'max', 'median']).dropna()
        if len(c):
            span = (f' Each configuration\'s training capacity spans only {c["min"].div(c["median"]).min():.2f}-'
                    f'{c["max"].div(c["median"]).max():.2f} x its median (with every other parameter varied at '
                    'the same time), so the within-configuration scale response is weakly identified.')
    big = t[t['mult'] >= 2]
    sh = big.groupby('mult')['share_of_simulated_change'].median()
    piv = t.pivot_table(index='config_id', columns='mult', values=['dMPSP_sim_pct', 'dMPSP_pred_pct'])
    cols = [m for m in (0.5, 2., 3.) if m in piv['dMPSP_sim_pct'].columns]
    tab = pd.DataFrame({'config_id': piv.index})
    for m in cols:
        tab[f'sim {m:g}x %'] = piv[('dMPSP_sim_pct', m)].to_numpy()
        tab[f'pred {m:g}x %'] = piv[('dMPSP_pred_pct', m)].to_numpy()
    rise = big[big['dMPSP_sim_pct'] > 0]
    out_flag = big['contextual_outside_range'].astype(str).str.contains('feedstock_capacity').mean()
    ind = big['in_domain'].mean()
    return ('**(ii-c) Capacity response** (report only): MPSP change from the configuration\'s baseline capacity '
            'to 0.5x / 2x / 3x, predicted (`models/release`, queried by `config_id` with the capacity) vs. '
            'simulated (BioSTEAM, every parameter at baseline, baseline yield, titer and productivity; '
            '`data_generation/sim_checks.py`, `data_generation/checks/capacity_response.csv`; '
            '`report_checks.capacity_response`, `tables/capacity_response_release.csv`). Median share of the '
            'simulated MPSP decrease that the model reproduces (configurations whose simulated MPSP falls): '
            + ', '.join(f'{m:g}x {v:.2f} (n = {int(big.loc[big["mult"] == m, "share_of_simulated_change"].notna().sum())})'
                        for m, v in sh.items() if np.isfinite(v)) + '.' + span
            + (' The simulated MPSP RISES at 2x or 3x for ' + ', '.join(sorted(rise['config_id'].unique()))
               + ' (their simulated fixed capital investment grows faster than the capacity'
               + _column_trays_note(sorted(rise['config_id'].unique())) + ').' if len(rise) else '')
            + f' The global applicability-domain check keeps {100 * ind:.0f} % of the 2-3x queries in_domain (the '
            'global capacity range is wide because the configurations differ in scale); the per-configuration '
            f'contextual range check flags {100 * out_flag:.0f} % of them (`contextual_outside_range`, with an '
            'ApplicabilityWarning). Scale-up effects should not be read from the model:\n\n'
            + md_table(tab.round(1)) + '\n')


def loco_structure_section(checks):
    """Section 5 paragraphs on the structure of the LOCO folds
    (``clustered_and_twin_checks``): the family-level comparison with the
    baselines, LOCO without the new-product folds, and LOCO per family."""
    import numpy as np
    import pandas as pd
    L = []
    bf = checks.get('baseline_family')
    if bf is not None and len(bf):
        t = bf.copy()
        for c in ('mean', 'median', 'se', 'dominant_delta'):
            t[c] = 100 * t[c]
        pf = checks.get('loco_per_family')
        cnt = (', '.join(f'{x.n_folds} {x.family}' for x in pf.sort_values('n_folds', ascending=False).itertuples())
               if pf is not None and len(pf) else '')
        L.append(f'**The LOCO folds are not independent**: its configurations cluster in {int(t["n_families"].max())} '
                 f'biorefinery families ({cnt}) that share inputs and errors, so the '
                 'fold-level p-values above (and the ± s.e. of the LOCO fold means) treat correlated folds as '
                 'independent and overstate the precision. Family-level comparison (the paired fold differences '
                 'averaged per family, then tested over the families; `report_checks.baseline_family_stats`, '
                 '`tables/cv_baselines_paired_family.csv`; percentage points):\n')
        bk = checks.get('baseline_family_known')
        if bk is not None and len(bk):
            k = bk.copy()
            for c in ('mean', 'median', 'se', 'dominant_delta'):
                k[c] = 100 * k[c]
            k['scheme'] = 'config, known-product families'
            t = pd.concat([t, k[t.columns.intersection(k.columns)]], ignore_index=True)
        L.append(md_table(t[['scheme', 'method', 'n_folds', 'n_families', 'mean', 'se', 'median', 'n_pinn_better',
                             'n_baseline_better', 'p_wilcoxon', 'p_sign', 'p_min_attainable', 'dominant_family',
                             'dominant_delta']]) + '\n')
        if bk is not None and len(bk):
            L.append('Rows `config, known-product families`: the same without the new-product folds '
                     f'({bk["excluded_folds"].iloc[0]}; each is a family of one configuration whose product is '
                     'nowhere else in training), i.e. for the use "a new configuration of a known product". With '
                     f'{int(bk["n_families"].max())} families the smallest attainable two-sided p is '
                     f'{bk["p_min_attainable"].max():.3f} (`p_min_attainable`), so the direction (`n_pinn_better`) '
                     'and the mean / median are the evidence there; '
                     '`tables/cv_baselines_paired_family_known_product.csv`.\n')
    kp, newp = checks.get('loco_known_product'), checks.get('loco_new_product_folds') or []
    if kp is not None and len(kp) and newp:
        r = kp.set_index('method')
        pin = r.loc['pinn'] if 'pinn' in r.index else None
        parts = []
        for m in ('plain', 'knn'):
            if m in r.index and np.isfinite(r.loc[m, 'diff_mean']):
                x = r.loc[m]
                parts.append(f'{METHOD_LABEL[m]} {100 * x["mean"]:.1f} % ({METHOD_LABEL[m]} - PI-NN '
                             f'{100 * x["diff_mean"]:+.1f} ± {100 * x["diff_se"]:.1f} pp, PI-NN better on '
                             f'{int(x["n_pinn_better"])}/{int(x["n_folds"])}, Wilcoxon p {x["p_wilcoxon"]:.2g})')
        cov = checks.get('known_product_coverage')
        ctxt = ''
        if cov is not None and len(cov):
            c = cov[cov['family'] == 'known-product folds'].set_index('key')['coverage']
            n = cov[cov['family'] == 'new-product folds']['n']
            ctxt = (f' {int(n.iloc[0]) if len(n) else 0} of the {int(cov[cov["key"] == "mpsp"]["n"].sum())} '
                    'out-of-fold scores that calibrate the shipped intervals come from these new-product folds; '
                    'coverage of the shipped quantiles on the known-product folds\' samples (pooled; '
                    '`tables/interval_coverage_known_product.csv`): '
                    + ', '.join(f'{k} {c[k]:.2f}' for k in c.index) + '.')
        if pin is not None:
            L.append(f'**LOCO without the new-product folds.** {len(newp)} of the {int(pin["n_folds"]) + len(newp)} '
                     f'LOCO folds ({", ".join(f"`{c}`" for c in newp)}) hold out the only development '
                     'configuration of their fermentation product, so they are new-product folds, not new '
                     f'configurations of a known product. Over the other {int(pin["n_folds"])} folds '
                     f'(`tables/cv_loco_known_product.csv`): PI-NN {100 * pin["mean"]:.1f} % (± '
                     f'{100 * pin["se"]:.1f} s.e.), ' + ', '.join(parts) + '.' + ctxt + '\n')
    pf = checks.get('loco_per_family')
    if pf is not None and len(pf):
        L.append('LOCO fold MAPE of the PI-NN per family (`tables/cv_loco_per_family.csv`; mean, range over the '
                 'family\'s folds): ' + '; '.join(
                     f'{x.family} {100 * x.mean:.1f} %' + (f' ({100 * x.min:.1f}-{100 * x.max:.1f}, {x.n_folds} folds)'
                                                           if x.n_folds > 1 else ' (1 fold)')
                     + (' [new-product fold]' if x.n_folds == 1 and x.min_config in newp else '')
                     for x in pf.itertuples()) + '.\n')
    return '\n'.join(L)


def test_without_twins_section(checks):
    """Section 6: the OOD test without the held-out configurations that are
    input twins of development configurations (section 2)."""
    import numpy as np
    oodt, tp = checks.get('ood_twins') or [], checks.get('test_no_twin')
    if not oodt or tp is None or not len(tp):
        return ''
    t = tp.drop(columns=[c for c in ('set',) if c in tp]).copy()
    for c in ('mpsp_mape', 'diff_reference_minus_method', 'ci_lo', 'ci_hi'):
        t[c] = 100 * t[c]
    t = t[['method', 'n_configs', 'n_samples', 'mpsp_mape', 'diff_reference_minus_method', 'ci_lo', 'ci_hi',
           'frac_boot_reference_better']]
    L = [f'**OOD test without the input twin** {", ".join(f"`{c}`" for c in oodt)} (it encodes identically to a '
         'development configuration, section 2; report only, from the saved test outputs; '
         '`tables/final_test_without_twins.csv`), same columns as above. Its rows in the per-configuration tables '
         'of this section (error, coverage, applicability domain) are those of a configuration whose input is in '
         'training:\n', md_table(t) + '\n']
    cv_ = checks.get('test_no_twin_cov')
    if cv_ is not None and len(cv_):
        c = cv_[cv_['method'] == 'pinn']
        if len(c):
            L.append('PI-NN interval coverage without it (sample-simultaneous for MPSP; pooled / mean over '
                     'configurations; `tables/final_intervals_without_twins.csv`): '
                     + ', '.join(f'{x.key} {x.coverage:.2f} / {x.coverage_cfg_mean:.2f}' for x in c.itertuples())
                     + '.')
    ad = checks.get('ad_no_twin')
    if ad is not None and len(ad):
        x = ad.iloc[0]
        f = lambda v: f'{100 * v:.1f} %' if np.isfinite(v) else '–'
        L.append(f'Applicability domain without it (`models/final_dev`; '
                 '`tables/applicability_domain_test_without_twins.csv`): in-domain fraction '
                 f'{x["frac_in_domain"]:.2f} of {int(x["n_samples"])} samples, MPSP MAPE {f(x["mape_in_domain"])} in '
                 f'/ {f(x["mape_out_of_domain"])} out of the domain.')
    return '\n'.join(L) + '\n'


def _test_without_twin_clause(checks):
    """'; without the input twin ...' clause of the OOD-test summary bullet."""
    oodt, tp = checks.get('ood_twins') or [], checks.get('test_no_twin')
    if not oodt or tp is None or not len(tp) or not (tp['method'] == 'pinn').any():
        return ''
    tw = checks.get('twins')
    mates = sorted({x.twin if x.config_id in oodt else x.config_id for x in tw.itertuples()
                    if x.crosses_split and (x.config_id in oodt or x.twin in oodt)}) if tw is not None else []
    r = tp[tp['method'] == 'pinn'].iloc[0]
    txt = (f'; without {", ".join(oodt)}, which is input-identical to the development configuration '
           f'{", ".join(mates)} (section 2) and so not new in input space: PI-NN {100 * r["mpsp_mape"]:.1f} % '
           f'(n = {int(r["n_samples"])} samples)')
    return txt + _test_vs_baseline(tp, tp['set'].iloc[0])


def _no_twin_cov(checks, key):
    """', test_ood without the input twin x / y' (PI-NN interval coverage)."""
    c = checks.get('test_no_twin_cov')
    if c is None or not len(c):
        return ''
    c = c[(c['method'] == 'pinn') & (c['key'] == key)]
    return (f', test_ood without the input twin {c["coverage"].iloc[0]:.2f} / {c["coverage_cfg_mean"].iloc[0]:.2f}'
            if len(c) else '')


def _known_family_text(checks, ms):
    """Clause of the family-level LOCO comparison over the known-product
    families only (without the new-product folds)."""
    bk = checks.get('baseline_family_known')
    if bk is None or not len(bk):
        return ''
    r = bk.set_index('method')
    ms = [m for m in ms if m in r.index]
    if not ms:
        return ''
    n = int(r['n_families'].max())
    return (f'. The mean leans on the two new-product families ({r["excluded_folds"].iloc[0]}); over the {n} '
            'known-product families (the use "a new configuration of a known product"): '
            + ', '.join(f'{METHOD_LABEL[m]} - PI-NN {100 * r.loc[m, "mean"]:+.1f} pp (median '
                        f'{100 * r.loc[m, "median"]:+.1f}), PI-NN better in {int(r.loc[m, "n_pinn_better"])}/{n}'
                        for m in ms)
            + (' -- the PI-NN is better in every known-product family' if all(int(r.loc[m, 'n_pinn_better']) == n
                                                                         for m in ms)
               else ' -- the direction is not the same in every known-product family')
            + f', but {n} families cannot reach p < 0.05 (smallest attainable two-sided p '
            f'{r["p_min_attainable"].max():.3f})')


def loco_structure_summary(checks):
    """Summary bullets on the LOCO fold structure (family clustering,
    new-product folds, per-family errors)."""
    import numpy as np
    L = []
    bf = checks.get('baseline_family')
    if bf is not None and len(bf):
        r = bf.set_index('method')
        ms = [m for m in ('plain', 'knn') if m in r.index]
        L.append(f'* The LOCO folds cluster in {int(r["n_families"].max())} families, so the fold-level p-values '
                 'above assume independent configurations; at the family level (fold differences averaged per '
                 'family, then paired over the families): '
                 + ', '.join(f'{METHOD_LABEL[m]} - PI-NN {100 * r.loc[m, "mean"]:+.1f} pp, PI-NN better in '
                             f'{int(r.loc[m, "n_pinn_better"])}/{int(r.loc[m, "n_families"])} families (Wilcoxon p '
                             f'{r.loc[m, "p_wilcoxon"]:.2g})' for m in ms)
                 + (' -- the LOCO advantage over these baselines holds at the family level'
                    if all(r.loc[m, 'mean'] > 0 and np.isfinite(r.loc[m, 'p_wilcoxon']) and r.loc[m, 'p_wilcoxon'] < 0.05
                           for m in ms) else ' -- an advantage is not established at the family level for all of them')
                 + '; median over the families ' + ', '.join(f'{METHOD_LABEL[m]} {100 * r.loc[m, "median"]:+.1f} pp'
                                                             for m in ms)
                 + _known_family_text(checks, ms) + ' (all baselines: section 5).')
    kp, newp = checks.get('loco_known_product'), checks.get('loco_new_product_folds') or []
    if kp is not None and len(kp) and newp and (kp['method'] == 'pinn').any():
        r = kp.set_index('method')
        pf = checks.get('loco_per_family')
        fam = ('; PI-NN per family ' + ', '.join(f'{x.family} {100 * x.mean:.1f}' for x in pf.itertuples()
                                                if x.min_config not in newp) + ' %'
               if pf is not None and len(pf) else '')
        L.append(f'* LOCO without its {len(newp)} new-product folds ({", ".join(newp)}: the only configuration of '
                 f'their product): PI-NN {100 * r.loc["pinn", "mean"]:.1f} % over {int(r.loc["pinn", "n_folds"])} '
                 'folds vs. ' + ', '.join(f'{METHOD_LABEL[m]} {100 * r.loc[m, "mean"]:.1f} %'
                                          for m in ('plain', 'knn') if m in r.index) + fam + ' (section 5).')
    return '\n'.join(L)


def twin_section(checks):
    """Section 2 paragraph on configurations with identical model inputs
    (``clustered_and_twin_checks``)."""
    tw = checks.get('twins')
    if tw is None:
        return ''
    near = checks.get('twin_nearest') or {}
    nd = (f' No other pair of configurations is closer than {near["distance"]:.2f} (`{near["config_id"]}` / '
          f'`{near["other"]}`; Euclidean distance of the standardized inputs).' if near.get('config_id') else '')
    if not len(tw):
        return ('**Input twins** (checked after modelling with the features of `models/final_dev`; report only): '
                'no two configurations encode to the same model input.' + nd + '\n')
    pairs = '; '.join(f'`{x.config_id}` ({x.set}) = `{x.twin}` ({x.twin_set})' for x in tw.itertuples())
    n = int(tw['n_features'].iloc[0])
    txt = ('**Input twins** (checked after modelling with the features of `models/final_dev`, baseline rows; '
           f'report only; `tables/input_twins.csv`): {len(tw)} pair(s) of configurations encode to the SAME model '
           f'input (0 of {n} features differ): {pairs}. They differ only in settings that are not inputs (the '
           'saccharification mode SHF / SSCF of the lactic configurations), so the model predicts the same '
           'coefficients for both.' + nd)
    oodt = checks.get('ood_twins') or []
    if oodt:
        txt += (f' **The held-out {", ".join(f"`{c}`" for c in oodt)} is therefore not new in input space**: its '
                'twin is a development configuration, so its OOD test result measures the unmodelled difference '
                'between the twins (a label shift at identical inputs), not generalization to a new configuration. '
                'Section 6 reports the OOD test also without it; the split rule (fixed before modelling) is '
                'unchanged.')
    lt = checks.get('loco_twin_folds') or []
    if lt:
        txt += (f' In leave-one-configuration-out, the folds {", ".join(f"`{c}`" for c in lt)} have their twin in '
                'training (section 5).')
    return txt + '\n'


def ad_section(ad, per):
    """Markdown: applicability-domain validation on the test sets."""
    if ad is None or not len(ad):
        return ''
    L = ['**Applicability domain on the test sets** (`models/final_dev`, features as in training; '
         '`in_domain` = kNN distance <= training threshold and no input outside the training range; report '
         'only; `tables/applicability_domain_test.csv`). Fold models of the CV are not saved, so the domain '
         'is validated on the test sets only:\n']
    t = ad.copy()
    for c in ('mape_in_domain', 'mape_out_of_domain'):
        t[c] = 100 * t[c]
    L.append(md_table(t) + '\n\n(MAPE in %.)\n')
    if per is not None and len(per):
        q = per[per['set'] == 'test_ood'].copy()
        q['mpsp_mape'] = 100 * q['mpsp_mape']
        L.append('Per held-out test configuration:\n')
        L.append(md_table(q) + '\n')
    L.append('The flag is informative (higher errors out of the domain, positive rank correlation of `ad_ratio` '
             'with the error on the OOD test) but not a guarantee: in-domain samples of a new configuration can '
             'still have large errors, and every training-like query is in the domain.\n')
    return '\n'.join(L)


def test_reuse_notice(ctx, plan=None, regen=None):
    """Markdown notice of the report: the fourth use of the frozen test split
    (THIS_RUN after PRIOR_RUNS), what changed in between and why the test
    numbers are optimistic. ``regen``: ``report_checks.regenerated_rows``."""
    mk = ctx.get('mk_extra') or {}
    model = [lab for on, lab in (
        ((ctx.get('lk_extra') or {}).get('ch_space') == 'surface',
         'the channel / price-intensity loss term in surface space with a fixed weight'),
        (bool(mk.get('power_gate')), 'the power gate (no electricity export without a turbogenerator)'),
        (bool(ctx.get('target_mask')), 'masking the known defective channel targets in the loss ('
         + '; '.join(f"{k}: {', '.join(v)}" for k, v in (ctx.get('target_mask') or {}).items()) + ')'),
        (bool(mk.get('clip_exempt')), 'never clipping the exact final-product stoichiometric ratio'),
        (float((plan or {}).get('hpo_lopo_weight', 0.) or 0.) > 0,
         'leave-one-product-out folds in the HPO criterion')) if on]
    dev3 = settings_vs_v3(ctx, plan)
    prior = [r for r, _ in PRIOR_RUNS]
    prov = '; '.join(f'{r}: {p}, `results/{r}/`' for r, p in PRIOR_RUNS)
    return (f'> **{TEST_USE.capitalize()} use of the frozen test split.** This is the {TEST_USE} run of the '
            f'pipeline that evaluates on the frozen test split (`split.json`, unchanged): {THIS_RUN}, after '
            f'{_and(prior)} ({prov}; each folder keeps that run\'s report and headline tables). Between the runs '
            'the pipeline was changed on development-set evidence only: v1 -> v2 the training protocol and the '
            'model (DESIGN.md section 6, round 4); v2 -> v3 the training data, by the data corrections found by '
            'the final review of v2 (D15: the succinic neutral configurations never simulated the salt formation, '
            'acidulation and gypsum filtration, regenerated; H1: 47 TAL_SA samples lost to a failed simulator '
            'reset, re-run; DESIGN.md section 1)'
            + (', and the model (round 5, DESIGN.md section 6: ' + '; '.join(model) + ')' if model else '')
            + f'; v3 -> {THIS_RUN} '
            + ('ONLY the training data' if not dev3 else
               'the training data, and the settings (' + ', '.join(f'`{f}`' for f in dev3) + ')')
            + ', by the simulation corrections found by the final review of v3 (D16: the TAL_SA_sugarcane system '
            'fed the juice recycle to the cane mill twice; D17: BioSTEAM\'s anaerobic IC reactor of the '
            'wastewater treatment counted its biogas twice, in every HP and TAL configuration, ethanol_sugarcane, '
            'IBO_EtOH_corn and oxalic_neutral_sugarcane; the 27 affected configurations regenerated; DESIGN.md '
            'section 1)' + ('; the model, training protocol and hyperparameter search are those of v3, the HPO '
                            're-run on the corrected rows' if not dev3 else '')
            + f'; summary in "{CHANGES_TITLE}" below. No test sample was used to fit, tune, select or calibrate '
            f'anything. But the test results of {_and(prior)} had been seen: the v1 results set the round-4 '
            'agenda, the review that found the v2 -> v3 data defects and set the round-5 agenda read the v2 '
            'report, and the final review that found D16 and D17 read the v3 report, test results included. '
            f'The test split is therefore an adaptively re-used hold-out: the {THIS_RUN} test MAPE and the '
            f'{" -> ".join(ALL_RUNS)} test changes are likely optimistic and are not an independent '
            'confirmation of the changes. The data corrections also replace the test rows of the regenerated '
            'configurations: v2 -> v3 changed the test sets (`n_test` below), and v3 -> ' + THIS_RUN
            + f' corrected the simulated values of {_regen_test_clause(regen)} (the same sample ids), so no two '
            'runs since v2 are scored on identical test targets. The primary evidence is the development CV and '
            'its matched LOCO estimate of the OOD test (themselves optimistic after selection). '
            f'{THIS_RUN} was fixed before its test evaluation (code version, protocol and settings in '
            f'`plan.json`) and is the released model whatever the comparison with {" / ".join(prior)} shows: '
            'that comparison is not used to choose a version, and the v3 checkpoints carry the D16 / D17 bias of '
            'their training data. The v1 release intervals had been calibrated partly on the test '
            'configurations, which v2 and later no longer do.\n')


def publication_notes(pub, samples=None):
    """Known causes of baseline-coefficient differences from the publication
    (data only): succinic acid at neutral pH (D15, by design: the
    publication's process did not simulate the salt formation) and the HP
    titer- and yield-independent cost a′ (fermentation capital)."""
    import numpy as np
    L = []
    cid = pub['config_id'].astype(str)
    neu, low = cid.str.startswith('succinic_neutral'), cid.str.startswith('succinic_') & \
        ~cid.str.startswith('succinic_neutral')
    if neu.any() and low.any():
        L.append(f'Succinic acid at neutral pH differs by design (D15, DESIGN.md section 1): the regenerated '
                 'configurations simulate the salt formation, acidulation and gypsum filtration that the '
                 "publication's TRY script (like the module default) did not, so their a′ is "
                 f'{pub.loc[neu, "a_ratio"].min():.2f}-{pub.loc[neu, "a_ratio"].max():.2f} x the publication\'s, '
                 f'against {pub.loc[low, "a_ratio"].min():.2f}-{pub.loc[low, "a_ratio"].max():.2f} x for the '
                 'low-pH configurations.')
    hp = cid.str.startswith('HP_')
    if hp.any():
        s = f'HP: a′ is {pub.loc[hp, "a_ratio"].min():.2f}-{pub.loc[hp, "a_ratio"].max():.2f} x the publication\'s'
        if samples is not None and 'ch__capex_fermentation__a' in samples:
            b = samples[(samples['config_id'].astype(str).isin(set(cid[hp])))
                        & (samples['sample_kind'].astype(str) == 'baseline')]
            if len(b):
                sh = (b['ch__capex_fermentation__a'] / b['a']).to_numpy(float)
                s += (f'; the fermentation capital is a median {np.nanmedian(sh):.0%} of the HP baseline a′ '
                      '(`ch__capex_fermentation__a`) and scales with 1/r, so the HP productivity response is '
                      "larger than the publication's as well (DESIGN.md section 1)")
        L.append(s + '. The module or BioSTEAM change behind the larger HP fermentation capital has not been '
                 'identified (the Bioindustrial-Park patches change no process value, '
                 '`data_generation/bip_patches/README.md`); the productivity conclusions for HP carry this '
                 'uncertainty.')
    from fermentation_insights.coeff_predictor import report_checks as RC
    d17 = RC.corrected_configs(samples)
    if d17 and 'b_ratio' in pub:
        on, rest = cid.isin(d17), ~cid.isin(d17)
        rng = lambda m: (f'{pub.loc[m, "b_ratio"].min():.2f}-{pub.loc[m, "b_ratio"].max():.2f} x (median '
                         f'{_nanmedian(pub.loc[m, "b_ratio"]):.2f})')
        if on.any() and pub.loc[on, 'b_ratio'].notna().any():
            L.append(f'The {int(on.sum())} of these configurations regenerated with the D17 correction (harness v4; '
                     'DESIGN.md section 1) no longer burn the biogas that BioSTEAM\'s anaerobic IC reactor counted '
                     "twice, which the publication's modules very likely did as well"
                     + (', and TAL_SA_sugarcane no longer feeds its mill the juice recycle twice (D16)'
                        if 'TAL_SA_sugarcane' in set(cid[on]) else '')
                     + ': their b′ is ' + rng(on) + " the publication's"
                     + (', the other configurations\' ' + rng(rest) if rest.any() and pub.loc[rest, 'b_ratio'].notna().any()
                        else '') + ' (the correction raises b′ most on corn and corn stover).')
    return ' '.join(L)



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
         'feedstock-input clipping `clip_groups`)' if ctx.get('mk_extra', {}).get('fs_prior') else '')
      + ('; power gate (`ModelConfig.power_gate`: net electricity intensity >= 0 where electricity is not a '
         'coproduct)' if ctx.get('mk_extra', {}).get('power_gate') else '')
      + ("; pure-b' feedstock intensity (`ModelConfig.fs_b_only`)" if ctx.get('mk_extra', {}).get('fs_b_only')
         else '')
      + ('; never clipped: ' + ', '.join(f'`{c}`' for c in ctx['mk_extra']['clip_exempt'])
         + ' (`ModelConfig.clip_exempt`)' if ctx.get('mk_extra', {}).get('clip_exempt') else '')
      + (f"; channel loss `LossConfig.ch_space={ctx['lk_extra'].get('ch_space')!r}` with fixed "
         f"`w_ch={ctx['lk_extra'].get('w_ch')}` (chosen on out-of-fold channel fidelity, not an HPO dimension)"
         if ctx.get('lk_extra') else '')
      + ('; defective channel targets masked in the PI-NN loss (' + ', '.join(
          f"{k}: {', '.join(v)}" for k, v in ctx['target_mask'].items()) + ')' if ctx.get('target_mask') else '')
      + '.\n')
    if D['recent_raw_files']:
        P(f'> **Partial data:** {len(D["recent_raw_files"])} raw file(s) were still being written when the tables '
          'were built; the numbers describe the data available at that time.\n')
    protocol = ctx.get('protocol') or PROTOCOL
    P('**Where numbers come from.** (1) *CV estimates*: models re-trained inside each fold on parts of the '
      'development set only (selected hyperparameters, '
      f'{plan["cv_members"]}-member ensembles for the networks), with the same training procedure as the final '
      f'models (protocol `{protocol["name"]}`, section 3). (2) *Final-test estimate*: the selected PI-NN '
      f'ensemble ({plan["final_members"]} members, `models/final_dev`) and the baselines, trained once on the whole '
      'development set and evaluated on the frozen test split. No test sample was used to fit, tune, select or '
      'calibrate anything, but the test results of the earlier runs ('
      + ', '.join(r for r, _ in PRIOR_RUNS) + ') were seen before this run was designed, so the test split is a '
      're-used hold-out (notice below). '
      f'(3) *Release model*: {plan["release_members"]}-member PI-NN ensemble with the same hyperparameters and '
      'procedure trained on ALL data '
      '(`models/release`); it has no independent test estimate -- the CV and final-test numbers are the estimates '
      'of its accuracy (slightly pessimistic: it is trained on more data). MAPE = mean absolute percentage error of '
      'MPSP over the held-out samples\' simulated TY points; '
      '"fit MAPE" is the irreducible error of the per-sample 4-coefficient fit on the same points.\n')
    v1 = checks.get('v1') or {}
    prior = checks.get('prior') or {'v1': v1}
    from fermentation_insights.coeff_predictor import report_checks as RC
    regen = RC.regenerated_rows(D.get('samples'), dev['sample_id'],
                                [c['config_id'] for c in sp.get('test_configs', [])])
    dev3 = settings_vs_v3(ctx, plan)
    P(test_reuse_notice(ctx, plan, regen))
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
          f'{100 * base[bm]:.1f} %; irreducible fit MAPE {100 * _m(sch, "pinn", "fit_mape"):.1f} %'
          + baseline_paired_clause(checks.get('baseline_paired'), sch, bm) + '.')
    P(loco_structure_summary(checks))
    lm = checks.get('loco_matched')
    if lm is not None and len(lm) and (lm['method'] == 'pinn').any():
        r = lm[lm['method'] == 'pinn'].iloc[0]
        lmn, lt = checks.get('loco_matched_no_twin'), checks.get('loco_twin_folds') or []
        tw = ''
        if lt and lmn is not None and len(lmn) and (lmn['method'] == 'pinn').any():
            r2 = lmn[lmn['method'] == 'pinn'].iloc[0]
            tw = (f'; {len(lt)} of its folds ({", ".join(lt)}) have an input-identical twin in training (section 2); '
                  f'without them {100 * r2["family_mean"]:.1f} %'
                  + (f' (lactic {100 * r2["lactic"]:.1f} % vs. {100 * r["lactic"]:.1f} % with them)'
                     if 'lactic' in r2 and np.isfinite(r2['lactic']) else ''))
        P(f'* LOCO restricted to the {", ".join(OOD_FAMILIES)} configurations (matched CV estimate of the OOD '
          f'test, family-balanced): PI-NN {100 * r["family_mean"]:.1f} % over {int(r["n_configs"])} '
          'configurations' + tw + '.')
    P(ablation_summary(checks.get('ablation_stats')))
    if len(fh):
        tp = checks.get('test_paired')
        for st, lab in (('test_ood', 'held-out configurations'), ('test_id', 'in-distribution samples')):
            r = fh[(fh['set'] == st) & (fh['method'] == 'pinn')]
            if len(r):
                P(f'* Final test, {lab} (PI-NN trained on dev): MPSP MAPE {100 * r["mpsp_mape"].iloc[0]:.1f} % '
                  f'(median {100 * r["mpsp_mape_median"].iloc[0]:.1f} %; fit MAPE {100 * r["fit_mape"].iloc[0]:.1f} %; '
                  f'n = {int(r["n_test"].iloc[0])} samples)' + _test_vs_baseline(tp, st)
                  + (_test_without_twin_clause(checks) if st == 'test_ood' else '') + '.')
    iv = fin.get('intervals')
    if iv is not None and len(iv):
        r = iv[(iv['method'] == 'pinn') & (iv['key'] == 'mpsp') & iv['calibration'].isin(['loco', 'holdout'])]
        if len(r):
            g = lambda st, c: float(r[r['set'] == st][c].iloc[0])
            sets = [st for st in ('test_ood', 'test_id') if (r['set'] == st).any()]
            src = ('dev LOCO out-of-fold scores' if r['calibration'].iloc[0] == 'loco'
                   else 'its own calibration hold-out (v1 protocol)')
            P(f'* PI-NN MPSP intervals (nominal {1 - ALPHA:.0%}, calibrated on {src}): realized '
              'sample-simultaneous coverage pooled over samples / mean over configurations '
              + ', '.join(f'{st} {g(st, "coverage"):.2f} / {g(st, "coverage_cfg_mean"):.2f}' for st in sets)
              + _no_twin_cov(checks, 'mpsp')
              + f'; cross-fitted estimate pooled {g(sets[0], "crossfit_pooled"):.2f} / mean over configurations '
              f'{g(sets[0], "crossfit_mean"):.2f}' + _family_undercoverage(checks, 'mpsp') + '.')
        rg = iv[(iv['method'] == 'pinn') & (iv['key'] == 'g') & iv['calibration'].isin(['loco', 'holdout'])]
        if len(rg):
            gg = lambda st, c: float(rg[rg['set'] == st][c].iloc[0])
            sets = [st for st in ('test_ood', 'test_id') if (rg['set'] == st).any()]
            fc = checks.get('family_coverage')
            low = (fc[(fc['key'] == 'g') & (fc['coverage'] < 1 - ALPHA - 0.1)] if fc is not None and len(fc)
                   else None)
            P('* PI-NN recovery (g) intervals UNDER-cover the new configurations of some families'
              + (' (coverage of the shipped quantile on the LOCO out-of-fold scores: '
                 + ', '.join(f'{x.family} {x.coverage:.2f}' for x in low.itertuples())
                 + '; the other families about 1.0)' if low is not None and len(low) else '')
              + ': realized coverage pooled / mean over configurations '
              + ', '.join(f'{st} {gg(st, "coverage"):.2f} / {gg(st, "coverage_cfg_mean"):.2f}' for st in sets)
              + _no_twin_cov(checks, 'g')
              + f'; cross-fitted estimate pooled {gg(sets[0], "crossfit_pooled"):.2f} / mean '
              f'{gg(sets[0], "crossfit_mean"):.2f}, 10 % quantile over configurations '
              f'{gg(sets[0], "crossfit_q10"):.2f} (section 6).')
    shc = checks.get('shipped_coverage')
    if shc is not None and len(shc):
        cv_ = lambda sch, k: shc[(shc['scheme'] == sch) & (shc['key'] == k)]
        parts = [f'{lab} MPSP {float(cv_(sch, "mpsp")["pooled"].iloc[0]):.2f}, g {float(cv_(sch, "g")["pooled"].iloc[0]):.2f}'
                 for sch, lab in (('feedstock', 'new feedstock'), ('product', 'new product'))
                 if len(cv_(sch, 'mpsp')) and len(cv_(sch, 'g'))]
        if parts:
            P('* The intervals are calibrated for new configurations of KNOWN feedstocks and products only: '
              'coverage of the shipped (release) quantiles on the leave-one-feedstock-out / leave-one-product-out '
              'out-of-fold scores (pooled over samples; nominal 0.90): ' + '; '.join(parts)
              + ' (section 6; `predict()` flags such rows `intervals_calibrated = False`).')
        summ = in_distribution_coverage(shc, checks.get('id_family_coverage'), checks.get('id_config_coverage'))[0]
        if summ:
            P(summ)
    P('')
    P(runs_comparison(fh, head, prior, regen, dev3))
    fmm = pd.to_numeric(D['samples'].get('fit_mape'), errors='coerce') if 'fit_mape' in D['samples'] else None
    P(changelog(ctx, plan, dict(
        regen=regen, data_hash=ctx.get('data_hash'), n_samples=len(D['samples']), n_points=len(D['points']),
        n_dev=len(dev), n_dev_configs=int(dev['config_id'].nunique()),
        n_configs=int(D['samples']['config_id'].nunique()),
        fit_mape_median=float(np.nanmedian(fmm)) if fmm is not None and fmm.notna().any() else None,
        selected=(hpo or {}).get('selected'), prior=prior.get('v3') or {})))
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
        rr = checks.get('random_cv_random_samples') or {}
        rtxt = ', '.join(f'{METHOD_LABEL[m]} {100 * rr[m]:.1f} %' for m in BASE_METHODS if m in rr)
        P(f'Within-configuration linear baseline for random samples (`evaluate.noise_floor`: per development '
          f'configuration with >= 10 random samples, {nf["n_configs"]} configurations, a 5-fold ridge on the '
          f'observable inputs): MPSP MAPE {100 * nf["floor_mape_cfg_gmean"]:.1f} % (configuration geo-mean; mean '
          f'{100 * nf["floor_mape_cfg_mean"]:.1f} %), vs. {100 * nf["const_mape_cfg_gmean"]:.1f} % for the '
          f'configuration mean and {100 * nf["fit_mape_cfg_gmean"]:.1f} % for the 4-coefficient fit itself. It is '
          'NOT a noise floor: it is the model-form error of a linear model fitted on 20-60 samples per '
          'configuration. Adding the unobserved biorefinery-specific `bp_*` draws to its inputs makes it worse '
          f'({100 * nf.get("with_bp_mape_cfg_gmean", float("nan")):.1f} %), and the networks predict random '
          'samples better in random CV (random samples, configuration geo-mean, mean over folds: '
          f'{rtxt}). The CV tables report spec-determined (baseline / productivity: every `bp_*` at baseline) and '
          'random samples separately (`*_spec`, `*_random`). `tables/noise_floor.csv`.\n')
    yt = checks.get('yt_check')
    if yt is not None and len(yt):
        dep = yt[yt['departs']]
        ok = yt[~yt['departs']]
        dev_ = lambda x: max(abs(x[c] - 1) for c in ('y_ratio_min', 'y_ratio_max', 't_ratio_min', 't_ratio_max'))
        P('**Measured vs. nominal yield and titer** (ADAPTERS.md checklist item 2; per simulated point '
          'F_ferm_product / F_sugar_fed / y and F_ferm_product / V_broth / t, i.e. the fermentation product leaving '
          'fermentation per sugar fed and per broth volume over the nominal y and t; `tables/yt_definition_check.csv`, '
          'min / median / max per configuration): '
          + (f'{len(ok)} configurations within ±{100 * max(dev_(x) for _, x in ok.iterrows()):.1f} %. ' if len(ok)
             else '')
          + (f'{len(dep)} depart by more than 1 %: ' + '; '.join(
              f'`{x.config_id}` y {x.y_ratio_min:.3f}-{x.y_ratio_max:.3f}, t {x.t_ratio_min:.3f}-{x.t_ratio_max:.3f}'
              for x in dep.itertuples()) + '. Known causes (DESIGN.md section 1, "Measured vs. nominal yield and '
             'titer"): TAL, D8 (the spiked acetate counts as substrate in y but not in F_sugar_fed); lactic yield, '
             'residual seed-train sugars recycled with the inoculum; ethanol titer, the ethanol of the CO2 vent counts '
             'in F_ferm_product but not in the module\'s broth titer. The lactic SSCF titer ratio below 1 is not '
             'explained. The coefficients are fitted on the nominal y and t (the modules\' own definitions), which '
             'are the inputs of the unifying equation.' if len(dep) else '') + '\n')
    sb = checks.get('sugar_check')
    if sb is not None and len(sb):
        ex = sb[sb['exceeds']]
        ok = sb[~sb['exceeds']]
        P('**Sugar fed vs. sugar available** (per simulated point F_sugar_fed / (F_feedstock x '
          'fs_w_fermentable_sugar_eq x sugar availability), the availability being the physics prior '
          '`features.sugar_availability`; `report_checks.sugar_balance_check`, `tables/sugar_balance_check.csv`): '
          + (f'{len(ok)} configurations at {ok["ratio_min"].min():.2f}-{ok["ratio_max"].max():.2f} (below 1: '
             'sugar losses before fermentation). ' if len(ok) else '')
          + ('MORE sugar fed than the feedstock holds (> 1.01): ' + '; '.join(
              f'`{x.config_id}` {x.ratio_min:.4f}-{x.ratio_max:.4f}' for x in ex.itertuples())
             + ' -- a mass-balance defect of the process model, like D16 (DESIGN.md section 1: the '
             'TAL_SA_sugarcane system of the v3 data appended the juice recycle to the crushing mill a second '
             'time, 1.088; corrected and regenerated since).'
             if len(ex) else 'None exceeds 1.01 (the v3 data had TAL_SA_sugarcane at 1.088: data defect D16, its '
             'system fed the juice recycle to the crushing mill twice; corrected and regenerated since, DESIGN.md '
             'section 1).') + '\n')
    cov = checks.get('data_coverage')
    if cov is not None and len(cov):
        short = cov[(cov['retained'] < cov['planned'].fillna(0)) | (cov['empty_reset_failed'] > 0)]
        P(f'Sample coverage per configuration (`tables/data_coverage.csv`; planned = 1 baseline + productivity + '
          'random samples of `data_generation/scripts/run_all.sh`, completed = `sample_end` written, '
          '`empty_reset_failed` = completed without any simulated point because the system state could not be '
          'restored after a convergence failure -- a defect of the harness of the earlier data, fixed before '
          'the v3 run (such samples are retried and were re-run, H1); retained = after the quality filters). '
          + (f'{len(short)} of {len(cov)} configurations keep fewer samples than planned:\n' if len(short)
             else 'Every configuration keeps its planned samples.\n'))
        if len(short):
            P(md_table(short) + '\n')
            one = short[short['retained'] <= 1]['config_id'].tolist()
            P((f'{", ".join(f"`{c}`" for c in one)} ' + ('keeps' if len(one) == 1 else 'keep') + ' only one '
               'sample (a LOCO fold with a single held-out sample). ' if one else '')
              + 'Samples rejected by the quality filters are not missing at random where a feasibility limit '
              'interacts with the sampled parameters (e.g. the lactic SSCF configurations lose the random samples '
              'of low feedstock carbohydrate content, which reach the 40 g/L minimum titer at too few yields). '
              'Known data defects and their effects: DESIGN.md section 1 ("Known data defects").\n')
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
    tw_txt = twin_section(checks)
    if tw_txt:
        P(tw_txt)
    # protocol
    P('## 3. Training protocol\n')
    if protocol['name'] == 'v2':
        P('Every network fit of this report -- HPO, every CV scheme (PI-NN, plain MLP, ablations), '
          '`models/final_dev` and `models/release` -- uses the same procedure (`experiments.PROTOCOL`): no '
          'calibration hold-out; each ensemble member early-stops on its own random 15 % sample-level hold-out '
          '(the hold-outs differ between members, so the ensemble sees every sample; patience '
          f'{plan["patience"]}, at least {plan["min_epochs"]} and at most {plan["epochs"]} epochs); no refit. '
          'The baselines fit exactly the rows the members train on. So the CV numbers estimate the training '
          'procedure that made the final models, with the selected hyperparameters: the HPO (section 4) chose '
          'them on the same development set and is not nested in the CV, so the CV is slightly optimistic for '
          'the whole procedure (HPO + fit). Intervals: section 6.\n')
    else:
        P('**v1 protocol** (`--v1-protocol`, `experiments.PROTOCOL_V1`): CV calibration (15 %) and early-stopping '
          f'hold-outs grouped per scheme by {INNER_GROUP_V1}; `final_dev` with a 15 % configuration-grouped '
          'calibration hold-out; HPO and release without one, members early-stopping on held-out configurations.\n')
    P(f'Ensemble members: CV {plan["cv_members"]}, final_dev {plan["final_members"]}, release '
      f'{plan["release_members"]} (equal: the conformal scores divide by the member spread); HPO '
      f'{plan.get("hpo_members", 1)}.\n')
    pt = checks.get('protocol')
    if pt is not None and len(pt):
        P('Best epochs (the restored early-stopping state) and epochs run per stage, network members:\n')
        P(md_table(pt) + '\n')
        P('(`members_at_max_epochs`: members that ran into the epoch limit; `members_best_le_20`: members whose '
          'best epoch was <= 20.) `tables/protocol_epochs.csv`.\n')
    prior = checks.get('prior') or {'v1': checks.get('v1') or {}}
    prev = next((r for r, _ in reversed(PRIOR_RUNS) if (prior.get(r) or {}).get('timings_fits') is not None), None)
    tf1 = (prior.get(prev) or {}).get('timings_fits') if prev else None
    tf2 = checks.get('timings_fits')
    if tf1 is not None and tf2 is not None and len(tf2):
        cmp = tf2.merge(tf1, on=['stage', 'method'], how='left', suffixes=('', f'_{prev}'))
        cmp = cmp[cmp['method'].isin(list(NN_BASES) + list(ABLATIONS))]
        cols = [c for c in ('stage', 'method', 'fits', 'mean_seconds', f'fits_{prev}', f'mean_seconds_{prev}',
                            'mean_member_epochs', f'mean_member_epochs_{prev}') if c in cmp]
        P(f'Cost per network fit, {THIS_RUN} vs {prev} (`results/{prev}/tables/timings_fits.csv`; seconds per fit, '
          'member epochs summed over the members):\n')
        P(md_table(cmp[cols]) + '\n')
    # hpo
    t = hpo['table']
    P('## 4. Hyperparameter selection (development set only)\n')
    P(f'Random search over {len(t)} candidates (candidate 0 = package defaults) with {plan["hpo_k"]}-fold '
      f'leave-configurations-out CV (folds stratified by biorefinery; {plan.get("hpo_members", 1)} network(s) per '
      'fit, trained with the protocol of section 3), criterion = fold mean of the configuration-balanced '
      'held-out MPSP MAPE (each held-out configuration counts once)'
      + (f', weighted {1 - hpo.get("lopo_weight", 0.):g} : {hpo.get("lopo_weight", 0.):g} between these LOCO folds '
         '(`criterion_loco`) and the leave-one-product-out folds (`criterion_lopo`, mean over the held-out products; '
         'DESIGN.md section 6, round 5)' if hpo.get('lopo_weight', 0.) > 0 else '') + ' (`criterion`). '
      + (f'Search size from the measured cost of candidate 0 and the {a.budget_hours} h budget '
         f'(estimate {plan["estimate"]["total_hours"]:.2f} h for the whole run on {a.threads} workers). '
         if plan.get('estimate') else '') + '\n')
    show = t[['candidate', 'width', 'depth', 'dropout', 'lr', 'weight_decay', 'w_mpsp', 'w_ch', 'w_phys',
              'colloc_ratio', 'early_stop', 'clip_inputs', 'criterion', 'criterion_se']
             + [c for c in ('criterion_loco', 'criterion_lopo') if hpo.get('lopo_weight', 0.) > 0 and c in t]
             + ['mpsp_mape', 'mpsp_mape_cfg_gmean', 'fit_mape', 'epochs', 'seconds_per_fit', 'selected']]
    P(md_table(show.sort_values('criterion'), pct_cols=('criterion', 'criterion_se', 'criterion_loco',
                                                        'criterion_lopo', 'mpsp_mape', 'mpsp_mape_cfg_gmean',
                                                        'fit_mape'))
      + '\n\n(MAPE columns in %.) **Selected:** ' + ', '.join(f'{k}={v}' for k, v in hpo['selected'].items())
      + '. `tables/hpo_candidates.csv`, `figures/hpo_search.png`.\n')
    P('**Only the PI-NN is tuned.** The plain MLP reuses the width, depth, learning rate, weight decay and '
      'dropout selected for the PI-NN; gradient boosting and kNN use fixed settings (`evaluate.py`), ridge its '
      'internal generalized-CV choice of the penalty. The PI-NN also profits from the non-nested selection on '
      'the development folds (section 5). Comparisons with the baselines are therefore tilted toward the '
      'PI-NN by an unknown amount.\n')
    # cv
    head = cv['headline']
    P('## 5. Generalization CV estimates (development set)\n')
    P('Selected hyperparameters; the same training procedure as the final models (section 3'
      + ('' if protocol['name'] == 'v2' else f'; v1 inner hold-outs grouped by {INNER_GROUP_V1}') + '). '
      f'Random K = {plan["random_k"]}; leave-one-configuration-out over '
      + (f'a stratified subset of {plan["loco_folds"]} configurations; ' if plan.get('loco_folds') else
         'all development configurations (PI-NN, plain MLP, baselines'
         + (f'; the PI-NN ablations on a stratified subset of {len(cv["ablation_folds"])}: '
            + ', '.join(f'`{c}`' for c in cv['ablation_folds']) if cv.get('ablation_folds') else '') + '); ')
      + 'productivity: train on the central 20-80 % quantile range of log(r / baseline r), test both tails. '
      'MAPE = mean over folds (± s.e.); pooled = over all out-of-fold samples; R² pooled out-of-fold. '
      'Caveats: the HPO selected the hyperparameters on this development set and is not nested in the CV; and '
      'the leave-one-configuration-out, -feedstock-out, -product-out and productivity folds are deterministic '
      f'(random CV: seed {a.seed}), so they are the folds on which rounds 4 and 5 diagnosed and selected the v2 and '
      f'v3 changes ({THIS_RUN} changed only the training data) -- post-selection estimates, optimistic, most of all '
      'leave-one-feedstock-out (4 folds) and leave-one-product-out (7 folds), whose v1 -> v2 gains are not '
      f'established at the fold level ("{CHANGES_TITLE}").\n')
    cols = ['scheme', 'method', 'n_folds_common', 'mpsp_mape', 'mpsp_mape_se', 'mpsp_mape_median', 'mpsp_mape_pooled',
            'fit_mape'] + [f'R2_{k}' for k in TARGETS]
    P(md_table(head[cols], pct_cols=('mpsp_mape', 'mpsp_mape_se', 'mpsp_mape_median', 'mpsp_mape_pooled',
                                     'fit_mape')) + '\n\n(MAPE columns in %.)\n')
    ab = checks.get('ablation_stats')
    if ab is not None and len(ab):
        P('**Ablations, paired over folds** (variant - full PI-NN fold MPSP MAPE, in percentage points; '
          '`n_variant_better` = folds where the variant has the lower error; `dominant_fold` = the fold with the '
          'largest |difference| and `mean_without_dominant` the mean without it; Wilcoxon signed-rank p; '
          '`tables/cv_ablation_paired.csv`, per fold `tables/cv_ablation_per_fold.csv`). The ablations ran on '
          'a stratified subset of the LOCO folds; the plain MLP is shown on all folds and on that subset:\n')
        t = ab.copy()
        for c in ('mean', 'median', 'se', 'dominant_delta', 'mean_without_dominant'):
            t[c] = 100 * t[c]
        P(md_table(t[['scheme', 'variant', 'subset', 'n_folds', 'mean', 'se', 't', 'median', 'n_variant_better',
                      'n_variant_worse', 'p_wilcoxon', 'dominant_fold', 'dominant_delta',
                      'mean_without_dominant']]) + '\n')
        P(ablation_summary(ab).replace('* ', '') + '\n')
    bp = checks.get('baseline_paired')
    if bp is not None and len(bp):
        P('**PI-NN vs. the baselines, paired over folds** (baseline - PI-NN fold MPSP MAPE, percentage points; '
          '> 0: the PI-NN is better; `n_pinn_better` / `n_baseline_better` = folds; Wilcoxon signed-rank p '
          '(>= 5 folds) and two-sided sign-test p; `dominant_fold` = the fold with the largest |difference|; '
          '`tables/cv_baselines_paired.csv`). The fold means of the summary are not a comparison by '
          'themselves: for a new product or feedstock one or two folds carry them:\n')
        t = bp.copy()
        for c in ('mean', 'median', 'se', 'dominant_delta', 'mean_without_dominant'):
            t[c] = 100 * t[c]
        P(md_table(t[['scheme', 'method', 'n_folds', 'mean', 'se', 't', 'median', 'n_pinn_better',
                      'n_baseline_better', 'p_wilcoxon', 'p_sign', 'dominant_fold', 'dominant_delta',
                      'mean_without_dominant']]) + '\n')
    P(loco_structure_section(checks))
    P('Per-coefficient MAE (pooled out-of-fold, natural units) and constraint-violation rates (`data (fit)` = '
      'the rates of the data themselves, i.e. of the per-sample fitted coefficients on the same held-out samples '
      'and with the same definitions; a method that reproduced the data exactly would show them):\n')
    vcols = ['scheme', 'method'] + [f'MAE_{k}' for k in TARGETS] + ['viol_mono_y', 'viol_mono_t', 'viol_positive',
                                                                     'viol_mono_r', 'viol_g_le_gmax']
    vcols += [c for c in ('coverage_MPSP',) if c in head and head[c].notna().any()]
    vcols += [c for c in ('crossfit_cov_mpsp', 'crossfit_cov_g') if c in head]
    P(md_table(with_data_rows(head[vcols], checks.get('data_viol'), 'cv', 'scheme')) + '\n')
    P(data_violation_note(checks.get('data_viol_by_config')))
    lm = checks.get('loco_matched')
    if lm is not None and len(lm):
        P(f'**LOCO restricted to the {", ".join(OOD_FAMILIES)} configurations** -- the families of the held-out '
          'test configurations, so this is the matched CV estimate of the OOD test (fold MAPE; `family_mean` = '
          'mean over the families of the mean over their configurations, as the OOD test holds out one '
          'configuration of each; `tables/cv_loco_ood_families.csv`):\n')
        P(md_table(lm, pct_cols=tuple(c for c in lm.columns if c not in ('method', 'n_configs'))) + '\n\n(%.)\n')
        lt, lmn = checks.get('loco_twin_folds') or [], checks.get('loco_matched_no_twin')
        if lt and lmn is not None and len(lmn):
            per = checks.get('loco_per_fold_pinn') or {}
            P(f'The folds {", ".join(f"`{c}`" for c in lt)} hold out a configuration whose input twin (section 2) '
              'is in training' + (' (PI-NN ' + ', '.join(f'{100 * per[c]:.1f} %' for c in lt if c in per) + ')'
                                   if per else '')
              + '; the held-out lactic test configuration has its twin in training as well, so these folds are '
              'the matched estimate of that one. Without them (`tables/cv_loco_ood_families_without_twin_folds.csv`'
              '; the lactic column is then the one lactic fold without a twin in training):\n')
            P(md_table(lmn, pct_cols=tuple(c for c in lmn.columns if c not in ('method', 'n_configs')))
              + '\n\n(%.)\n')
    cfx = cv.get('crossfit')
    if cfx is not None and len(cfx):
        sub = cfx[(cfx['weighting'] == 'sample') & cfx['method'].isin(NN_BASES)
                  & cfx['scheme'].isin(['config', 'random'])]
        if len(sub):
            P(f'**Cross-fitted interval coverage** (nominal {1 - ALPHA:.2f}): for every CV fold, the conformal '
              'quantile of the out-of-fold scores of all OTHER folds and the fraction of the fold\'s own samples '
              'covered (`train.crossfit_coverage`; a sample counts as covered for MPSP only if all its TY points '
              'are). `mean` = fold-balanced, `min` / `q10` = worst / 10 % quantile over folds, `frac_ok` = '
              'fraction of folds with coverage >= nominal. The final intervals use the LOCO (`config`) PI-NN '
              'scores; `tables/cv_coverage_crossfit.csv` (also with configuration weighting), out-of-fold scores '
              '`tables/cv_oof_scores.csv.gz`.\n')
            P(md_table(sub[['scheme', 'method', 'key', 'mean', 'pooled', 'min', 'q10', 'frac_ok', 'n_folds']])
              + '\n')
    lp = checks.get('lopo_per_product')
    if lp is not None and len(lp):
        P('Leave-one-product-out, per held-out product (MAPE %):\n')
        P(md_table(lp, pct_cols=tuple(c for c in lp.columns if c not in ('product',))) + '\n')
    fs_ = checks.get('fidelity')
    if fs_ is not None and len(fs_):
        P('**Out-of-fold channel and price-intensity fidelity** (PI-NN ensembles of the CV; '
          '`evaluate.channel_fidelity`, configuration-balanced; DESIGN.md section 6, round 5): `ch_<c>` = RMS '
          'error of the channel numerator over the sample\'s points / RMS target total numerator (price channels '
          'as price x intensity; `ch_rms` = root sum of squares over the channels); per price key the error of '
          'dMPSP/dp = intensity surface / g against the simulated per-point intensity: `sens` relative RMS, '
          '`bias` |median ratio - 1|, `resp` RMS (dMPSP/dp predicted - simulated) x p / MPSP; `*_valid` without '
          'the configurations whose feedstock targets are masked in the loss (DEFECTIVE_TARGETS). In %:\n')
        cols = ['scheme'] + [c for c in FIDELITY_REPORT_COLS + ('ch_feedstock_valid', 'sens_fs_valid')
                             if c in fs_]
        P(md_table(fs_[cols], pct_cols=tuple(cols[1:])) + '\n')
        pc = checks.get('fidelity_per_config')
        if pc is not None and len(pc):
            q = pc[pc['scheme'] == 'config']
            cc = ['config_id', 'n_samples'] + [c for c in ('bias_fs', 'resp_fs', 'bias_el', 'resp_el', 'share_el',
                                                            'bias_ng', 'resp_ng', 'ch_electricity')
                                               if c in q]
            P('Per held-out configuration (leave-one-configuration-out; medians over its samples; `bias` signed, '
              '`share_el` = simulated p_el dMPSP/dp_el / MPSP; the relative `bias_*` / `sens_*` of a price are '
              f'shown only where its simulated share of MPSP is at least {100 * FIDELITY_MIN_SHARE:.0f} % '
              '(`–` otherwise: the ratio is undefined; `resp_*` is the error in % of MPSP), in %; all schemes in '
              '`tables/cv_fidelity_per_config.csv`, per sample `tables/cv_channel_fidelity.csv.gz`:\n')
            P(md_table(q[cc].reset_index(drop=True), pct_cols=tuple(cc[2:])) + '\n')
    P('Full tidy results: `tables/cv_results.csv`; summary: `tables/cv_headline.csv`; out-of-fold predictions: '
      '`tables/cv_oof_predictions.csv.gz`.\n')
    # final
    fh = fin['headline']
    P('## 6. Final test estimate (models trained on the development set, evaluated once)\n')
    if len(fh):
        cols = ['set', 'method', 'n_train', 'n_test', 'mpsp_mape', 'mpsp_mape_median', 'fit_mape'] + \
               [f'R2_{k}' for k in TARGETS] + ['coverage_MPSP', 'viol_mono_y', 'viol_mono_t', 'viol_mono_r']
        P(md_table(with_data_rows(fh[fh['set'] != 'train_all'][cols], checks.get('data_viol'), 'test', 'set'),
                   pct_cols=('mpsp_mape', 'mpsp_mape_median', 'fit_mape'))
          + '\n\n(test_ood = held-out configurations, test_id = in-distribution samples; MAPE in %; '
          '`coverage_MPSP` = pointwise coverage of the MPSP intervals; `data (fit)` = violation rates of the data '
          'themselves, section 5.)\n')
        pc = checks.get('test_per_config')
        if pc is not None and len(pc):
            P('PI-NN per held-out test configuration:\n')
            P(md_table(pc, pct_cols=('mpsp_mape', 'mpsp_mape_median', 'fit_mape')) + '\n')
        tp = checks.get('test_paired')
        if tp is not None and len(tp):
            t = tp.drop(columns=['n_samples']).copy()
            for c in ('mpsp_mape', 'diff_reference_minus_method', 'ci_lo', 'ci_hi'):
                t[c] = 100 * t[c]
            P('**Paired comparison on the test sets** (MPSP MAPE pooled over the TY points, %; '
              '`diff_reference_minus_method` = PI-NN - method in percentage points with a 95 % percentile bootstrap '
              'CI that resamples the test configurations, the same draws for every method; '
              '`frac_boot_reference_better` = fraction of draws in which the PI-NN has the lower error; '
              '`tables/final_test_paired.csv`). With 4 held-out configurations the OOD test cannot rank the '
              'methods: its CIs are those of 4 clusters.\n')
            P(md_table(t) + '\n')
        P(test_without_twins_section(checks))
        P(intervals_section(fin, protocol))
        P(interval_extrapolation_section(checks.get('shipped_coverage'), cv.get('crossfit'),
                                         checks.get('family_coverage'), checks.get('config_coverage')))
        P(in_distribution_coverage(checks.get('shipped_coverage'), checks.get('id_family_coverage'),
                                   checks.get('id_config_coverage'))[1])
        P(ad_section(checks.get('ad'), checks.get('ad_per_config')))
    else:
        P('Final models failed (see log).\n')
    # release
    P('## 7. Release model\n')
    rel = fin['release']
    if rel['error']:
        P(f'Release training FAILED: {rel["error"]}\n')
    else:
        conf = fin.get('release_conformal') or {}
        rc = fin.get('release_check') or {}
        ins = fh[(fh['set'] == 'train_all')]
        ins_txt = f'{100 * ins["mpsp_mape"].iloc[0]:.1f} %' if len(ins) else 'n/a'
        P(f'`{os.path.relpath(rel["saved_model"], PKG_DIR)}`: {plan["release_members"]} members trained on all '
          f'{len(D["samples"])} samples (development + test) with the protocol of section 3; epochs run '
          f'{rel["epochs"]}, best epochs {rel.get("best_epochs")}. In-sample MPSP MAPE (NOT a generalization '
          f'estimate): {ins_txt}.\n')
        if conf:
            cfm = (conf.get('crossfit') or {}).get('mpsp') or {}
            P(f'Intervals: split-conformal quantiles from {conf.get("source")} ({conf.get("n")} samples of '
              f'{conf.get("n_calib_groups")} configurations, weighting `{conf.get("weighting")}`, cross-fitted MPSP '
              f'coverage {_fmt(cfm.get("pooled", float("nan")))} pooled over samples (the quantity the intervals '
              f'target), {_fmt(conf.get("mpsp_coverage_cv", float("nan")))} averaged over configurations '
              '(`conformal["mpsp_coverage_cv"]` of the checkpoint)): the same scores as `final_dev`, '
              'none of them from a test sample. The release ensemble trains on more data than the CV fold '
              'ensembles; both its errors and its member spread shrink with more data, so the direction of the '
              'resulting coverage error is unknown; flagged `coverage_guaranteed = False` in its checkpoint.\n')
        else:
            P('Intervals: none (no usable out-of-fold scores; see log).\n')
        P(f'Loads with `CoefficientPredictor.load`: {rc.get("ok")} '
          + (f'(e.g. {rc.get("example_config")} baseline: ' + ', '.join(f'{k}={v:.4g}' for k, v in
                                                                         rc.get('example', {}).items())
             + f'; intervals reported: {rc.get("intervals")}).' if rc.get('ok') else f'({rc.get("error")}).') + '\n')
        if rc.get('ok') and 'rows_max_rel' in rc:
            P(f'User API vs. the validated path (`experiments.api_consistency`): `CoefficientPredictor.predict` on all '
              f'{rc["n_rows"]} training rows reproduces `Ensemble.predict_samples` (behind every CV and test number '
              f'of this report) to a max relative difference of {rc["rows_max_rel"]:.1g}; the by-name query '
              '(`config_id`, names, productivity and prices only) of each of the '
              f'{rc["n_configs"]} configurations reproduces the MPSP of its baseline sample at its baseline yield '
              f'and titer to {rc["by_name_max_rel_mpsp"]:.1g} (e.g. {rc.get("example_config")} by name: '
              + ', '.join(f'{k}={v:.4g}' for k, v in rc.get('example_by_name', {}).items()) + ')'
              + (f', with the applicability domain of its baseline row (max |ad_ratio difference| '
                 f'{rc["by_name_ad_max_abs"]:.1g}, {rc["by_name_in_domain_mismatch"]} in_domain mismatches)'
                 if 'by_name_ad_max_abs' in rc else '')
              + (f'; without prices it uses the configuration\'s baseline electricity and natural-gas prices '
                 f'(MPSP difference {rc["by_name_no_prices_max_rel_mpsp"]:.1g})'
                 if 'by_name_no_prices_max_rel_mpsp' in rc else '')
              + (f'. New configurations: `models/final_dev` on the {rc["new_configs"]["n_rows"]} stored rows of '
                 f'the {rc["new_configs"]["n_configs"]} held-out test configurations (passed as stored, '
                 'not in its training data) reproduces the validated path to '
                 f'{rc["new_configs"]["max_rel"]:.1g} ({rc["new_configs"].get("in_domain_mismatch", 0)} '
                 'in_domain mismatches)' if isinstance(rc.get('new_configs'), dict) else '') + '.\n')
    tb = checks.get('twin_bias')
    if tb is not None and len(tb):
        t = tb.copy()
        for c in ('err_baseline_sample', 'err_mean'):
            t[c] = 100 * t[c]
        P('**Input twins in the release model** (section 2: identical inputs, so identical coefficients, e.g. '
          'a′ below): in-sample signed MPSP error, predicted / fitted - 1 at the sample\'s baseline yield and titer '
          '(fitted = the sample\'s 4-coefficient fit, the training target), for the baseline sample and as the '
          'mean over the configuration\'s samples, in % (`tables/release_twin_bias.csv`). The model averages the '
          'twins, so it over-predicts one and under-predicts the other:\n')
        P(md_table(t) + '\n')
    # checks
    P('## 8. Physical-consistency and interpretability checks\n')
    pub = checks['pub']
    P('**(i) Baseline coefficients vs. publication** (data only, no model; `sample_kind == baseline` of each '
      'configuration vs. `results/data/<config>_coefficients.npy`; differences are expected from the newer '
      f'BioSTEAM / harness). {len(pub)} of 32 publication configurations present.\n')
    if len(pub):
        med = {k: _nanmedian(pub[f'{k}_ratio']) for k in TARGETS}
        und = {k: pub.loc[pub[f'{k}_ratio'].isna(), 'config_id'].tolist() for k in TARGETS}
        P('Median ratio ours/publication (over the defined ratios): '
          + ', '.join(f"{k}′ {v:.3f}" if k != 'g' else f'g {v:.3f}' for k, v in med.items()) + '. '
          + (f'Ratios shown as `–` are undefined: one of the two values is near zero (below {PUB_RATIO_MIN:g} x '
             'the median |publication value| of the coefficient) or the signs differ -- compare the values '
             'themselves (' + '; '.join(f"{k}′: {', '.join(v)}" for k, v in und.items() if v) + ').'
             if any(und.values()) else '') + '\n')
        note = publication_notes(pub, D.get('samples'))
        if note:
            P(note + '\n')
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
    pre = checks.get('price_response')
    if pre:
        P('**(ii-b) Price response per configuration**: error of the predicted MPSP response to a 100 % change '
          'of each price, (dMPSP/dp predicted - simulated) x p / MPSP, in % of MPSP, at the configuration\'s '
          'baseline yield and titer (median over its samples; simulated = decomposed simulation intensities with '
          'the fitted g; `report_checks.price_response_error`, `tables/price_response_<model>.csv`). The '
          'predictions are exactly linear in the prices, but the intensities they multiply are learned, and are '
          'biased for some configurations even in-sample:\n')
        for name, t in pre.items():
            q = t.set_index('config_id')
            big = q[(q.abs() > 0.10).any(axis=1)]
            P(f'*{name}*: median |error| over configurations ' + ', '.join(
                f'{c.replace("_response_err", "").replace("_sensitivity_rel_err", " dMPSP/dp (relative)")} '
                f'{100 * q[c].abs().median():.1f} %' for c in q.columns)
              + (f'; configurations with an error above 10 % ({len(big)}):\n' if len(big) else '.\n'))
            if len(big):
                b = (100 * big).round(1).reset_index()
                P(md_table(b) + '\n')
    P(capacity_response_text(checks.get('capacity_response'), D['samples']))
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
    P('## 9. Figures\n')
    P('\n'.join(f'* `{os.path.relpath(f, a.out)}`' + (f' -- {FIGURE_CAPTIONS[os.path.basename(f)]}'
                                                        if os.path.basename(f) in FIGURE_CAPTIONS else '')
                 for f in figs) + '\n')
    P('## 10. Reproduction and timings\n')
    P(f'`{" ".join(ctx["argv"])}` (PYTHONPATH=/home/user/fermentation_insights). Python {platform.python_version()}, '
      f'{ctx["versions"]}. Every fit is cached under `cache/` (keyed by data hash, split, job and code version).\n')
    P(md_table(timings) + '\n')
    tc = checks.get('timings_compute_run')
    if tc is not None and len(tc):
        tot = tc[tc['stage'].astype(str).str.startswith('total')]['minutes']
        P('Every fit of this run was cached (the report was regenerated from the cached fits); the stage '
          'timings of the run that computed them are in `tables/timings_compute_run.csv`'
          + (f' (total {float(tot.iloc[0]):.1f} min).' if len(tot) else '.') + '\n')
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
    ap.add_argument('--loco-folds', type=int, default=0,
                    help='quick mode only: number of LOCO folds of every method (0 = all)')
    ap.add_argument('--members', type=int, default=5,
                    help='full mode: ensemble members of the CV, final_dev and release ensembles (equal, so that '
                         'the CV out-of-fold scores calibrate the final intervals)')
    ap.add_argument('--ablation-loco-folds', type=int, default=None,
                    help='LOCO folds (stratified subset of configurations) of the PI-NN ablations (default: '
                         f'{ABLATION_LOCO_FOLDS} in full mode, 4 in quick mode; 0 = all; the PI-NN and the plain MLP '
                         'always run on every LOCO fold; full mode: the budget may shrink it)')
    ap.add_argument('--hpo-members', type=int, default=2, help='full mode: ensemble members per HPO fit')
    ap.add_argument('--hpo-lopo-weight', type=float, default=None,
                    help='weight of the leave-one-product-out folds in the HPO criterion (default '
                         f'HPO_LOPO_WEIGHT = {HPO_LOPO_WEIGHT:g}: the v2 criterion, LOCO folds only; quick '
                         'mode: 0; DESIGN.md section 6, round 5)')
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
    # round 4 (v2; DESIGN.md section 6): each flag restores the v1 behaviour of one component
    ap.add_argument('--v1-g-inputs', action='store_true',
                    help="masking preset 'physics' (v1: the recovery g also sees the reagent prices) instead of "
                         "'physics_v2'")
    ap.add_argument('--no-r-structure', action='store_true',
                    help="PI-NN without the hard productivity structure (ModelConfig.r_structure='physics': "
                         "fermentation capital and electricity intensity = theta0(x) + theta1(x) r_ref/r, "
                         "theta1 >= 0, so dMPSP/dr <= 0 by construction)")
    ap.add_argument('--no-extrapolation-safeguards', action='store_true',
                    help='PI-NN without the sugar-availability prior of the feedstock intensity '
                         '(ModelConfig.fs_prior), the neutralization gate (ModelConfig.neutralization_gate) '
                         "and the clipping of the feedstock inputs (ModelConfig.clip_groups)")
    # round 5 (DESIGN.md section 6): each flag restores the v2 behaviour of one component
    ap.add_argument('--no-power-gate', action='store_true',
                    help='PI-NN without the power gate (ModelConfig.power_gate: net electricity intensity >= 0 '
                         'where electricity is not among the coproducts)')
    ap.add_argument('--no-clip-exempt', action='store_true',
                    help='clip the final-product stoichiometric ratio like any other input when clip_inputs >= 0 '
                         '(v2; instead of CLIP_EXEMPT)')
    ap.add_argument('--fs-b-only', action='store_true',
                    help="opt-in: pure-b' feedstock intensity (ModelConfig.fs_b_only; not adopted, DESIGN.md "
                         "section 6, round 5)")
    ap.add_argument('--v2-channel-loss', action='store_true',
                    help="v2 channel loss: channels and intensities in coefficient space with the HPO's w_ch, "
                         "instead of CHANNEL_LOSS (LossConfig.ch_space='surface', fixed w_ch)")
    ap.add_argument('--no-target-mask', action='store_true',
                    help='keep the known defective channel targets (DEFECTIVE_TARGETS: the TAL feedstock channel, '
                         'D8) in the PI-NN loss')
    ap.add_argument('--v2-model', action='store_true',
                    help='the v2 model settings: --no-power-gate --v2-channel-loss --no-target-mask '
                         '--no-clip-exempt (and no --fs-b-only)')
    ap.add_argument('--v1-protocol', action='store_true',
                    help='v1 training protocol (PROTOCOL_V1): configuration-grouped calibration and '
                         'early-stopping hold-outs in CV (per scheme INNER_GROUP_V1) and final_dev, instead '
                         'of sample-level early stopping without a calibration hold-out everywhere')
    ap.add_argument('--split-seed', type=int, default=None, help='only when the split is created (default 2026)')
    ap.add_argument('--replan', action='store_true', help='recompute plan.json from the pilot timings')
    ap.add_argument('--fresh', action='store_true', help='delete the fit cache first')
    ap.add_argument('--plan-only', action='store_true',
                    help='run data, split and the HPO pilot (candidate 0), write plan.json with the budget '
                         'estimate and stop')
    ap.add_argument('--force-pool', action='store_true', help=argparse.SUPPRESS)  # debug: pool even with 1 worker
    ap.add_argument('--nice', type=int, default=0, help='increase the process niceness')
    a = ap.parse_args(argv)
    if a.v2_model:
        a.no_power_gate = a.v2_channel_loss = a.no_target_mask = a.no_clip_exempt = True
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
        if POWER_GATE and not a.no_power_gate:
            mk_extra['power_gate'] = True
    if (FS_B_ONLY or a.fs_b_only) and not a.v2_model:
        mk_extra['fs_b_only'] = True
    if CLIP_EXEMPT and not a.no_clip_exempt:
        mk_extra['clip_exempt'] = list(CLIP_EXEMPT)
    if not a.v2_channel_loss and CHANNEL_LOSS:
        ctx['lk_extra'] = dict(CHANNEL_LOSS)
    if not a.no_target_mask and DEFECTIVE_TARGETS:
        ctx['target_mask'] = {k: list(v) for k, v in DEFECTIVE_TARGETS.items()}
    if not a.no_masking:
        mk_extra['channel_inputs'] = 'physics' if a.v1_g_inputs else MASKING_PRESET
    if not a.no_r_structure:
        mk_extra['r_structure'] = R_STRUCTURE
    if mk_extra:
        ctx['mk_extra'] = mk_extra
    ctx['protocol'] = dict(PROTOCOL_V1 if a.v1_protocol else PROTOCOL)
    T = lambda name, t0: stage_t.__setitem__(name, time.time() - t0)

    t0 = time.time()
    D = stage_data(a)
    T('data', t0)
    ctx['data_hash'] = D['hash']
    t0 = time.time()
    sp, dev = stage_split(a, D)
    ctx['split_hash'] = sp['hash']
    T('split', t0)
    plan = check_plan(base_plan(a))
    if a.hpo_candidates:
        plan['hpo_candidates'] = a.hpo_candidates
    if a.loco_folds:
        if not a.quick:
            raise SystemExit('--loco-folds subsets the LOCO folds of every method; the final intervals need the '
                             'out-of-fold scores of every development configuration (use --ablation-loco-folds '
                             'to subset only the ablations)')
        plan['loco_folds'] = a.loco_folds
    runner = Runner(a.threads, a.processed, D['hash'], os.path.join(a.out, 'cache'), a.force_pool)
    try:
        t0 = time.time()
        hpo = stage_hpo(a, plan, runner, dev, ctx, sp['n_test_ood'] + sp['n_test_id'], len(D['samples']))
        T('hpo', t0)
        _dump_json(dict(plan, data_hash=ctx['data_hash'], split_hash=ctx['split_hash'], threads=a.threads,
                        code_version=CODE_VERSION, protocol=ctx['protocol']['name'],
                        selected=hpo['selected'] if hpo else None), os.path.join(a.out, 'plan.json'))
        if hpo is None:
            log(f'plan-only: plan written to {os.path.join(a.out, "plan.json")}; stopping')
            return 0
        t0 = time.time()
        cv = stage_cv(a, plan, runner, dev, ctx, hpo['selected'], D)
        T('cv', t0)
        t0 = time.time()
        # final after CV: the final intervals are calibrated from the CV out-of-fold scores
        fin = stage_final(a, plan, runner, dev, ctx, hpo['selected'], sp, D, a.models, cv)
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
    cf = cv['crossfit']
    if len(cf):
        # the CV ensembles have no intervals of their own (no calibration hold-out): cross-fitted coverage
        m = cf[(cf['weighting'] == 'sample') & cf['key'].isin(['mpsp', 'g'])].pivot_table(
            index=['scheme', 'method'], columns='key', values='mean').add_prefix('crossfit_cov_').reset_index()
        cv['headline'] = cv['headline'].merge(m, on=['scheme', 'method'], how='left')
    cv['headline'] = cv['headline'].drop(columns=[c for c in cv['headline'] if c.startswith('coverage_')
                                                  and cv['headline'][c].isna().all()])
    cv['headline'].to_csv(os.path.join(tdir, 'cv_headline.csv'), index=False)
    if len(cv['oof']):
        cv['oof'].to_csv(os.path.join(tdir, 'cv_oof_predictions.csv.gz'), index=False)
    if len(cv['oof_scores']):
        cv['oof_scores'].to_csv(os.path.join(tdir, 'cv_oof_scores.csv.gz'), index=False)
    cv['crossfit'].to_csv(os.path.join(tdir, 'cv_coverage_crossfit.csv'), index=False)
    checks['fidelity'], checks['fidelity_per_config'] = fidelity_tables(cv.get('fidelity'))
    if len(checks['fidelity']):
        cv['fidelity'].to_csv(os.path.join(tdir, 'cv_channel_fidelity.csv.gz'), index=False)
        checks['fidelity'].to_csv(os.path.join(tdir, 'cv_fidelity_summary.csv'), index=False)
        checks['fidelity_per_config'].to_csv(os.path.join(tdir, 'cv_fidelity_per_config.csv'), index=False)
    checks['loco_matched'] = loco_matched(cv['results'], S)
    checks['loco_matched'].to_csv(os.path.join(tdir, 'cv_loco_ood_families.csv'), index=False)
    checks['protocol'] = protocol_table(cv, fin, hpo, plan)
    checks['protocol'].to_csv(os.path.join(tdir, 'protocol_epochs.csv'), index=False)
    fin['table'].to_csv(os.path.join(tdir, 'final_test_results.csv'), index=False)
    fin['headline'] = final_headline(fin['table'])
    fin['headline'].to_csv(os.path.join(tdir, 'final_test_headline.csv'), index=False)
    fin['intervals'], fin['intervals_per_config'] = final_intervals(fin, S)
    fin['intervals'].to_csv(os.path.join(tdir, 'final_intervals.csv'), index=False)
    fin['intervals_per_config'].to_csv(os.path.join(tdir, 'final_intervals_per_config.csv'), index=False)
    checks['prior'] = {run: load_run_tables(a.out, run, fallback=not a.quick) for run, _ in PRIOR_RUNS}
    checks['v1'] = checks['prior']['v1']
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
    checks['data_viol'] = data_violation_table(dev, plan, a.seed, fin, S, Pts)
    checks['data_viol'].to_csv(os.path.join(tdir, 'data_constraint_violations.csv'), index=False)
    from fermentation_insights.coeff_predictor import report_checks as _RC
    checks['data_viol_by_config'] = _RC.data_violations_by_config(S, Pts, list(D['samples']['sample_id'].astype(str)))
    checks['data_viol_by_config'].to_csv(os.path.join(tdir, 'data_constraint_violations_by_config.csv'), index=False)
    report_only_checks(a, D, cv, fin, checks, S, Pts, tdir, fp)
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
        # ablations vs the PI-NN / plain MLP on the same folds (the ablations may skip LOCO folds)
        if checks.get('ablation_per_fold') is not None and len(checks['ablation_per_fold']):
            _fig('ablations.png', fig_ablation, checks['ablation_stats'], checks['ablation_per_fold'])
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
    checks['timings_fits'] = tim_fit
    timings['note'] = ''
    timings.loc[len(timings)] = dict(stage='median s per (member-epoch x training sample)', minutes=np.nan,
                                     note=f'{sec_per:.3g}')
    # a run whose fits were all cached (e.g. a report regenerated after a text
    # change) keeps the stage timings of the run that computed them
    tpath, cpath = os.path.join(tdir, 'timings.csv'), os.path.join(tdir, 'timings_compute_run.csv')
    checks['timings_compute_run'] = None
    if runner is not None and getattr(runner, 'n_run', 1) == 0:
        if not os.path.exists(cpath) and os.path.exists(tpath):
            old = pd.read_csv(tpath)
            tot = old[old['stage'].astype(str).str.startswith('total')]['minutes']
            if len(tot) and float(tot.iloc[0]) > float(timings.iloc[-2]['minutes']):
                old.to_csv(cpath, index=False)
        if os.path.exists(cpath):
            checks['timings_compute_run'] = pd.read_csv(cpath)
    elif os.path.exists(cpath):
        os.remove(cpath)
    timings.to_csv(tpath, index=False)

    rpath = write_report(a, ctx, D, sp, dev, plan, hpo, cv, fin, checks, made, timings)
    log(f'figures: {len(made)} written{"; errors: " + str(figure_errors) if figure_errors else ""}')
    log(f'report: {rpath}')
    if not a.quick and os.path.realpath(a.out) == os.path.realpath(DEFAULT_RESULTS):
        # the generated accuracy block of README.md quotes these tables (tests/test_readme_numbers.py);
        # the figures quoted in the README prose are checked by that test, not rewritten
        try:
            from fermentation_insights.coeff_predictor import report_checks as RC
            log('README.md accuracy block '
                + ('updated' if RC.update_readme(os.path.join(PKG_DIR, 'README.md'), tdir) else 'unchanged'))
        except Exception as e:
            log(f'README.md accuracy block NOT updated ({type(e).__name__}: {e}); run '
                'python -m fermentation_insights.coeff_predictor.report_checks')
    log(f'===== done in {(time.time() - t_start) / 60:.1f} min (stages: '
        + ', '.join(f'{k} {v / 60:.1f}' for k, v in stage_t.items()) + ' min)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
