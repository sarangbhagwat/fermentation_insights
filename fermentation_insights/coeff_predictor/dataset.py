# -*- coding: utf-8 -*-
"""
Build the processed training tables (DESIGN.md section 7) from raw harness
JSON-lines files.

    data = build(['/home/user/work/data/raw'])     # -> ProcessedData
    data.save(out_dir)                            # samples/points/configs .csv.gz
    data = ProcessedData.load(out_dir)            # or load(out_dir)

CLI::

    python -m fermentation_insights.coeff_predictor.dataset \
        --raw /home/user/work/data/raw --out <dir>

Raw reading
-----------
* every ``*.jsonl`` below each raw dir (recursively) is read line by line;
  unparsable lines (e.g. a truncated last line of a running campaign) are
  skipped and counted;
* a point belongs to the latest ``sample_start`` of its ``sample_id`` in the
  same file (resumed campaigns re-run incomplete samples); an attempt is
  complete when its ``sample_end`` is seen. If a sample_id is complete in
  several files (e.g. a ``_smoke`` file and a production file with the same
  seed), non-smoke files win, then the most recently modified file;
* each sample uses the latest ``config`` record of its configuration that
  precedes its chosen attempt in that attempt's own file (a resumed campaign
  appends a new config record, possibly of a newer adapter version), else
  the last such record of the file; otherwise the record ranked like attempts (non-smoke first, then the most
  recently modified file). Keys whose records differ between files are
  counted (``stats['config_conflicts']``); samples whose chosen attempt comes
  from a smoke file while a production file has a config record for the
  same configuration are dropped (``smoke_superseded``);
* columns absent in some adapters simply become NaN; see
  :func:`physics.opex_channels` for how missing cost columns are imputed
  (0 -> they fall into the remainder channels; recorded in
  ``diag_imputed``);
* feedstock descriptors ``fs_*`` of a sample (:func:`feedstock_descriptors`):
  median over the sample's ok points (harness v2 describes the feedstock
  stream at every point), else the ``sample_start`` record (feedstock
  described after applying the sample's parameters), else the config's
  baseline description.

Data corrections (:mod:`data_patches`, ``build(..., patches=True)``): recorded
rows of harness/adapter versions that predate a fix are corrected right after
reading (counts in ``report['data_patches']``; samples carry
``harness_version``, ``adapter_version`` and the applied patch ids
``data_patches``); samples that cannot be corrected (ethanol_corn adapter v1,
D2) are rejected (``unpatchable_D2``).

Feasible points: ``status == 'ok'``, finite positive y, t, MPSP, recovery,
and ``y <= inflection_yield`` (fitted domain of the publication) with the
sample's OWN inflection yield (min of the adapter's ``x_inflection_yield``
over its points, else the config value; points flagged
``x_above_inflection`` are dropped too); the samples table stores this
sample-level ``inflection_yield`` (config value in
``inflection_yield_config``). ``g`` is the median recovery over all ok
points of the sample (before inflection clipping), as in the publication.
Natural gas: harness v2 samples use the recorded per-pool NG flows/costs
(``physics.natural_gas_pools``, ``diag_ng_pool = 'split'``); for legacy (v1)
samples the pool (``physics.detect_ng_pool``) is decided once per
configuration (majority of the per-sample votes, ``diag_ng_pool_vote``).
Fixed channel: the TEA labor cost recorded by harness v2 (``diag_fixed_source``
'labor' | 'labor_foc'), else the FOC intercept ('foc_intercept').
``report['per_biorefinery']`` (:func:`biorefinery_summary`) and
``report['rejections_by_biorefinery']`` summarize the build per biorefinery.

Targets (per retained sample; see :mod:`physics`)
-------------------------------------------------
``a, b, c, d`` (primed, x g), ``g`` (median recovery), ``A_mpsp..D_mpsp``
(unprimed), ``R2``, ``fit_mape`` (fraction, i.e. 0.01 = 1 %), ``fit_rmse``
(USD/kg), ``n_points``, ``ch__<channel>__<k>``, ``int__<fs|el|ng>__<k>``,
``phi`` [1/yr per USD of FCI, or of installed cost if FCI is unavailable --
see ``diag_capital_basis``], ``F_fixed`` [USD/yr], ``closure_rel``, ``sum_check``,
``channels_ok`` (bool: phi identified without fallback, closure below the
threshold, natural-gas pool admissible, fixed costs >= 0 and measured labor
not above FOC -- use it to mask channel losses), ``diag_*`` diagnostics.

Quality filters (rejection reason in ``report['rejections']``), in order:
``smoke_superseded``, ``incomplete``, ``unpatchable_<id>``, ``param_error``, ``too_few_points`` (< min_points ok
points), ``too_few_yields`` / ``too_few_titers`` (< 3 distinct values),
``missing_columns`` (no F_final_product / operating_hours -> no
decomposition), ``nonfinite_targets``, ``bad_g`` (g <= 0 or > 1.05 x
stoichiometric mass ratio), ``low_R2``, ``high_closure``.

Samples whose phi/F regression is degenerate are re-decomposed with a prior
phi = median phi of the identified samples of the same configuration (else
biorefinery, else all) and flagged ``diag_phi_flag = 'prior'``.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from collections import Counter, defaultdict
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import data_patches, downstream, physics

__all__ = ('DEFAULT_RAW_DIR', 'DEFAULT_PROCESSED_DIR', 'QualityFilters',
           'ProcessedData', 'read_raw', 'build', 'load', 'TOTAL_TARGETS',
           'CHANNEL_TARGETS', 'INTENSITY_TARGETS', 'SPEC_COLUMNS', 'feedstock_descriptors',
           'summary', 'biorefinery_summary', 'main')

DEFAULT_RAW_DIR = '/home/user/work/data/raw'
DEFAULT_PROCESSED_DIR = os.path.join(os.path.dirname(__file__), 'data', 'processed')

TOTAL_TARGETS = ('a', 'b', 'c', 'd', 'g')
CHANNEL_TARGETS = tuple(f'ch__{c}__{k}' for c in physics.CHANNELS for k in physics.COEFF_KEYS)
INTENSITY_TARGETS = tuple(f'int__{q}__{k}' for q in physics.INTENSITY_KEYS
                          for k in physics.COEFF_KEYS)

#: specification columns copied from the config record into ``samples``
SPEC_COLUMNS = (
    'feedstock', 'feedstock_type', 'fermentation_product', 'final_product',
    'neutralization', 'neutralizing_agent', 'n_acid_groups', 'pKa1',
    'final_product_stoich_mass_ratio', 'upgrading', 'separation_variant',
    'separation_principle', 'aerobic', 'fermentation_mode', 'T_fermentation',
    'product_water_solubility', 'theoretical_yield', 'inflection_yield',
    'y_min', 'y_max', 't_min', 't_max', 'baseline_productivity',
    'baseline_yield', 'baseline_titer', 'ferm_product_chemical',
    'final_product_chemical',
) + downstream.SAMPLE_FIELDS   # structured downstream specification (data_patches S1)
_DESC_PREFIXES = ('fp_', 'pp_', 'fs_')
_POINT_KEEP = ('sample_id', 'point_id', 'y', 't', 'r', 'MPSP', 'recovery',
               'F_final_product', 'F_ferm_product', 'F_sugar_fed', 'F_feedstock',
               'V_broth', 'V_fermenters', 'GWP', 'FEC')


@dataclass
class QualityFilters:
    min_points: int = 8
    min_yields: int = 3
    min_titers: int = 3
    min_R2: float = 0.9
    max_closure_rel: float = 0.05
    max_g_over_stoich: float = 1.05
    require_complete: bool = True
    clip_to_inflection: bool = True
    drop_superseded_smoke: bool = True


# %% Raw reading

def _iter_jsonl(path, stats):
    with open(path, 'r', errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except Exception:
                stats['bad_lines'] += 1
                continue
            if not isinstance(rec, dict):
                stats['bad_lines'] += 1
                continue
            yield rec


def _find_files(raw_dirs):
    files = []
    for d in ([raw_dirs] if isinstance(raw_dirs, str) else raw_dirs):
        if os.path.isfile(d):
            files.append(d)
        else:
            files.extend(glob.glob(os.path.join(d, '**', '*.jsonl'), recursive=True))
    return sorted(set(files))


#: config-record fields ignored when checking two records of a key for conflicts
_CONFIG_VOLATILE = ('load_seconds', 'versions', '_file', 'record', 'timestamp', 'time', 'seconds')


def _config_signature(rec):
    return json.dumps({k: v for k, v in rec.items() if k not in _CONFIG_VOLATILE},
                      sort_keys=True, default=str)


def _sid_key(sid):
    """(biorefinery, config_id) from a sample id ``bio|config|s<seed>|idx``."""
    parts = str(sid).rsplit('|', 2)[0].split('|', 1)
    return (parts[0], parts[1]) if len(parts) == 2 else (parts[0], None)


def read_raw(raw_dirs, require_complete=True):
    """Read raw JSONL files.

    Returns ``(configs, starts, points, stats)``: DataFrame of config records
    (one per biorefinery/config_id, ranked like sample attempts: non-smoke
    files first, then the most recently modified, then the later file),
    DataFrame of chosen sample attempts (sample_start fields + sample_end
    fields + ``_file``, ``_complete``, ``_config`` = the config record of the
    attempt's own file if it has one, ``_smoke_superseded`` = the attempt is
    from a smoke file although a production file has a config record for its
    configuration), DataFrame of their point records, and a stats dict
    (``config_conflicts``: keys whose config records differ between files)."""
    stats = Counter()
    files = _find_files(raw_dirs)
    stats['files'] = len(files)
    configs = {}       # key -> (rank, record)
    file_configs = {}  # (file, key) -> record
    conflicts = set()
    attempts = {}   # sid -> list of attempt dicts
    for fi, path in enumerate(files):
        mtime = os.path.getmtime(path)
        smoke = 'smoke' in os.path.basename(path).lower()
        current = {}  # sid -> attempt dict in this file
        current_cfg = {}  # key -> latest config record read so far in this file
        for rec in _iter_jsonl(path, stats):
            kind = rec.get('record')
            if kind == 'config':
                key = (rec.get('biorefinery'), rec.get('config_id'))
                rec = dict(rec, _file=path)
                file_configs[(path, key)] = rec
                current_cfg[key] = rec
                rank = (not smoke, mtime, fi)
                old = configs.get(key)
                if old is not None and _config_signature(old[1]) != _config_signature(rec):
                    conflicts.add(key)
                if old is None or rank >= old[0]:
                    configs[key] = (rank, rec)
            elif kind == 'sample_start':
                sid = rec.get('sample_id')
                if sid is None:
                    stats['bad_records'] += 1
                    continue
                if sid in current and not current[sid]['_complete']:
                    stats['restarted_attempts'] += 1
                att = {'start': rec, 'end': None, 'points': [], '_file': path,
                       '_complete': False, '_rank': (not smoke, mtime, fi),
                       '_config': current_cfg.get(_sid_key(sid))}
                current[sid] = att
                attempts.setdefault(sid, []).append(att)
            elif kind == 'point':
                sid = rec.get('sample_id')
                att = current.get(sid)
                if att is None:  # point without sample_start (e.g. file cut)
                    att = {'start': {'sample_id': sid}, 'end': None, 'points': [],
                           '_file': path, '_complete': False, '_rank': (not smoke, mtime, fi)}
                    current[sid] = att
                    attempts.setdefault(sid, []).append(att)
                    stats['orphan_points'] += 1
                att['points'].append(rec)
            elif kind == 'sample_end':
                sid = rec.get('sample_id')
                att = current.get(sid)
                if att is not None:
                    att['end'] = rec
                    att['_complete'] = True
            else:
                stats['unknown_records'] += 1
    prod_keys = {k for k, (rank, _) in configs.items() if rank[0]}
    rows, pts = [], []
    for sid, atts in attempts.items():
        complete = [a for a in atts if a['_complete']]
        if complete:
            if len({a['_file'] for a in complete}) > 1:
                stats['duplicate_samples'] += 1
            chosen = max(complete, key=lambda a: a['_rank'])
        else:
            chosen = atts[-1]
        start = {k: v for k, v in chosen['start'].items() if k != 'record'}
        end = chosen['end'] or {}
        row = dict(start)
        key = _sid_key(sid)
        row.update({'n_ok_reported': end.get('n_ok'), 'param_error': end.get('param_error') or '',
                    'sample_seconds': end.get('seconds'), '_file': chosen['_file'],
                    '_complete': chosen['_complete'],
                    '_config': chosen.get('_config') or file_configs.get((chosen['_file'], key)),
                    '_smoke_superseded': bool(not chosen['_rank'][0] and key in prod_keys)})
        rows.append(row)
        pts.extend(chosen['points'])
    starts = pd.DataFrame(rows)
    points = pd.DataFrame(pts)
    if len(points):
        points = points.drop(columns=['record'], errors='ignore')
    cfg = pd.DataFrame([rec for _, rec in configs.values()]).drop(columns=['record'], errors='ignore')
    stats['config_records'] = len(cfg)
    stats['samples_seen'] = len(starts)
    stats['points_seen'] = len(points)
    stats['config_conflicts'] = len(conflicts)
    stats = dict(stats)
    if conflicts:
        stats['config_conflict_keys'] = sorted(f'{b}|{c}' for b, c in conflicts)
    return cfg, starts, points, stats


# %% Building

@dataclass
class ProcessedData:
    samples: pd.DataFrame
    points: pd.DataFrame
    configs: pd.DataFrame
    report: dict = field(default_factory=dict)

    def save(self, out_dir=DEFAULT_PROCESSED_DIR):
        os.makedirs(out_dir, exist_ok=True)
        self.samples.to_csv(os.path.join(out_dir, 'samples.csv.gz'), index=False)
        self.points.to_csv(os.path.join(out_dir, 'points.csv.gz'), index=False)
        self.configs.to_csv(os.path.join(out_dir, 'configs.csv.gz'), index=False)
        with open(os.path.join(out_dir, 'report.json'), 'w') as fh:
            json.dump(self.report, fh, indent=1, default=str)
        return out_dir

    @classmethod
    def load(cls, out_dir=DEFAULT_PROCESSED_DIR):
        r = lambda n: pd.read_csv(os.path.join(out_dir, f'{n}.csv.gz'), low_memory=False)
        rep = {}
        p = os.path.join(out_dir, 'report.json')
        if os.path.exists(p):
            with open(p) as fh:
                rep = json.load(fh)
        return cls(samples=r('samples'), points=r('points'), configs=r('configs'), report=rep)


def load(out_dir=DEFAULT_PROCESSED_DIR) -> ProcessedData:
    return ProcessedData.load(out_dir)


def _num(s):
    return pd.to_numeric(s, errors='coerce')


def _n_distinct(v, rel=1e-6):
    v = np.sort(np.asarray(v, float))
    if not len(v):
        return 0
    return int(1 + np.sum(np.diff(v) > rel * np.maximum(np.abs(v[1:]), 1e-12)))


def _parse_sid(sid):
    parts = str(sid).rsplit('|', 2)
    seed = idx = np.nan
    if len(parts) == 3:
        try:
            seed = int(parts[1].lstrip('s'))
        except ValueError:
            pass
        try:
            idx = int(parts[2])
        except ValueError:
            pass
    return seed, idx


def _feasible(P, inflection, clip):
    st = P['status'] if 'status' in P.columns else pd.Series('ok', index=P.index)
    st = st.fillna('ok')
    m = (st == 'ok')
    for c in ('y', 't', 'MPSP', 'recovery'):
        v = _num(P[c]) if c in P.columns else pd.Series(np.nan, index=P.index)
        m &= np.isfinite(v) & (v > 0)
    n_beyond = 0
    if clip:
        beyond = pd.Series(False, index=P.index)
        if inflection is not None and np.isfinite(inflection):
            beyond |= _num(P['y']) > inflection * (1 + 1e-6)
        if 'x_above_inflection' in P.columns:   # flagged by the adapter itself
            beyond |= P['x_above_inflection'].map(lambda v: v is True or str(v).lower() in ('true', '1', '1.0'))
        beyond &= m
        n_beyond = int(beyond.sum())
        m &= ~beyond
    return m.to_numpy(), n_beyond


def _sample_inflection(P, config_value):
    """The sample's own inflection yield: min of the adapter-recorded
    ``x_inflection_yield`` (depends on the sampled byproduct/cell-mass
    yields) over its points, else the config (baseline) value."""
    if 'x_inflection_yield' in P.columns:
        v = _num(P['x_inflection_yield']).to_numpy(float)
        v = v[np.isfinite(v) & (v > 0)]
        if len(v):
            return float(v.min())
    return config_value


def _present(v):
    return v is not None and not (isinstance(v, float) and not np.isfinite(v))


def feedstock_descriptors(P_ok, srow, crow):
    """``fs_*`` descriptors of one sample: the median over its ok points
    (harness v2 records the feedstock stream at every point, i.e. the
    sample-specific composition), else the ``sample_start`` record (feedstock
    after applying the sample's parameters), else the config record.
    Non-numeric fields take the most frequent value."""
    out = {k: v for k, v in (crow or {}).items() if k.startswith('fs_')}
    out.update({k: v for k, v in (srow or {}).items() if k.startswith('fs_') and _present(v)})
    for c in P_ok.columns:
        if not c.startswith('fs_'):
            continue
        v = _num(P_ok[c])
        if v.notna().any():
            out[c] = float(v.median())
        else:
            s = P_ok[c].dropna()
            s = s[s.map(lambda x: isinstance(x, str) and x != '')]
            if len(s):
                out[c] = s.mode().iloc[0]
    return out


def _decompose(P, phi_prior, ng_pool, capital_basis='auto', g=None):
    return physics.decompose_sample(P, g=g, phi_prior=phi_prior, ng_pool=ng_pool,
                                    capital_basis=capital_basis)


def build(raw_dirs=DEFAULT_RAW_DIR, out_dir=None, filters: QualityFilters = None,
          ng_pool='auto', capital_basis='auto', verbose=False, patches=True) -> ProcessedData:
    """Read raw dirs, correct recorded data of older harness/adapter versions
    (:mod:`data_patches`, unless ``patches=False``), fit/decompose every
    sample, filter, and (optionally) save the three tables to ``out_dir``."""
    f = filters or QualityFilters()
    cfg, starts, points, stats = read_raw(raw_dirs, f.require_complete)
    patch_counts = {}
    if patches and len(starts):
        cfg, starts, points, patch_counts = data_patches.apply(cfg, starts, points)
    rejections = Counter()
    rej_bio = defaultdict(Counter)

    def _rej(reason, sid):
        rejections[reason] += 1
        rej_bio[str(_sid_key(sid)[0])][reason] += 1

    n_beyond_total = 0
    if not len(starts) or not len(points):
        empty = ProcessedData(pd.DataFrame(), pd.DataFrame(), cfg,
                              {'stats': stats, 'rejections': {}, 'n_retained': 0,
                               'data_patches': patch_counts})
        if out_dir:
            empty.save(out_dir)
        return empty
    cfg_by_key = {(r['biorefinery'], r['config_id']): r for r in cfg.to_dict('records')} if len(cfg) else {}
    groups = {sid: P for sid, P in points.groupby('sample_id', sort=False)}
    candidates = []   # (sid, start_row, P_feasible, config_row)
    for srow in starts.to_dict('records'):
        sid = srow['sample_id']
        P = groups.get(sid)
        if P is None or not len(P):
            _rej('no_points', sid)
            continue
        bio = P['biorefinery'].iloc[0] if 'biorefinery' in P.columns else str(sid).split('|')[0]
        cid = P['config_id'].iloc[0] if 'config_id' in P.columns else None
        # config record of the attempt's own file, else the ranked global one
        crow = srow.get('_config') if isinstance(srow.get('_config'), dict) else cfg_by_key.get((bio, cid))
        if crow is None:
            stats['samples_without_config'] = stats.get('samples_without_config', 0) + 1
            crow = {'biorefinery': bio, 'config_id': cid}
        if f.drop_superseded_smoke and srow.get('_smoke_superseded', False):
            _rej('smoke_superseded', sid)
            continue
        if f.require_complete and not srow.get('_complete', False):
            _rej('incomplete', sid)
            continue
        unpatchable = srow.get('_unpatchable')
        if isinstance(unpatchable, str) and unpatchable:
            _rej(f'unpatchable_{unpatchable}', sid)
            continue
        if srow.get('param_error'):
            _rej('param_error', sid)
            continue
        infl = crow.get('inflection_yield')
        infl = float(infl) if infl is not None and pd.notna(infl) else None
        infl = _sample_inflection(P, infl)
        mask, nb = _feasible(P, infl, f.clip_to_inflection)
        n_beyond_total += nb
        # recovery g as in the publication: median over all ok points
        # (before clipping to the fitted domain)
        g_all = physics.median_recovery(_num(P.loc[_feasible(P, None, False)[0], 'recovery']))
        Pf = P.loc[mask].copy()
        if len(Pf) < f.min_points:
            _rej('too_few_points', sid)
            continue
        if _n_distinct(_num(Pf['y'])) < f.min_yields:
            _rej('too_few_yields', sid)
            continue
        if _n_distinct(_num(Pf['t'])) < f.min_titers:
            _rej('too_few_titers', sid)
            continue
        need = [c for c in ('F_final_product', 'operating_hours')
                if c not in Pf.columns or not np.isfinite(_num(Pf[c])).all()
                or not (_num(Pf[c]) > 0).all()]
        if need:
            _rej('missing_columns', sid)
            continue
        candidates.append((sid, srow, Pf, crow, bio, cid, infl, g_all))

    # pass 1: decompose without prior
    dec = {}
    pool_of = {}
    for sid, srow, Pf, crow, bio, cid, infl, g_all in candidates:
        pool_of[sid] = ng_pool
        try:
            dec[sid] = _decompose(Pf, None, ng_pool, capital_basis, g_all)
        except Exception as e:  # pragma: no cover - defensive
            dec[sid] = e
    # legacy (v1) data only: natural-gas pool decided once per configuration
    # (majority of the per-sample votes; ties -> 'utility'), so that the
    # materials_other / utilities_other targets do not jump between samples
    # of one config. Harness v2 samples (pool 'split': NG cost recorded per
    # pool) neither vote nor are re-decomposed.
    if ng_pool == 'auto':
        votes = defaultdict(Counter)
        for sid, srow, Pf, crow, bio, cid, infl, g_all in candidates:
            d = dec[sid]
            if isinstance(d, physics.SampleDecomposition) and d.diagnostics['ng_pool'] in ('material', 'utility'):
                votes[(bio, cid)][d.diagnostics['ng_pool']] += 1
        for sid, srow, Pf, crow, bio, cid, infl, g_all in candidates:
            d = dec[sid]
            if not isinstance(d, physics.SampleDecomposition):
                continue
            vote = d.diagnostics['ng_pool']
            if vote in ('material', 'utility'):
                c = votes[(bio, cid)]
                pool_of[sid] = 'material' if c['material'] > c['utility'] else 'utility'
                if pool_of[sid] != vote:
                    dec[sid] = _decompose(Pf, None, pool_of[sid], capital_basis, g_all)
            dec[sid].diagnostics['ng_pool_vote'] = vote
    # phi priors from identified samples
    phis = defaultdict(list)
    for sid, srow, Pf, crow, bio, cid, infl, g_all in candidates:
        d = dec[sid]
        if isinstance(d, physics.SampleDecomposition) and d.diagnostics.get('phi_flag') == 'ok' \
                and d.closure_rel <= f.max_closure_rel:
            phis[('c', bio, cid)].append(d.phi)
            phis[('b', bio)].append(d.phi)
            phis[('all',)].append(d.phi)
    for sid, srow, Pf, crow, bio, cid, infl, g_all in candidates:
        d = dec[sid]
        if isinstance(d, physics.SampleDecomposition) and d.diagnostics.get('phi_flag') != 'ok':
            for key in (('c', bio, cid), ('b', bio), ('all',)):
                if phis.get(key):
                    d2 = _decompose(Pf, float(np.median(phis[key])), pool_of[sid], capital_basis, g_all)
                    d2.diagnostics['phi_prior_level'] = key[0]
                    d2.diagnostics['ng_pool_vote'] = d.diagnostics.get('ng_pool_vote', d2.diagnostics['ng_pool'])
                    dec[sid] = d2
                    break

    rows, prow_frames = [], []
    for sid, srow, Pf, crow, bio, cid, infl, g_all in candidates:
        d = dec[sid]
        if not isinstance(d, physics.SampleDecomposition):
            _rej('decomposition_error', sid)
            continue
        row = d.to_row()
        key_vals = [row[k] for k in ('a', 'b', 'c', 'd', 'g', 'R2')] + \
                   [row[c] for c in CHANNEL_TARGETS]
        if not np.all(np.isfinite(key_vals)):
            _rej('nonfinite_targets', sid)
            continue
        stoich = _num(pd.Series([crow.get('final_product_stoich_mass_ratio')])).iloc[0]
        if d.g <= 0 or (np.isfinite(stoich) and stoich > 0 and d.g > f.max_g_over_stoich * stoich):
            _rej('bad_g', sid)
            continue
        if d.fit.R2 < f.min_R2:
            _rej('low_R2', sid)
            continue
        if not d.closure_rel <= f.max_closure_rel:
            _rej('high_closure', sid)
            continue
        seed, idx = _parse_sid(sid)
        out = {'sample_id': sid, 'config_id': cid, 'biorefinery': bio,
               'sample_kind': srow.get('sample_kind', Pf['sample_kind'].iloc[0]
                                       if 'sample_kind' in Pf.columns else np.nan),
               'seed': seed, 'sample_index': idx, 'source_file': os.path.basename(srow['_file'])}
        for c in SPEC_COLUMNS:
            out[c] = crow.get(c, np.nan)
        # the sample's own inflection yield (fitted domain), config value kept
        out['inflection_yield_config'] = out['inflection_yield']
        out['inflection_yield'] = infl if infl is not None else np.nan
        for k, v in crow.items():
            if k.startswith(_DESC_PREFIXES):
                out[k] = v
        # the sample's own feedstock description (points > sample_start > config)
        Pall = groups[sid]
        out.update(feedstock_descriptors(Pall.loc[_feasible(Pall, None, False)[0]], srow, crow))
        stc = Counter(Pall['status'].fillna('ok')) if 'status' in Pall.columns else Counter()
        out.update(n_points_raw=len(Pall), n_ok_raw=stc.get('ok', 0),
                   n_infeasible=stc.get('infeasible', 0), n_failed=stc.get('failed', 0))
        out['r'] = srow.get('r', _num(Pf['r']).median() if 'r' in Pf.columns else np.nan)
        out['r_mult'] = _num(Pf['r_mult']).median() if 'r_mult' in Pf.columns else np.nan
        for k, v in srow.items():
            if k.startswith(('p_', 'bp_')):
                out[k] = v
        out['harness_version'] = srow.get('harness_version', np.nan)
        out['adapter_version'] = srow.get('adapter_version', np.nan)
        out['data_patches'] = srow.get('data_patches', '') or ''
        out.update(row)
        out['channels_ok'] = bool(d.diagnostics.get('phi_flag') == 'ok'
                                  and d.closure_rel <= f.max_closure_rel
                                  and d.diagnostics.get('ng_pool_admissible', True)
                                  and not d.diagnostics.get('F_negative', False)
                                  and not d.diagnostics.get('fixed_exceeds_foc', False))
        rows.append(out)
        keep = [c for c in _POINT_KEEP if c in Pf.columns]
        pp = Pf[keep].reset_index(drop=True)
        for c in keep:
            if c != 'sample_id':
                pp[c] = _num(pp[c])
        pp['MPSP_fit'] = physics.basis_matrix(pp['y'], pp['t']) @ d.fit.coef
        pp = pd.concat([pp, d.point_contrib.reset_index(drop=True)], axis=1)
        prow_frames.append(pp)
    samples = pd.DataFrame(rows)
    pts = pd.concat(prow_frames, ignore_index=True) if prow_frames else pd.DataFrame()

    # configs table
    configs = cfg.copy()
    if len(configs):
        configs = configs.drop(columns=['versions', '_file'], errors='ignore')
        n_raw = Counter()
        for srow in starts.to_dict('records'):
            n_raw[str(srow['sample_id']).rsplit('|', 2)[0]] += 1
        keys = configs['biorefinery'].astype(str) + '|' + configs['config_id'].astype(str)
        configs['n_samples_raw'] = [n_raw.get(k, 0) for k in keys]
        n_ret = Counter((samples['biorefinery'].astype(str) + '|' + samples['config_id'].astype(str))) \
            if len(samples) else Counter()
        configs['n_samples_retained'] = [n_ret.get(k, 0) for k in keys]

    report = {'stats': stats, 'rejections': dict(rejections),
              'points_beyond_inflection_dropped': n_beyond_total,
              'n_candidates': len(candidates), 'n_retained': len(samples),
              'filters': f.__dict__, 'ng_pool': ng_pool, 'capital_basis': capital_basis,
              'patches': bool(patches), 'data_patches': patch_counts}
    if len(samples):
        report.update(
            median_R2=float(samples['R2'].median()),
            median_fit_mape=float(samples['fit_mape'].median()),
            median_closure_rel=float(samples['closure_rel'].median()),
            max_sum_check=float(samples['sum_check'].max()),
            phi_flags=dict(Counter(samples['diag_phi_flag'])),
            ng_pools=dict(Counter(samples['diag_ng_pool'])),
            ng_pool_inadmissible=int((~samples['diag_ng_pool_admissible'].astype(bool)).sum()),
            F_negative=int(samples['diag_F_negative'].astype(bool).sum()),
            capital_bases=dict(Counter(samples['diag_capital_basis'])),
            channels_ok=int(samples['channels_ok'].sum()),
            median_phi=float(samples['phi'].median()),
            median_area_coverage=float(samples['diag_area_coverage'].median()),
            median_el_area_mismatch=float(samples['diag_el_area_mismatch_rel'].median()),
            min_materials_other_rel=float(samples['diag_materials_other_min_rel'].min()),
            min_utilities_other_rel=float(samples['diag_utilities_other_min_rel'].min()),
            imputed=dict(Counter(x for s in samples['diag_imputed'].fillna('')
                                 for x in str(s).split(';') if x)),
            per_config={f'{b}|{c}': int(n) for (b, c), n in
                        samples.groupby(['biorefinery', 'config_id']).size().items()},
        )
    report['rejections_by_biorefinery'] = {k: dict(v) for k, v in rej_bio.items()}
    report['per_biorefinery'] = biorefinery_summary(samples, pts, report['rejections_by_biorefinery'])
    data = ProcessedData(samples=samples, points=pts, configs=configs, report=report)
    if out_dir:
        data.save(out_dir)
    if verbose:
        print(summary(data))
    return data


def _counts(v):
    return {str(k): int(n) for k, n in Counter(v.fillna('')).items()}


def biorefinery_summary(samples, points=None, rejections_by_biorefinery=None) -> dict:
    """Per-biorefinery quality/physics diagnostics of the processed tables
    (dict biorefinery -> dict; also stored in ``report['per_biorefinery']``):
    retained/rejected counts, median R2 / fit MAPE / closure_rel, phi (median,
    range, max coefficient of variation over the baseline+productivity
    samples of a configuration = same parameters, different productivity),
    fixed channel (median ``F_fixed``, median ``F_fixed / labor_cost``,
    negative / above-FOC counts, source), fraction ``channels_ok``, natural
    gas (median price, priced NG per kg final product over the points, price
    source and mixed pools), capex weighting and the max channel-sum error."""
    rej = rejections_by_biorefinery or {}
    bios = sorted((set(samples['biorefinery'].astype(str)) if len(samples) else set()) | set(rej))
    out = {}
    for bio in bios:
        s = samples[samples['biorefinery'].astype(str) == bio] if len(samples) else samples
        r = {'n_kept': int(len(s)), 'rejected': dict(rej.get(bio, {}))}
        if len(s):
            med = lambda c: float(np.nanmedian(_num(s[c]))) if c in s and _num(s[c]).notna().any() else np.nan
            cv = []
            if 'sample_kind' in s:
                for _, gq in s[s['sample_kind'].isin(['baseline', 'productivity'])].groupby('config_id'):
                    if len(gq) > 1:
                        cv.append(float(np.std(gq['phi']) / abs(np.mean(gq['phi']))))
            with np.errstate(divide='ignore', invalid='ignore'):
                fl = _num(s['F_fixed']) / _num(s['diag_labor_cost']) if 'diag_labor_cost' in s else pd.Series(np.nan)
            r.update(
                R2=med('R2'), fit_mape=med('fit_mape'), closure_rel=med('closure_rel'),
                closure_rel_max=float(_num(s['closure_rel']).max()),
                phi=med('phi'), phi_min=float(s['phi'].min()), phi_max=float(s['phi'].max()),
                phi_cv_productivity=max(cv) if cv else np.nan,
                F_fixed=med('F_fixed'),
                fixed_over_labor=float(np.nanmedian(fl)) if fl.notna().any() else np.nan,
                F_negative=int(s['diag_F_negative'].astype(bool).sum()),
                fixed_exceeds_foc=int(s.get('diag_fixed_exceeds_foc', pd.Series(False)).astype(bool).sum()),
                fixed_source=_counts(s.get('diag_fixed_source', pd.Series(dtype=object))),
                channels_ok_frac=float(s['channels_ok'].astype(bool).mean()),
                ng_pool=_counts(s['diag_ng_pool']),
                ng_price=med('natural_gas_price'),
                ng_price_source=_counts(s.get('diag_ng_price_source', pd.Series(dtype=object))),
                ng_mixed_pool=_counts(s.get('diag_ng_mixed_pool', pd.Series(dtype=object))),
                ng_unpriced_flow_share=med('diag_ng_unpriced_flow_share'),
                capex_weights=_counts(s['diag_capex_weights']),
                foc_alloc=_counts(s.get('diag_foc_alloc', pd.Series(dtype=object))),
                sum_check_max=float(s['sum_check'].max()))
            if points is not None and len(points) and 'q__ng' in points:
                q = _num(points.loc[points['sample_id'].isin(s['sample_id']), 'q__ng'])
                r['ng_kg_per_kg'] = float(q.median()) if q.notna().any() else np.nan
        out[bio] = r
    return out


def summary(data: ProcessedData) -> str:
    rep = data.report
    L = []
    st = rep.get('stats', {})
    L.append(f"files={st.get('files', 0)} bad_lines={st.get('bad_lines', 0)} "
             f"samples_seen={st.get('samples_seen', 0)} points_seen={st.get('points_seen', 0)} "
             f"duplicates={st.get('duplicate_samples', 0)} restarted={st.get('restarted_attempts', 0)} "
             f"config_conflicts={st.get('config_conflicts', 0)}")
    L.append(f"candidates={rep.get('n_candidates', 0)} retained={rep.get('n_retained', 0)} "
             f"rejections={rep.get('rejections', {})} "
             f"points_beyond_inflection_dropped={rep.get('points_beyond_inflection_dropped', 0)}")
    s = data.samples
    if len(s):
        L.append('per biorefinery/config (retained / raw):')
        cfg = data.configs
        for _, r in cfg.iterrows() if len(cfg) else []:
            L.append(f"  {r['biorefinery']:>10s} {str(r['config_id']):<40s} "
                     f"{int(r.get('n_samples_retained', 0)):4d} / {int(r.get('n_samples_raw', 0)):4d}")
        L.append(f"median R2={rep['median_R2']:.5f}  median fit_mape={rep['median_fit_mape']:.4f}  "
                 f"median closure_rel={rep['median_closure_rel']:.2e}  "
                 f"max sum_check={rep['max_sum_check']:.1e}")
        L.append(f"phi: median={rep['median_phi']:.4f}/yr flags={rep['phi_flags']}  "
                 f"channels_ok={rep['channels_ok']}/{len(s)}  ng_pool={rep['ng_pools']}  "
                 f"capital_basis={rep.get('capital_bases')}  "
                 f"median area_coverage={rep['median_area_coverage']:.4f}")
        L.append(f"median |net el (columns) - net el (areas)|/consumption="
                 f"{rep['median_el_area_mismatch']:.2e}  min materials_other/R="
                 f"{rep['min_materials_other_rel']:.3g}  min utilities_other/R="
                 f"{rep['min_utilities_other_rel']:.3g}")
        if rep.get('imputed'):
            L.append(f"imputed (missing->0) columns: {rep['imputed']}")
        pb = rep.get('per_biorefinery') or biorefinery_summary(s, data.points, rep.get('rejections_by_biorefinery'))
        L.append('per biorefinery: kept rejected | R2 fit_mape closure | phi (cv prod.) | F_fixed F/labor | '
                 'ok-frac | NG price, kg/kg, source | capex weights')
        for b_, r in pb.items():
            if not r['n_kept']:
                L.append(f"  {b_:>10s} kept=0 rejected={r['rejected']}")
                continue
            L.append(f"  {b_:>10s} {r['n_kept']:3d} {r['rejected']} | {r['R2']:.4f} {r['fit_mape']:.4f} "
                     f"{r['closure_rel']:.2e} | {r['phi']:.4f} ({r['phi_cv_productivity']:.1e}) | "
                     f"{r['F_fixed']:.3g} {r['fixed_over_labor']:.3g} | {r['channels_ok_frac']:.2f} | "
                     f"{r['ng_price']:.4g} {r.get('ng_kg_per_kg', np.nan):.3g} {r['ng_price_source']} "
                     f"{r['ng_mixed_pool']} | {r['capex_weights']} {r['foc_alloc']}")
        L.append('median total coefficients (primed) by config:')
        med = s.groupby(['biorefinery', 'config_id'])[['a', 'b', 'c', 'd', 'g', 'R2']].median()
        for (b, c), r in med.iterrows():
            L.append(f"  {b}|{c}: a={r['a']:.4g} b={r['b']:.4g} c={r['c']:.4g} d={r['d']:.4g} "
                     f"g={r['g']:.4f} R2={r['R2']:.4f}")
    return '\n'.join(L)


def main(argv=None):
    ap = argparse.ArgumentParser(description='Build processed coeff_predictor tables.')
    ap.add_argument('--raw', nargs='+', default=[DEFAULT_RAW_DIR])
    ap.add_argument('--out', default=DEFAULT_PROCESSED_DIR)
    ap.add_argument('--min-points', type=int, default=8)
    ap.add_argument('--min-R2', type=float, default=0.9)
    ap.add_argument('--max-closure', type=float, default=0.05)
    ap.add_argument('--include-incomplete', action='store_true')
    ap.add_argument('--no-inflection-clip', action='store_true')
    ap.add_argument('--ng-pool', default='auto', choices=('auto', 'split', 'utility', 'material'))
    ap.add_argument('--capital-basis', default='auto', choices=('auto', 'FCI', 'installed_cost'))
    ap.add_argument('--no-patches', action='store_true',
                    help='do not correct recorded data of older harness/adapter versions')
    a = ap.parse_args(argv)
    f = QualityFilters(min_points=a.min_points, min_R2=a.min_R2, max_closure_rel=a.max_closure,
                       require_complete=not a.include_incomplete,
                       clip_to_inflection=not a.no_inflection_clip)
    data = build(a.raw, out_dir=a.out, filters=f, ng_pool=a.ng_pool,
                 capital_basis=a.capital_basis, patches=not a.no_patches)
    print(summary(data))
    print(f'-> {a.out}')
    return data


if __name__ == '__main__':
    main()
