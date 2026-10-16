# -*- coding: utf-8 -*-
"""
Post-hoc corrections of recorded raw data (DESIGN.md, "Data corrections").

The adapter review (defects D*, minor findings M*) found errors in data that
were already recorded. Those whose correct value follows from recorded
columns are corrected here, when :func:`dataset.build` loads the raw JSONL
(``build(..., patches=True)``, the default). Forward fixes are versioned:

* harness fixes by ``harness_version`` (``sample_start`` record; v3 fixes D3);
* adapter fixes by ``adapter_version`` (configuration dict, written into the
  ``config`` record; v2 fixes D2/D4/D6/D7/D11/D12/M1/M3/M12/M23). A sample's
  configuration is the latest ``config`` record preceding its
  ``sample_start`` in its own file (:func:`dataset.read_raw`), so resumed files
  that mix versions are patched row by row.

A patch applies to rows of the affected biorefinery/configuration whose
version predates the fix (``harness_version`` < 3, ``adapter_version``
missing or < 2), and only where the recorded values carry the defect's
signature (e.g. the boiler gas is still in the material pool). Each patch is
idempotent: patched rows/records get the patch id in their ``data_patches``
field (';'-joined) and are skipped on a second pass. :func:`apply` returns
counts per patch (``report.json['data_patches']``).

Patches
-------
D3  (harness v2, ``collect_natural_gas``): the boiler's natural gas -- the
    ``BoilerTurbogenerator`` 'Fuel' stream utility, charged by the TEA as a
    utility at ``bst.stream_prices['Fuel']`` = ``BT.natural_gas_price`` -- was
    booked as an unpriced material feed. The boiler gas ``B`` is moved to the
    utility pool at the recorded ``natural_gas_price`` (= ``BT.natural_gas_price``
    in all affected adapters): ``ng_flow_utility += B``,
    ``ng_cost_utility += B x natural_gas_price``, ``ng_flow_material -= B``
    (``ng_cost_material`` unchanged: the gas was unpriced there). ``B`` =
    ``ng_flow_material`` (HP, oxalic, isobutanol: all their material NG is the
    boiler's; requires ``ng_flow_utility + ng_flow_material = natural_gas_flow``
    and ``ng_cost_material = 0``), ``x_natural_gas_BT`` (TAL, ethanol cornstover
    and sugarcane, succinic; requires ``B <= ng_flow_material``; succinic also
    ``ng_flow_material = x_natural_gas_BT + x_natural_gas_drying`` and
    ``ng_flow_utility = 0``). Equivalent to the reviewer's table (HP/oxalic/
    isobutanol: utility = ``natural_gas_flow`` at ``natural_gas_price``
    [oxalic: ``x_BT_natural_gas_cost``], material 0; succinic: material =
    ``x_natural_gas_drying``). Not affected: lactic (priced boiler gas),
    ethanol corn (no boiler; 'Natural gas' utilities were recognised).
D4  (ethanol_sugarcane): MPSP was net of the cane module's RIN credit
    ('Ethanol RIN D5'), recorded as ``coproduct_sales``. Applied only if the
    credit per kg of product stream, ``coproduct_sales/(operating_hours x
    x_F_product_stream)``, is constant over the rows (i.e. the coproduct sales
    are the credit alone; 0.2593 USD/kg in the recorded data):
    ``x_RIN_credit = coproduct_sales/operating_hours`` [USD/h],
    ``MPSP += coproduct_sales/(operating_hours x F_final_product)``,
    ``x_MPSP_denatured += credit/kg``, ``x_MESP_USD_per_gal_denatured +=
    credit/kg x kg/gal``, ``coproduct_sales = 0`` (the credit and the product
    price enter the DCF identically, so the solved price shifts exactly by it).
D6  (HP): ``base_cost`` used the intermediate ``fermentation_lime`` stream
    (price frozen at 0.286335); the TEA charges ``lime_fresh`` (linked to it,
    same flow) at the sampled ``p_base_price``: ``base_cost =
    x_lime_fermentation x p_base_price``.
D7  (HP): per-sample inflection: ``x_inflection_yield = (1 - p_cell_mass_yield)
    x theoretical_yield`` (``p_cell_mass_yield`` sets
    ``R302.regular_biomass_conversion``; the succinic/xylitol regular
    conversions are 0: config ``inflection_yield / theoretical_yield`` = 1 -
    baseline 0.245; checked per sample against ``x_cell_mass_conversion`` at
    low yield) and ``x_above_inflection = y > x_inflection_yield``.
D2  (ethanol_corn, adapter v1): molar yeast growth destroyed sugar mass; not
    patchable -- samples are flagged ``_unpatchable = 'D2'`` and rejected by
    :func:`dataset.build` (``unpatchable_D2``); the configuration was
    regenerated with adapter v2.
M1  config: TAL ``neutralizing_agent = 'NaOH'`` (pH-control base; its cost is
    ``base_cost``), succinic low-pH ``neutralizing_agent = 'Ca(OH)2'`` (partial
    pH control); ``neutralization`` unchanged (False).
M3  TAL: ``bp_fermentation_fermentation_y_lipolytica_yield`` ->
    ``p_cell_mass_yield`` (sample values and the config's parameter list).
M12 TAL ``product_water_solubility = 9.8`` g/L water at 25 C (experimental data
    of the TAL module; v1: 5.41 from the fitted model).
M23 HP non-neutral: ``separation_principle`` without the 'acidulation+' prefix.
M11 (resumed samples appending duplicate rows) needs no patch:
    :func:`dataset.read_raw` keeps only the points of the chosen (complete)
    attempt of each sample.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

__all__ = ('MARK', 'FIXED_HARNESS_VERSION', 'FIXED_ADAPTER_VERSION', 'apply',
           'sample_versions', 'POINT_PATCHES', 'CONFIG_PATCHES')

#: field holding the ';'-joined ids of the patches applied to a row/record
MARK = 'data_patches'
#: harness version that contains the D3 fix
FIXED_HARNESS_VERSION = 3
#: adapter version that contains the adapter-level fixes
FIXED_ADAPTER_VERSION = 2

POINT_PATCHES = ('D3', 'D4', 'D6', 'D7')
CONFIG_PATCHES = ('M1', 'M3', 'M12', 'M23')

KG_PER_GAL_ETHANOL = 3.78541 * 0.789
TAL_WATER_SOLUBILITY_25C = 9.8
TAL_OLD_CELL_MASS = 'bp_fermentation_fermentation_y_lipolytica_yield'
_RTOL = 1e-6


def _num(s):
    return pd.to_numeric(s, errors='coerce')


def _marks(v):
    return set(str(v).split(';')) - {'', 'nan', 'None'} if isinstance(v, str) else set()


def _join(marks):
    return ';'.join(sorted(marks))


def _version(v, default):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return default
    return v if np.isfinite(v) else default


def _close(a, b, scale=None):
    scale = np.maximum(np.abs(b) if scale is None else scale, 1.)
    return np.abs(a - b) <= _RTOL * scale


# %% Versions

def sample_versions(cfg, starts):
    """DataFrame indexed by sample_id: biorefinery, config_id,
    harness_version (0 if absent), adapter_version (1 if absent) and the
    sample's config record (dict)."""
    cfg_by_key = {(r.get('biorefinery'), r.get('config_id')): r
                  for r in (cfg.to_dict('records') if len(cfg) else [])}
    rows = []
    for srow in starts.to_dict('records'):
        sid = str(srow['sample_id'])
        parts = sid.rsplit('|', 2)[0].split('|', 1)
        bio, cid = (parts[0], parts[1]) if len(parts) == 2 else (parts[0], None)
        crow = srow.get('_config') if isinstance(srow.get('_config'), dict) else cfg_by_key.get((bio, cid))
        rows.append({'sample_id': srow['sample_id'], 'biorefinery': bio, 'config_id': cid,
                     'harness_version': _version(srow.get('harness_version'), 0.),
                     'adapter_version': _version((crow or {}).get('adapter_version'), 1.),
                     'config': crow})
    return pd.DataFrame(rows).set_index('sample_id') if rows else \
        pd.DataFrame(columns=['biorefinery', 'config_id', 'harness_version', 'adapter_version', 'config'])


# %% Config-record patches

def _patch_config(rec):
    """Patch one config record in place; returns the ids applied now."""
    if _version(rec.get('adapter_version'), 1.) >= FIXED_ADAPTER_VERSION:
        return set()
    bio, cid = rec.get('biorefinery'), str(rec.get('config_id'))
    done = _marks(rec.get(MARK))
    new = set()
    if bio == 'TAL':
        if 'M1' not in done and not rec.get('neutralizing_agent'):
            rec['neutralizing_agent'] = 'NaOH'
            new.add('M1')
        if 'M3' not in done and isinstance(rec.get('parameters'), str) \
                and TAL_OLD_CELL_MASS in rec['parameters']:
            try:
                params = json.loads(rec['parameters'])
                if not any(p[0] == 'p_cell_mass_yield' for p in params):
                    params = [['p_cell_mass_yield', *p[1:]] if p[0] == TAL_OLD_CELL_MASS else p
                              for p in params]
                    rec['parameters'] = json.dumps(params)
                    new.add('M3')
            except Exception:
                pass
        if 'M12' not in done and _num(pd.Series([rec.get('product_water_solubility')])).iloc[0] \
                != TAL_WATER_SOLUBILITY_25C:
            rec['product_water_solubility'] = TAL_WATER_SOLUBILITY_25C
            new.add('M12')
    elif bio == 'succinic':
        if 'M1' not in done and not rec.get('neutralizing_agent'):
            rec['neutralizing_agent'] = 'Ca(OH)2'
            new.add('M1')
    elif bio == 'HP':
        sp = rec.get('separation_principle')
        neutral = rec.get('neutralization')
        neutral = neutral is True or str(neutral).lower() in ('true', '1', '1.0')
        if 'M23' not in done and not neutral and isinstance(sp, str) and sp.startswith('acidulation+'):
            rec['separation_principle'] = sp[len('acidulation+'):]
            new.add('M23')
    if new:
        rec[MARK] = _join(done | new)
    return new


# %% Point patches

def _rows_without(points, mask, pid):
    if MARK in points.columns:
        mask = mask & ~points[MARK].map(lambda v: pid in _marks(v)).to_numpy()
    return mask


def _mark_rows(points, mask, pid):
    if MARK not in points.columns:
        points[MARK] = ''
    idx = points.index[mask]
    points.loc[idx, MARK] = [_join(_marks(v) | {pid}) for v in points.loc[idx, MARK]]


def _col(points, c):
    return _num(points[c]).to_numpy(float) if c in points.columns else np.full(len(points), np.nan)


def _d3_source(bio, cid):
    if bio in ('HP', 'oxalic', 'isobutanol'):
        return 'material'
    if bio in ('TAL', 'succinic') or (bio == 'ethanol' and cid in ('ethanol_cornstover', 'ethanol_sugarcane')):
        return 'x_natural_gas_BT'
    return None


def _patch_D3(points, hv, counts):
    bio = points['biorefinery'].astype(str).to_numpy()
    cid = points['config_id'].astype(str).to_numpy()
    src = [_d3_source(b, c) for b, c in zip(bio, cid)]
    has_src = np.array([v is not None for v in src])
    is_mat = np.array([v == 'material' for v in src])
    fu, cu = _col(points, 'ng_flow_utility'), _col(points, 'ng_cost_utility')
    fm, cm = _col(points, 'ng_flow_material'), _col(points, 'ng_cost_material')
    ngf, ngp = _col(points, 'natural_gas_flow'), _col(points, 'natural_gas_price')
    BT, dry = _col(points, 'x_natural_gas_BT'), _col(points, 'x_natural_gas_drying')
    cand = (hv < FIXED_HARNESS_VERSION) & has_src & np.isfinite(fu + cu + fm + cm + ngp)
    cand = _rows_without(points, cand, 'D3')
    B = np.where(is_mat, fm, BT)
    with np.errstate(invalid='ignore'):
        ok = cand & np.isfinite(B) & (B >= 0) & (B <= fm + _RTOL * np.maximum(fm, 1.))
        ok &= ~is_mat | (_close(fu + fm, ngf) & (np.abs(cm) <= 1e-9))
        is_suc = bio == 'succinic'
        ok &= ~is_suc | (_close(fm, BT + np.nan_to_num(dry)) & (np.abs(fu) <= 1e-9))
    B = np.minimum(np.where(ok, B, 0.), np.where(ok, fm, 0.))
    if ok.any():
        idx = points.index[ok]
        points.loc[idx, 'ng_flow_utility'] = (fu + B)[ok]
        points.loc[idx, 'ng_cost_utility'] = (cu + B * ngp)[ok]
        points.loc[idx, 'ng_flow_material'] = (fm - B)[ok]
        _mark_rows(points, ok, 'D3')
    _count(counts, 'D3', points, ok, skipped=cand & ~ok,
           moved_kg_per_h_total=float(B[ok].sum()))
    return ok


def _patch_D4(points, av, counts):
    m = (points['config_id'].astype(str).to_numpy() == 'ethanol_sugarcane') & \
        (points['biorefinery'].astype(str).to_numpy() == 'ethanol') & (av < FIXED_ADAPTER_VERSION)
    cop, h = _col(points, 'coproduct_sales'), _col(points, 'operating_hours')
    Ff, Fs = _col(points, 'F_final_product'), _col(points, 'x_F_product_stream')
    with np.errstate(invalid='ignore', divide='ignore'):
        m &= np.isfinite(cop + h + Ff + Fs) & (cop > 0) & (h > 0) & (Ff > 0) & (Fs > 0)
        m = _rows_without(points, m, 'D4')
        c = cop / (h * Fs)   # credit per kg product stream (USD/kg denatured)
    info = {}
    if m.any():
        cm = c[m]
        spread = float((cm.max() - cm.min()) / np.median(cm))
        info = dict(credit_USD_per_kg_product_stream=float(np.median(cm)), credit_rel_spread=spread)
        if spread > 1e-6:   # coproduct sales are not the credit alone: do not patch
            _count(counts, 'D4', points, np.zeros(len(points), bool), skipped=m,
                   reason='credit per kg not constant', **info)
            return m & False
        idx = points.index[m]
        points.loc[idx, 'x_RIN_credit'] = (cop / h)[m]
        points.loc[idx, 'MPSP'] = (_col(points, 'MPSP') + cop / (h * Ff))[m]
        for k, f in (('x_MPSP_denatured', 1.), ('x_MESP_USD_per_gal_denatured', KG_PER_GAL_ETHANOL)):
            if k in points.columns:
                points.loc[idx, k] = (_col(points, k) + c * f)[m]
        points.loc[idx, 'coproduct_sales'] = 0.
        _mark_rows(points, m, 'D4')
    _count(counts, 'D4', points, m, **info)
    return m


def _patch_D6(points, av, starts_by_sid, counts):
    m = (points['biorefinery'].astype(str).to_numpy() == 'HP') & (av < FIXED_ADAPTER_VERSION)
    pb = points['sample_id'].map(starts_by_sid.get('p_base_price', {})).pipe(_num).to_numpy(float)
    lime = _col(points, 'x_lime_fermentation')
    m &= np.isfinite(pb + lime)
    m = _rows_without(points, m, 'D6')
    if m.any():
        points.loc[points.index[m], 'base_cost'] = (lime * pb)[m]
        _mark_rows(points, m, 'D6')
    _count(counts, 'D6', points, m)
    return m


def _patch_D7(points, av, starts_by_sid, theo_by_sid, counts):
    m = (points['biorefinery'].astype(str).to_numpy() == 'HP') & (av < FIXED_ADAPTER_VERSION)
    pc = points['sample_id'].map(starts_by_sid.get('p_cell_mass_yield', {})).pipe(_num).to_numpy(float)
    theo = points['sample_id'].map(theo_by_sid).pipe(_num).to_numpy(float)
    y = _col(points, 'y')
    m &= np.isfinite(pc + theo + y)
    m = _rows_without(points, m, 'D7')
    infl = (1. - pc) * theo
    # basis check: at low yield the cell-mass conversion is the sampled value
    xcc = _col(points, 'x_cell_mass_conversion')
    low = m & np.isfinite(xcc) & (y < 0.5 * infl)
    bad = pd.Series(low & ~_close(xcc, pc, 1.)).groupby(points['sample_id'].to_numpy()).any()
    bad_sids = set(bad.index[bad.to_numpy()])
    skipped = m & points['sample_id'].isin(bad_sids).to_numpy()
    m &= ~skipped
    if m.any():
        idx = points.index[m]
        points.loc[idx, 'x_inflection_yield'] = infl[m]
        points.loc[idx, 'x_above_inflection'] = (y > infl + 1e-9)[m]
        _mark_rows(points, m, 'D7')
    with np.errstate(invalid='ignore'):
        n_above = int((m & (y > infl + 1e-9)).sum())
    _count(counts, 'D7', points, m, skipped=skipped, points_above_inflection=n_above)
    return m


def _count(counts, pid, points, mask, skipped=None, **extra):
    c = counts.setdefault(pid, {'points': 0, 'samples': 0, 'by_config': {}})
    mask = np.asarray(mask, bool)
    if mask.any():
        sub = points.loc[mask]
        c['points'] += int(mask.sum())
        c['samples'] += int(sub['sample_id'].nunique())
        for k, n in Counter(sub['biorefinery'].astype(str) + '|' + sub['config_id'].astype(str)).items():
            c['by_config'][k] = c['by_config'].get(k, 0) + int(n)
    if skipped is not None:
        c['points_skipped'] = c.get('points_skipped', 0) + int(np.asarray(skipped, bool).sum())
    c.update(extra)


# %% Entry point

def apply(cfg, starts, points):
    """Apply all patches to the tables returned by :func:`dataset.read_raw`.

    Returns ``(cfg, starts, points, counts)``: patched copies (config dicts
    referenced by ``starts['_config']`` are copied and patched too) and a
    dict patch id -> counts (``points``/``samples``/``by_config`` patched,
    ``points_skipped`` = rows of a predating version without the defect's
    signature; config patches: ``config_records`` and ``configs``)."""
    counts = {}
    cfg = cfg.copy()
    starts = starts.copy()
    points = points.copy()
    if not len(starts):
        return cfg, starts, points, counts
    # config records: the global table and each sample's own record
    cp = defaultdict(set)
    n_rec = Counter()
    if len(cfg):
        recs = cfg.to_dict('records')
        for r in recs:
            for pid in _patch_config(r):
                cp[pid].add(f"{r.get('biorefinery')}|{r.get('config_id')}")
                n_rec[pid] += 1
        cfg = pd.DataFrame(recs)
    if '_config' in starts.columns:
        copies = {}
        new_cfgs = []
        for c in starts['_config']:
            if isinstance(c, dict):
                if id(c) not in copies:
                    copies[id(c)] = dict(c)
                    for pid in _patch_config(copies[id(c)]):
                        cp[pid].add(f"{c.get('biorefinery')}|{c.get('config_id')}")
                        n_rec[pid] += 1
                new_cfgs.append(copies[id(c)])
            else:
                new_cfgs.append(c)
        starts['_config'] = new_cfgs
    for pid in CONFIG_PATCHES:
        if pid in cp or pid not in counts:
            counts[pid] = {'config_records': int(n_rec.get(pid, 0)), 'configs': sorted(cp.get(pid, ()))}
    ver = sample_versions(cfg, starts)
    smark = {sid: set() for sid in starts['sample_id']}
    for sid, crow in ver['config'].items():
        if isinstance(crow, dict):
            smark[sid] |= _marks(crow.get(MARK))
    # sample-level fields
    is_tal = (ver['biorefinery'] == 'TAL') & (ver['adapter_version'] < FIXED_ADAPTER_VERSION)
    if TAL_OLD_CELL_MASS in starts.columns:
        sids = set(ver.index[is_tal.to_numpy()])
        old = _num(starts[TAL_OLD_CELL_MASS])
        cur = _num(starts['p_cell_mass_yield']) if 'p_cell_mass_yield' in starts.columns \
            else pd.Series(np.nan, index=starts.index)
        m = starts['sample_id'].isin(sids) & old.notna() & cur.isna()
        if m.any():
            starts.loc[m, 'p_cell_mass_yield'] = old[m]
            starts.loc[m, TAL_OLD_CELL_MASS] = np.nan
            for sid in starts.loc[m, 'sample_id']:
                smark[sid].add('M3')
        counts['M3']['samples'] = int(m.sum())
    unpatchable = (ver['biorefinery'] == 'ethanol') & (ver['config_id'] == 'ethanol_corn') & \
        (ver['adapter_version'] < FIXED_ADAPTER_VERSION)
    up = set(ver.index[unpatchable.to_numpy()])
    starts['_unpatchable'] = ['D2' if sid in up else '' for sid in starts['sample_id']]
    counts['D2'] = {'samples_unpatchable': len(up)}
    # point-level fields
    if len(points):
        points = points.reset_index(drop=True)
        hv = points['sample_id'].map(ver['harness_version']).fillna(0.).to_numpy(float)
        av = points['sample_id'].map(ver['adapter_version']).fillna(1.).to_numpy(float)
        sb = {c: dict(zip(starts['sample_id'], starts[c]))
              for c in ('p_base_price', 'p_cell_mass_yield') if c in starts.columns}
        theo = {sid: (c or {}).get('theoretical_yield') if isinstance(c, dict) else np.nan
                for sid, c in ver['config'].items()}
        _patch_D3(points, hv, counts)
        _patch_D4(points, av, counts)
        _patch_D6(points, av, sb, counts)
        _patch_D7(points, av, sb, theo, counts)
        if MARK in points.columns:
            for sid, v in zip(points['sample_id'], points[MARK]):
                if sid in smark:
                    smark[sid] |= _marks(v)
    starts[MARK] = [_join(smark.get(sid, set())) for sid in starts['sample_id']]
    starts['adapter_version'] = starts['sample_id'].map(ver['adapter_version']).to_numpy()
    return cfg, starts, points, counts
