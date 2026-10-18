# -*- coding: utf-8 -*-
"""
Simulation checks of the training data and of the model's contextual
response (report-only: nothing here feeds the training data). One
configuration per call, in the biorefinery's own environment, with the same
PYTHONPATH as the campaigns (Bioindustrial-Park worktree + this repository)::

    python -m fermentation_insights.coeff_predictor.data_generation.sim_checks \
        --adapter HP --config HP_glucose --out checks_raw/HP_glucose.json

All uncertain parameters are at their baselines and the configuration's
baseline productivity is used. Two checks:

``capacity``
    MPSP at the baseline yield and titer with ``p_feedstock_capacity``
    scaled by ``--capacity-mults`` (the training data sample each
    configuration's capacity only about +/-10-20 % around its baseline).
    Compared with the model by ``report_checks.capacity_response``.

``ic_biogas``
    Known data defect D17 (``biosteam_fixes``): every point of a 3 x 3
    titer-yield grid (15/50/85 % of the yield range; 15/50/85 % of the
    titer range below the configuration's maximum titer) plus the baseline
    point, simulated as in the campaigns and again with the corrected
    ``InternalCirculationRx._run``. Skipped (recorded as such) when the
    system has no InternalCirculationRx with the lumped method.

``duplicate_inlets``
    Known data defect D16: units whose inlet list holds the same stream
    object twice (the stream then enters the unit twice and the unit
    creates mass). Where found, the same grid is simulated as in the
    campaigns and again with the duplicate references removed (in memory).

``--collect DIR`` (no BioSTEAM needed) merges the per-configuration JSON
files into ``capacity_response.csv``, ``ic_biogas_points.csv`` /
``ic_biogas_bias.csv`` and ``duplicate_inlets_points.csv`` /
``duplicate_inlets_bias.csv`` (per configuration: MPSP change at the
baseline and at the lowest-yield grid point, and the change of the refitted
coefficients, at the same recovery) in ``--csv-dir``.
``scripts/run_sim_checks.sh`` runs every configuration.
"""
from __future__ import annotations

import argparse
import glob
import importlib
import json
import os
import sys
import time
import traceback
import warnings

FRACS = (0.15, 0.5, 0.85)
CAPACITY_MULTS = (0.5, 0.85, 1.0, 1.15, 2.0, 3.0)


def _grid(adapter, dom):
    pts = [('baseline', dom['baseline_yield'], dom['baseline_titer'])]
    for i, fy in enumerate(FRACS):
        y = dom['y_min'] + fy * (dom['y_max'] - dom['y_min'])
        t_hi = dom['t_max']
        try:
            mt = adapter.max_titer(y)
        except Exception:
            mt = None
        if mt is not None and mt == mt:
            t_hi = min(t_hi, mt)
        for j, ft in enumerate(FRACS):
            t = dom['t_min'] + ft * (t_hi - dom['t_min'])
            pts.append((f'grid_{i}{j}', float(y), float(t)))
    return pts


def _ic_units(adapter):
    try:
        from biosteam.wastewater.high_rate import InternalCirculationRx as IC
    except Exception:
        return []
    return [u for u in adapter.system.units if isinstance(u, IC)]


def _duplicate_inlets(system):
    """[(unit, [indices of repeated inlet references])]."""
    out = []
    for u in system.units:
        seen, dup = set(), []
        for i, s in enumerate(u.ins):
            if id(s) in seen:
                dup.append(i)
            seen.add(id(s))
        if dup:
            out.append((u, dup))
    return out


def _remove_duplicate_inlets(found):
    for u, dup in found:
        for i in reversed(dup):
            s = u.ins[i]
            u.ins._streams.pop(i)
            s._sink = u          # still the unit's inlet (the first reference)


def _simulate(adapter, y, t, r, ics=()):
    from .harness import InfeasibleSpecification, collect_natural_gas
    out = dict(y=float(y), t=float(t), r=float(r))
    t0 = time.time()
    try:
        adapter.set_try(y, t, r)
        adapter.simulate()
        res = adapter.collect()
        out.update(status='ok', MPSP=float(res['MPSP']),
                   recovery=float(res.get('recovery', float('nan'))),
                   F_sugar_fed=float(res.get('F_sugar_fed', float('nan'))),
                   F_feedstock=float(res.get('F_feedstock', float('nan'))))
        tea = adapter.tea
        out['FCI'] = float(getattr(tea, 'FCI', float('nan')))
        out.update({k: float(v) for k, v in collect_natural_gas(adapter.system).items()})
        if ics:
            out['biogas_kg_h'] = float(sum(u.outs[0].F_mass for u in ics))
            out['biogas_CH4_kg_h'] = float(sum(u.outs[0].imass['CH4'] for u in ics
                                               if 'CH4' in u.outs[0].chemicals))
            out['biogas_LHV_MW'] = float(sum(u.outs[0].LHV for u in ics) / 3.6e6)
            out['ic_mass_gain_kg_h'] = float(sum(
                sum(s.F_mass for s in u.outs) - sum(s.F_mass for s in u.ins) for u in ics))
    except InfeasibleSpecification as e:
        out.update(status='infeasible', error=str(e)[:200])
    except Exception as e:
        out.update(status='failed', error=f'{type(e).__name__}: {e}'[:200])
        try:
            adapter.reset()
        except Exception:
            pass
    out['seconds'] = round(time.time() - t0, 2)
    return out


def run_config(adapter_name, config_id, checks, capacity_mults=CAPACITY_MULTS):
    warnings.filterwarnings('ignore')
    mod = importlib.import_module(
        f'fermentation_insights.coeff_predictor.data_generation.adapters.{adapter_name}')
    A = mod.Adapter()
    cfg = [c for c in A.configurations() if c['config_id'] == config_id][0]
    A.load(cfg)
    params = A.parameters()
    for p in params:
        p.setter(p.baseline)
    dom = A.domain()
    y0, t0, r0 = dom['baseline_yield'], dom['baseline_titer'], dom['baseline_productivity']
    rec = dict(adapter=adapter_name, config_id=config_id, domain={k: float(v) for k, v in dom.items()
                                                                if isinstance(v, (int, float))},
               versions=A.versions())
    if 'capacity' in checks:
        cap = [p for p in params if p.name == 'p_feedstock_capacity']
        if not cap:
            rec['capacity'] = dict(skipped='no p_feedstock_capacity parameter')
        else:
            cap = cap[0]
            rows = []
            for m in capacity_mults:
                cap.setter(cap.baseline * m)
                row = _simulate(A, y0, t0, r0)
                row.update(mult=float(m), p_feedstock_capacity=float(cap.baseline * m))
                rows.append(row)
                print('capacity', config_id, m, row.get('MPSP'), row['status'], flush=True)
            cap.setter(cap.baseline)
            rec['capacity'] = dict(baseline=float(cap.baseline), points=rows)
    if 'ic_biogas' in checks:
        ics = _ic_units(A)
        lumped = [u for u in ics if str(getattr(u, 'method', 'lumped')).lower() != 'separate']
        rec_ic = dict(units=[(u.ID, str(getattr(u, 'method', ''))) for u in ics])
        if not lumped:
            rec_ic['skipped'] = 'no InternalCirculationRx with the lumped method'
        else:
            from . import biosteam_fixes
            pts = _grid(A, dom)
            for mode in ('as_is', 'fixed'):
                if mode == 'fixed':
                    biosteam_fixes.apply('ic_biogas')
                rows = []
                for name, y, t in pts:
                    row = _simulate(A, y, t, r0, ics)
                    row['point'] = name
                    rows.append(row)
                    print('ic_biogas', mode, config_id, name, row.get('MPSP'), row['status'], flush=True)
                rec_ic[mode] = rows
            biosteam_fixes.remove('ic_biogas')
        rec['ic_biogas'] = rec_ic
    if 'duplicate_inlets' in checks:
        found = _duplicate_inlets(A.system)
        rec_d = dict(units=[(u.ID, [u.ins[i].ID for i in dup]) for u, dup in found])
        if not found:
            rec_d['skipped'] = 'no unit with a repeated inlet'
        else:
            pts = _grid(A, dom)
            for mode in ('as_is', 'fixed'):
                if mode == 'fixed':
                    _remove_duplicate_inlets(found)
                rows = []
                for name, y, t in pts:
                    row = _simulate(A, y, t, r0)
                    row['point'] = name
                    rows.append(row)
                    print('duplicate_inlets', mode, config_id, name, row.get('MPSP'), row['status'], flush=True)
                rec_d[mode] = rows
        rec['duplicate_inlets'] = rec_d
    return rec


# %% collection (no BioSTEAM)

def collect(raw_dir, csv_dir):
    import numpy as np
    import pandas as pd
    from fermentation_insights.coeff_predictor import physics
    recs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(raw_dir, '*.json')))]
    cap_rows, ic_pts, ic_bias, dup_pts, dup_bias = [], [], [], [], []
    for rec in recs:
        cid = rec['config_id']
        c = rec.get('capacity') or {}
        for row in c.get('points', ()):
            cap_rows.append(dict(config_id=cid, capacity_baseline=c['baseline'],
                                 **{k: row.get(k) for k in ('mult', 'p_feedstock_capacity', 'y', 't', 'r',
                                                            'status', 'MPSP', 'recovery', 'FCI')}))
        for key, pts_, bias_, extra in (('ic_biogas', ic_pts, ic_bias, ('biogas_kg_h', 'ic_mass_gain_kg_h')),
                                        ('duplicate_inlets', dup_pts, dup_bias, ('F_sugar_fed',))):
            if key in rec:
                _paired_effect(rec[key], cid, pts_, bias_, extra)
    os.makedirs(csv_dir, exist_ok=True)
    written = []
    for name, rows in (('capacity_response.csv', cap_rows), ('ic_biogas_points.csv', ic_pts),
                       ('ic_biogas_bias.csv', ic_bias), ('duplicate_inlets_points.csv', dup_pts),
                       ('duplicate_inlets_bias.csv', dup_bias)):
        if rows:
            path = os.path.join(csv_dir, name)
            pd.DataFrame(rows).to_csv(path, index=False, float_format='%.6g')
            written.append(path)
    return written


def _paired_effect(ic, cid, ic_pts, ic_bias, extra):
    """Rows of one configuration's as-is / corrected check (``ic``) for the
    points and the per-configuration effect tables."""
    import numpy as np
    import pandas as pd
    from fermentation_insights.coeff_predictor import physics
    units = ';'.join(str(u[0]) for u in ic.get('units', ()))
    if 'as_is' not in ic:
        ic_bias.append(dict(config_id=cid, units=units, affected=False, note=ic.get('skipped', '')))
        return
    a = pd.DataFrame(ic['as_is']).set_index('point')
    f = pd.DataFrame(ic['fixed']).set_index('point')
    for mode, d in (('as_is', a), ('fixed', f)):
        for p, row in d.iterrows():
            ic_pts.append(dict(config_id=cid, mode=mode, point=p, **row.to_dict()))
    ok = (a['status'] == 'ok') & (f['status'] == 'ok')
    rel = (f['MPSP'] / a['MPSP'] - 1)[ok]
    out = dict(config_id=cid, units=units, affected=True, n_points_ok=int(ok.sum()))
    if ok.get('baseline', False):
        out.update(y_baseline=a.loc['baseline', 'y'], t_baseline=a.loc['baseline', 't'],
                   MPSP_baseline_as_is=a.loc['baseline', 'MPSP'],
                   dMPSP_baseline_pct=100 * rel['baseline'])
        for c in extra:
            if c in a and c in f:
                out[f'{c}_as_is'] = a.loc['baseline', c]
                out[f'{c}_fixed'] = f.loc['baseline', c]
    grid = [p for p in rel.index if p.startswith('grid_')]
    if grid:
        ylo = min(a.loc[grid, 'y'])
        low = [p for p in grid if a.loc[p, 'y'] == ylo]
        out.update(y_low=ylo, dMPSP_low_yield_pct=100 * float(rel[low].mean()),
                   dMPSP_grid_min_pct=100 * float(rel[grid].min()),
                   dMPSP_grid_max_pct=100 * float(rel[grid].max()))
    if len(grid) >= 6:
        g = physics.median_recovery(a.loc[grid, 'recovery'].to_numpy(float))
        fa = physics.fit_ols(a.loc[grid, 'y'], a.loc[grid, 't'], a.loc[grid, 'MPSP'])
        ff = physics.fit_ols(f.loc[grid, 'y'], f.loc[grid, 't'], f.loc[grid, 'MPSP'])
        for k, ca, cf in zip(('a', 'b', 'c', 'd'), g * fa.coef, g * ff.coef):
            out[f'{k}_as_is'] = ca
            out[f'{k}_fixed'] = cf
            out[f'd{k}_pct'] = 100 * (cf / ca - 1) if ca != 0 else np.nan
        out['fit_R2_as_is'] = fa.R2
        out['fit_R2_fixed'] = ff.R2
    ic_bias.append(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--adapter')
    ap.add_argument('--config')
    ap.add_argument('--out', help='JSON file for one configuration')
    ap.add_argument('--checks', nargs='*', default=['capacity', 'ic_biogas'])
    ap.add_argument('--capacity-mults', type=float, nargs='*', default=list(CAPACITY_MULTS))
    ap.add_argument('--collect', help='directory of per-configuration JSON files')
    ap.add_argument('--csv-dir', help='output directory of --collect')
    args = ap.parse_args(argv)
    if args.collect:
        for p in collect(args.collect, args.csv_dir or args.collect):
            print('wrote', p)
        return 0
    if not (args.adapter and args.config and args.out):
        ap.error('--adapter, --config and --out are required (or --collect)')
    t0 = time.time()
    try:
        rec = run_config(args.adapter, args.config, args.checks, args.capacity_mults)
    except Exception as e:
        traceback.print_exc()
        rec = dict(adapter=args.adapter, config_id=args.config,
                   error=f'{type(e).__name__}: {e}'[:500])
    rec['seconds'] = round(time.time() - t0, 1)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, 'w') as fh:
        json.dump(rec, fh, indent=1, default=float)
    print('wrote', args.out, rec['seconds'], 's')
    return 0


if __name__ == '__main__':
    sys.exit(main())
