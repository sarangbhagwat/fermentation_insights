# -*- coding: utf-8 -*-
"""
Simulation checks of the training data and of the model's contextual
response (report-only: nothing here feeds the training data). One
configuration per call, in the biorefinery's own environment, with the same
PYTHONPATH as the campaigns (Bioindustrial-Park worktree + this repository)::

    python -m fermentation_insights.coeff_predictor.data_generation.sim_checks \
        --adapter HP --config HP_glucose --out checks_raw/HP_glucose.json

All uncertain parameters are at their baselines and the configuration's
baseline productivity is used. The system is loaded and simulated as the
campaigns of the current harness do: the BioSTEAM corrections of
``harness.apply_biosteam_fixes`` are applied first (harness v4: D17; recorded
as ``harness_version`` / ``harness_fixes``); the ``ic_biogas`` check toggles
them itself. Five checks:

``capacity``
    MPSP at the baseline yield and titer with ``p_feedstock_capacity``
    scaled by ``--capacity-mults`` (the training data sample each
    configuration's capacity only about +/-10-20 % around its baseline).
    Compared with the model by ``report_checks.capacity_response``.

``ic_biogas``
    Known data defect D17 (``biosteam_fixes``): every point of a 3 x 3
    titer-yield grid (15/50/85 % of the yield range; 15/50/85 % of the
    titer range below the configuration's maximum titer) plus the baseline
    point, simulated three times: ``as_is`` (BioSTEAM's original method, as
    the campaigns of harness <= 3), ``harness`` (the BioSTEAM corrections
    the current harness applies, ``harness.apply_biosteam_fixes``, as its
    campaigns simulate; recorded with ``harness_version`` and
    ``harness_fixes``) and ``fixed`` (the corrected
    ``InternalCirculationRx._run``). Since harness v4 ``harness`` equals
    ``fixed``: no double count in the campaigns. Skipped (recorded as such)
    when the system has no InternalCirculationRx with the lumped method.

``duplicate_inlets``
    Known data defect D16: units whose inlet list holds the same stream
    object twice (the stream then enters the unit twice and the unit
    creates mass). Where found, the same grid is simulated as in the
    campaigns and again with the duplicate references removed (in memory).

``unit_bounds``
    Known data defect D18: design results outside the unit class's own
    ``_bounds`` at the baseline point (BioSTEAM does not warn; e.g. its
    distillation columns are bounded at 3-24 ft diameter and 27-170 ft
    height), with the unit's installed cost. Where a column is wider than
    :data:`TRAY_D_MAX`, its trays are priced by Seider's exponential
    correlation ``412.7 exp(0.1482 D)`` per tray (``biosteam.units.design_tools
    .column_design``, fitted to about 2-16 ft) far beyond its range, so the
    same grid is simulated as in the campaigns (``as_is``) and with such a
    column's trays priced as N = ceil((D / 24 ft)^2) parallel tray stacks of
    the same total area (``fixed``, in memory: a sensitivity, not BioSTEAM's
    method), and the capacity response at the baseline yield and titer
    (:data:`COLUMN_CAPACITY_MULTS`) in both modes.

``tea_settings``
    No simulation: the TEA conventions of the loaded configuration
    (``bst.CE``, the TEA class, IRR, income tax, operating days, plant life
    (``duration``), depreciation, construction and startup schedule, working
    capital, financing, Lang factor, labor cost: :data:`TEA_ATTRS`) and
    which TEA settings among the contextual parameters (:data:`TEA_PARAMS`)
    the adapter samples. The others are fixed at the module's value and are
    NaN in the training data (e.g. the IRR of 0.15 of ethanol_corn and
    IBO_EtOH_corn), which the model's input does not carry.

``--collect DIR`` (no BioSTEAM needed) merges the per-configuration JSON
files into ``capacity_response.csv``, ``ic_biogas_points.csv`` /
``ic_biogas_bias.csv``, ``duplicate_inlets_points.csv`` /
``duplicate_inlets_bias.csv`` (per configuration: MPSP change at the
baseline and at the lowest-yield grid point, and the change of the refitted
coefficients, at the same recovery), ``unit_bounds.csv`` and
``column_trays_points.csv`` / ``column_trays_bias.csv`` /
``column_trays_capacity.csv``, ``tea_settings.csv`` in ``--csv-dir`` (and,
with ``--registry-dir``, the TEA settings as :data:`TEA_SETTINGS_JSON`,
which ``predict.py`` reads from the package's ``registry/``).
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
#: [ft] upper diameter bound of BioSTEAM's distillation columns
#: (``Distillation._bounds['Diameter']``); check ``unit_bounds``
TRAY_D_MAX = 24.
#: capacity multipliers of the ``unit_bounds`` check's tray-cost sensitivity
COLUMN_CAPACITY_MULTS = (0.85, 1.0, 1.15, 2.0)


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


def _out_of_bounds(system, total_installed):
    """Design results outside the unit class's ``_bounds`` (a bound named
    'Diameter' also covers 'Rectifier diameter' / 'Stripper diameter'):
    [dict(unit, unit_class, key, value, lo, hi, installed_cost,
    installed_share)]."""
    rows = []
    for u in system.units:
        b = getattr(u, '_bounds', None) or {}
        dr = getattr(u, 'design_results', None) or {}
        for name, lim in b.items():
            try:
                lo, hi = (float(v) for v in lim)
            except (TypeError, ValueError):
                continue
            for k, v in dr.items():
                if not isinstance(v, (int, float)) or not (k == name or k.lower().endswith(' ' + name.lower())):
                    continue
                if v > hi * (1 + 1e-4) or v < lo * (1 - 1e-4):
                    ic = float(getattr(u, 'installed_cost', float('nan')))
                    rows.append(dict(unit=u.ID, unit_class=type(u).__name__, key=k, value=float(v), lo=lo, hi=hi,
                                     side='above' if v > hi else 'below', installed_cost=ic,
                                     installed_share=ic / total_installed if total_installed else float('nan')))
    return rows


def _wide_columns(system, dmax=TRAY_D_MAX):
    """{unit ID: max tray diameter [ft]} of the units with tray costs whose
    (rectifier / stripper) diameter exceeds ``dmax``."""
    out = {}
    for u in system.units:
        cost = getattr(u, 'baseline_purchase_costs', None) or {}
        if not any(k.endswith('rays') for k in cost):        # 'Trays', 'Rectifier trays', 'Stripper trays'
            continue
        dr = getattr(u, 'design_results', None) or {}
        d = [float(v) for k, v in dr.items() if isinstance(v, (int, float)) and k.lower().endswith('diameter')]
        if d and max(d) > dmax:
            out[u.ID] = max(d)
    return out


def _parallel_trays(orig, dmax=TRAY_D_MAX):
    """Tray cost of a column wider than ``dmax`` as N = ceil((D / dmax)^2)
    parallel stacks of diameter D / sqrt(N) (same total area)."""
    import math

    def cost(N_T, Di):
        if Di <= dmax:
            return orig(N_T, Di)
        n = math.ceil((Di / dmax) ** 2)
        return n * orig(N_T, Di / math.sqrt(n))
    return cost


def _unit_bounds(adapter, dom, params, y0, t0, r0, capacity_mults=COLUMN_CAPACITY_MULTS):
    from biosteam.units import design_tools
    base = _simulate(adapter, y0, t0, r0)
    rec = dict(baseline=base)
    if base.get('status') != 'ok':
        rec['skipped'] = f"baseline point not simulated ({base.get('status')})"
        return rec
    tot = float(getattr(adapter.tea, 'installed_equipment_cost', float('nan')))
    rec['installed_equipment_cost'] = tot
    rec['out_of_bounds'] = _out_of_bounds(adapter.system, tot)
    wide = _wide_columns(adapter.system)
    rec['wide_columns'] = wide
    if not wide:
        rec['skipped'] = f'no column wider than {TRAY_D_MAX:g} ft'
        return rec
    orig = design_tools.compute_purchase_cost_of_trays
    cap = [p for p in params if p.name == 'p_feedstock_capacity']
    rec['units'] = [[u, f'D {d:.1f} ft'] for u, d in wide.items()]
    rec['capacity'] = []
    try:
        for mode in ('as_is', 'fixed'):
            design_tools.compute_purchase_cost_of_trays = orig if mode == 'as_is' else _parallel_trays(orig)
            rows = []
            for name, y, t in _grid(adapter, dom):
                row = _simulate(adapter, y, t, r0)
                row['point'] = name
                rows.append(row)
                print('unit_bounds', mode, name, row.get('MPSP'), row['status'], flush=True)
            rec[mode] = rows
            if cap:
                for m in capacity_mults:
                    cap[0].setter(cap[0].baseline * m)
                    row = _simulate(adapter, y0, t0, r0)
                    row.update(mode=mode, mult=float(m))
                    rec['capacity'].append(row)
                    print('unit_bounds capacity', mode, m, row.get('MPSP'), row['status'], flush=True)
                cap[0].setter(cap[0].baseline)
    finally:
        design_tools.compute_purchase_cost_of_trays = orig
        if cap:
            cap[0].setter(cap[0].baseline)
    return rec


#: TEA attributes recorded by the ``tea_settings`` check (where the TEA class
#: has them): the financial convention, the construction / startup schedule
#: and the capital-cost basis of the configuration's module
TEA_ATTRS = ('IRR', 'income_tax', 'operating_days', 'duration', 'depreciation', 'construction_schedule',
             'startup_months', 'startup_FOCfrac', 'startup_VOCfrac', 'startup_salesfrac', 'WC_over_FCI',
             'finance_interest', 'finance_years', 'finance_fraction', 'lang_factor', 'labor_cost')
#: contextual parameters of the training data that are TEA settings: their
#: TEA attribute (same unit)
TEA_PARAMS = {'p_IRR': 'IRR', 'p_income_tax': 'income_tax', 'p_operating_days': 'operating_days'}
#: file name of the collected TEA settings in ``--registry-dir``
TEA_SETTINGS_JSON = 'tea_settings.json'


def _tea_settings(adapter, params):
    """TEA conventions of the loaded configuration (uncertain parameters at
    their baselines): ``bst.CE``, the TEA class and :data:`TEA_ATTRS`, and
    which of :data:`TEA_PARAMS` the adapter samples (the others are fixed at
    the module's value and are NaN in the training data)."""
    import biosteam as bst
    tea = adapter.tea
    out = dict(CE=float(bst.CE), tea_class=type(tea).__name__)
    for a in TEA_ATTRS:
        try:
            v = getattr(tea, a)
        except Exception:
            continue
        if v is None:
            out[a] = None
        elif hasattr(v, 'tolist'):
            out[a] = v.tolist()
        elif isinstance(v, (list, tuple)):
            out[a] = [x if isinstance(x, (int, float)) else str(x) for x in v]
        elif isinstance(v, (bool, int, float, str)):
            out[a] = v
        else:
            out[a] = str(v)
    names = {p.name for p in params}
    out['sampled'] = sorted(p for p in TEA_PARAMS if p in names)
    out['fixed'] = {p: out.get(a) for p, a in TEA_PARAMS.items() if p not in names and out.get(a) is not None}
    return out


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
    from . import harness
    from .schema import HARNESS_VERSION
    fixes = harness.apply_biosteam_fixes(getattr(A, 'biosteam_fixes', harness.BIOSTEAM_FIXES))
    A.load(cfg)
    params = A.parameters()
    for p in params:
        p.setter(p.baseline)
    dom = A.domain()
    y0, t0, r0 = dom['baseline_yield'], dom['baseline_titer'], dom['baseline_productivity']
    rec = dict(adapter=adapter_name, config_id=config_id, domain={k: float(v) for k, v in dom.items()
                                                                if isinstance(v, (int, float))},
               versions=A.versions(), harness_version=HARNESS_VERSION, harness_fixes=fixes)
    if 'tea_settings' in checks:
        rec['tea_settings'] = _tea_settings(A, params)
        print('tea_settings', config_id, json.dumps(rec['tea_settings'], default=str), flush=True)
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
            for mode in ('as_is', 'harness', 'fixed'):
                if mode == 'as_is':
                    biosteam_fixes.remove()
                elif mode == 'harness':
                    biosteam_fixes.remove()
                    rec_ic['harness_version'] = HARNESS_VERSION
                    rec_ic['harness_fixes'] = harness.apply_biosteam_fixes()
                else:
                    biosteam_fixes.apply('ic_biogas')
                rows = []
                for name, y, t in pts:
                    row = _simulate(A, y, t, r0, ics)
                    row['point'] = name
                    rows.append(row)
                    print('ic_biogas', mode, config_id, name, row.get('MPSP'), row['status'], flush=True)
                rec_ic[mode] = rows
            biosteam_fixes.remove()
            harness.apply_biosteam_fixes(fixes.split(';') if fixes else ())   # back to the campaign state
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
    if 'unit_bounds' in checks:
        rec['unit_bounds'] = _unit_bounds(A, dom, params, y0, t0, r0)
    return rec


# %% collection (no BioSTEAM)

def collect(raw_dir, csv_dir, registry_dir=None):
    import numpy as np
    import pandas as pd
    from fermentation_insights.coeff_predictor import physics
    recs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(raw_dir, '*.json')))]
    cap_rows, ic_pts, ic_bias, dup_pts, dup_bias = [], [], [], [], []
    ub_rows, col_pts, col_bias, col_cap = [], [], [], []
    tea_rows, tea_reg = [], {}
    for rec in recs:
        cid = rec['config_id']
        ts = rec.get('tea_settings')
        if ts is not None:
            tea_reg[cid] = dict(ts, adapter=rec.get('adapter'))
            row = dict(config_id=cid, adapter=rec.get('adapter'))
            for k, v in ts.items():
                if k == 'fixed':
                    row[k] = '; '.join(f'{p}={x:g}' for p, x in v.items())
                elif k == 'duration' and isinstance(v, list) and len(v) == 2:
                    row['duration'] = f'{v[0]}-{v[1]}'
                    try:
                        row['plant_life_years'] = int(v[1]) - int(v[0])
                    except (TypeError, ValueError):
                        pass
                elif isinstance(v, list):
                    row[k] = ' '.join(f'{x:g}' if isinstance(x, float) else str(x) for x in v)
                else:
                    row[k] = v
            tea_rows.append(row)
        c = rec.get('capacity') or {}
        for row in c.get('points', ()):
            cap_rows.append(dict(config_id=cid, harness_version=rec.get('harness_version', 3),
                                 capacity_baseline=c['baseline'],
                                 **{k: row.get(k) for k in ('mult', 'p_feedstock_capacity', 'y', 't', 'r',
                                                            'status', 'MPSP', 'recovery', 'FCI')}))
        for key, pts_, bias_, extra in (('ic_biogas', ic_pts, ic_bias, ('biogas_kg_h', 'ic_mass_gain_kg_h')),
                                        ('duplicate_inlets', dup_pts, dup_bias, ('F_sugar_fed',))):
            if key in rec:
                _paired_effect(rec[key], cid, pts_, bias_, extra)
        ub = rec.get('unit_bounds')
        if ub is not None:
            for row in ub.get('out_of_bounds', ()):
                ub_rows.append(dict(config_id=cid, harness_version=rec.get('harness_version'), **row))
            if not ub.get('out_of_bounds'):
                ok = (ub.get('baseline') or {}).get('status') == 'ok'
                ub_rows.append(dict(config_id=cid, harness_version=rec.get('harness_version'), unit='', key='',
                                    note='no design result outside the unit bounds' if ok
                                    else ub.get('skipped', 'not checked')))
            _paired_effect(ub, cid, col_pts, col_bias, ('FCI',))
            for row in ub.get('capacity', ()):
                col_cap.append(dict(config_id=cid, **{k: row.get(k) for k in ('mode', 'mult', 'y', 't', 'r', 'status',
                                                                              'MPSP', 'FCI')}))
    if col_cap:
        cc = pd.DataFrame(col_cap)
        ref = cc[cc['mult'] == 1.].set_index(['config_id', 'mode'])
        key = list(zip(cc['config_id'], cc['mode']))
        cc['MPSP_rel_to_1x'] = cc['MPSP'].to_numpy() / ref['MPSP'].reindex(key).to_numpy()
        cc['FCI_rel_to_1x'] = cc['FCI'].to_numpy() / ref['FCI'].reindex(key).to_numpy()
        col_cap = cc.to_dict('records')
    os.makedirs(csv_dir, exist_ok=True)
    written = []
    for name, rows in (('capacity_response.csv', cap_rows), ('ic_biogas_points.csv', ic_pts),
                       ('ic_biogas_bias.csv', ic_bias), ('duplicate_inlets_points.csv', dup_pts),
                       ('duplicate_inlets_bias.csv', dup_bias), ('unit_bounds.csv', ub_rows),
                       ('column_trays_points.csv', col_pts), ('column_trays_bias.csv', col_bias),
                       ('column_trays_capacity.csv', col_cap), ('tea_settings.csv', tea_rows)):
        if rows:
            path = os.path.join(csv_dir, name)
            pd.DataFrame(rows).to_csv(path, index=False, float_format='%.6g')
            written.append(path)
    if tea_reg and registry_dir:
        # read by predict.py (TEA settings the model input does not carry)
        path = os.path.join(registry_dir, TEA_SETTINGS_JSON)
        note = ('TEA settings of each training configuration at its baseline (data_generation/sim_checks.py, '
                'check tea_settings); "fixed": contextual parameters the adapter does not sample (NaN in the '
                'training data) and their value in the module')
        with open(path, 'w') as fh:       # one line per configuration
            fh.write('{"note": ' + json.dumps(note) + ',\n "configs": {\n')
            fh.write(',\n'.join(f'  {json.dumps(k)}: {json.dumps(tea_reg[k])}' for k in sorted(tea_reg)))
            fh.write('\n }\n}\n')
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
    h = pd.DataFrame(ic['harness']).set_index('point') if ic.get('harness') else None
    if h is not None:
        for p, row in h.iterrows():
            ic_pts.append(dict(config_id=cid, mode='harness', point=p, **row.to_dict()))
    if ok.get('baseline', False):
        out.update(y_baseline=a.loc['baseline', 'y'], t_baseline=a.loc['baseline', 't'],
                   MPSP_baseline_as_is=a.loc['baseline', 'MPSP'],
                   dMPSP_baseline_pct=100 * rel['baseline'])
        for c in extra:
            if c in a and c in f:
                out[f'{c}_as_is'] = a.loc['baseline', c]
                out[f'{c}_fixed'] = f.loc['baseline', c]
            if h is not None and c in h and 'baseline' in h.index:
                out[f'{c}_harness'] = h.loc['baseline', c]
    if h is not None:
        # the current harness's campaign state vs. the corrected method
        okh = ok & (h['status'] == 'ok').reindex(ok.index, fill_value=False)
        out.update(harness_version=ic.get('harness_version'), harness_fixes=ic.get('harness_fixes', ''),
                   n_points_ok_harness=int(okh.sum()),
                   dMPSP_harness_vs_fixed_max_abs_pct=100 * float(
                       (h['MPSP'].reindex(ok.index) / f['MPSP'] - 1)[okh].abs().max()) if okh.any() else np.nan,
                   dMPSP_harness_vs_as_is_baseline_pct=100 * float(h.loc['baseline', 'MPSP'] / a.loc['baseline', 'MPSP'] - 1)
                   if okh.get('baseline', False) else np.nan)
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
    ap.add_argument('--registry-dir', help='--collect: also write the tea_settings records (if any) as '
                                           f'{TEA_SETTINGS_JSON} into this directory (the package registry/, '
                                           'read by predict.py)')
    args = ap.parse_args(argv)
    if args.collect:
        for p in collect(args.collect, args.csv_dir or args.collect, args.registry_dir):
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
