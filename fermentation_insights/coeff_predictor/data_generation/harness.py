# -*- coding: utf-8 -*-
"""
Generic simulation harness for generating training data for the coefficient
predictor.

The harness is run *inside* a biorefinery-specific conda environment (it only
depends on numpy, pandas and biosteam/thermosteam). For each "sample" (one
configuration x one set of technological/contextual parameter values x one
productivity) it simulates a small titer-yield design, recording per-point
system-level (TEA, per-process-area costs and utilities) and process-level
(mass and energy flows) data in JSON-lines files with the canonical schema
defined in :mod:`schema`.

Coefficient fitting and all machine learning happen downstream, outside the
simulation environments.

Usage (from a biorefinery env)::

    python -m fermentation_insights.coeff_predictor.data_generation.run \
        --adapter HP --out <dir> --n-random 40 --seed 1
"""
from __future__ import annotations

import json
import math
import os
import platform
import signal
import subprocess
import time
import traceback
from dataclasses import dataclass, field
from typing import Callable, Iterable, Optional

import numpy as np

from .schema import CANONICAL_AREAS, AREA_QUANTITIES

__all__ = (
    'InfeasibleSpecification', 'Parameter', 'BiorefineryAdapter',
    'run_campaign', 'collect_area_metrics', 'collect_tea_metrics',
    'describe_chemical', 'describe_feedstock', 'latin_hypercube',
)


class InfeasibleSpecification(Exception):
    """Raised by adapters when a TRY combination is physically infeasible
    (e.g., the sugar concentration needed for the titer cannot be reached)."""


class PointTimeout(Exception):
    pass


# %% Parameters

@dataclass
class Parameter:
    """An uncertain technological/contextual parameter.

    Parameters
    ----------
    name : str
        Column name. Use a canonical ``p_*`` name (see schema) when the
        parameter has a clear cross-biorefinery meaning, else ``bp_*``.
    baseline : float
        Baseline value.
    inverse_cdf : callable
        Maps a uniform variate u in (0, 1) to a parameter value.
    setter : callable
        Applies a value to the loaded system.
    description : str
        Original parameter name/units from the biorefinery model.
    """
    name: str
    baseline: float
    inverse_cdf: Callable[[float], float]
    setter: Callable[[float], None]
    description: str = ''

    @classmethod
    def from_biosteam(cls, parameter, name=None):
        """Create from a biosteam.evaluation.Parameter with a chaospy
        distribution."""
        dist = parameter.distribution
        return cls(
            name=name or 'bp_' + _slug(parameter.name),
            baseline=float(parameter.baseline),
            inverse_cdf=lambda u, d=dist: float(d.inv(u)),
            setter=parameter.setter,
            description=f'{parameter.element_name}: {parameter.name} [{parameter.units}]'
                        if hasattr(parameter, 'element_name') else str(parameter.name),
        )


def _slug(s):
    s = ''.join(c if c.isalnum() else '_' for c in str(s).lower())
    while '__' in s:
        s = s.replace('__', '_')
    return s.strip('_')[:60]


def latin_hypercube(n, d, rng):
    """Return an (n, d) Latin hypercube sample in (0, 1)."""
    if n == 0 or d == 0:
        return np.zeros((n, d))
    u = (rng.random((n, d)) + np.arange(n)[:, None]) / n
    for j in range(d):
        u[:, j] = u[rng.permutation(n), j]
    return np.clip(u, 1e-6, 1 - 1e-6)


# %% Adapter interface

class BiorefineryAdapter:
    """Base class for biorefinery adapters.

    Subclasses must implement every method that raises NotImplementedError.
    Attributes set in ``load`` (``system``, ``tea``, ``unit_groups``) are used
    by the generic collectors.
    """
    #: Short biorefinery label, e.g. 'HP'.
    name: str = None
    #: Maps UnitGroup names (as in the biorefinery) to CANONICAL_AREAS.
    area_map: dict = {}

    system = None
    tea = None
    unit_groups = ()
    #: Main product stream (set in ``load``); used to exclude its own sales
    #: from ``coproduct_sales``.
    product_stream = None

    # ---- configurations ----
    def configurations(self) -> list[dict]:
        """Return a list of configuration dicts. Each must contain the keys
        config_id, feedstock, fermentation_product, final_product,
        neutralization (bool), separation_variant (str), aerobic (bool),
        ferm_product_chemical (thermosteam ID), final_product_chemical
        (thermosteam ID)."""
        raise NotImplementedError

    def load(self, config: dict) -> None:
        """Load (or reload) the system for the given configuration and set
        ``self.system``, ``self.tea``, ``self.unit_groups``."""
        raise NotImplementedError

    # ---- domain ----
    def domain(self) -> dict:
        """Return dict with theoretical_yield, inflection_yield, y_min, y_max,
        t_min, t_max (g/g sugar and g/L) and baseline_productivity (g/L/h),
        baseline_yield, baseline_titer for the loaded configuration."""
        raise NotImplementedError

    def max_titer(self, y: float) -> Optional[float]:
        """Optional analytical upper bound on feasible titer at yield y."""
        return None

    # ---- parameters ----
    def parameters(self) -> list[Parameter]:
        """Return the uncertain parameters for the loaded configuration
        (excluding yield, titer and productivity)."""
        raise NotImplementedError

    # ---- simulation ----
    def set_try(self, y: float, t: float, r: float) -> None:
        """Set yield [g/g sugar], titer [g/L], productivity [g/L/h]. Raise
        InfeasibleSpecification if infeasible."""
        raise NotImplementedError

    def simulate(self) -> None:
        """Simulate the system robustly (with retries/solver switching).
        Raise InfeasibleSpecification for infeasible specifications."""
        raise NotImplementedError

    def reset(self) -> None:
        """Restore a converged state (e.g., baseline TRY) after a failure."""
        pass

    def collect(self) -> dict:
        """Return a dict with the adapter-specific POINT_COLUMNS (MPSP,
        recovery, F_* flows, V_broth, etc.). Generic TEA and area metrics are
        collected by the harness."""
        raise NotImplementedError

    # ---- descriptors ----
    def descriptors(self, config: dict) -> dict:
        """Return descriptor dict for the loaded configuration. The default
        implementation describes the fermentation and final product chemicals
        and the feedstock stream (``self.feedstock_stream``)."""
        out = {}
        chems = _get_chemicals(self)
        for key, prefix in (('ferm_product_chemical', 'fp_'),
                            ('final_product_chemical', 'pp_')):
            ID = config.get(key)
            if ID and chems is not None and ID in chems:
                out.update({prefix + k: v for k, v in describe_chemical(chems[ID]).items()})
        fs = getattr(self, 'feedstock_stream', None)
        if fs is not None:
            out.update(describe_feedstock(fs))
        return out

    def versions(self) -> dict:
        return _versions()


# %% Generic collectors

def _safe(f, default=float('nan')):
    try:
        v = f()
        if v is None:
            return default
        return float(v)
    except Exception:
        return default


def collect_tea_metrics(tea, product_stream=None) -> dict:
    """TEA metrics (annual basis). ``coproduct_sales`` excludes the sales of
    ``product_stream`` (the main product) at its current price."""
    d = {}
    for k, attr in (('AOC', 'AOC'), ('TCI', 'TCI'), ('FCI', 'FCI'),
                    ('installed_cost', 'installed_equipment_cost'),
                    ('material_cost', 'material_cost'),
                    ('utility_cost', 'utility_cost'),
                    ('coproduct_sales', 'sales'), ('FOC', 'FOC'), ('VOC', 'VOC'),
                    ('operating_hours', 'operating_hours'), ('IRR', 'IRR'),
                    ('income_tax', 'income_tax')):
        d[k] = _safe(lambda: getattr(tea, attr))
    if product_stream is not None:
        own = _safe(lambda: product_stream.price * product_stream.F_mass, 0.) * d['operating_hours']
        d['product_stream_sales'] = own
        d['coproduct_sales'] = d['coproduct_sales'] - own
    return d


def _utility_cost(units):
    cost = 0.
    for u in units:
        try:
            cost += u.utility_cost
        except Exception:
            pass
    return cost


def _group_material_cost(group, counted):
    """Material cost [USD/h] of the system feeds entering ``group``'s units.

    Unlike ``UnitGroup.get_material_cost`` (which, with the default
    ``extend_feed_ends=True``, traces feeds back through upstream storage/
    handling units and can count the same feed in several groups), each feed
    is attributed once, to the group whose units it enters first.
    """
    import biosteam as bst
    inlets = set(bst.utils.feeds_from_units(group.units))
    bst.utils.filter_out_missing_streams(inlets)
    cost = 0.
    for s in bst.utils.feeds(inlets):
        if id(s) in counted:
            continue
        counted.add(id(s))
        try:
            cost += s.cost
        except Exception:
            pass
    return cost


def _get_chemicals(adapter):
    for getter in (lambda: adapter.system.thermo.chemicals,
                   lambda: adapter.feedstock_stream.chemicals,
                   lambda: adapter.product_stream.chemicals):
        try:
            c = getter()
            if c is not None:
                return c
        except Exception:
            pass
    try:
        import thermosteam as tmo
        return tmo.settings.chemicals
    except Exception:
        return None


def collect_area_metrics(unit_groups: Iterable, area_map: dict) -> dict:
    """Sum UnitGroup metrics into canonical areas. Unmapped groups go to
    'facilities'. Material costs attribute each system feed to exactly one
    group (see :func:`_group_material_cost`)."""
    acc = {a: {q: 0. for q in AREA_QUANTITIES} for a in CANONICAL_AREAS}
    counted = set()
    for g in unit_groups:
        area = area_map.get(g.name, 'facilities')
        if area is None:  # explicitly ignored group
            continue
        a = acc[area]
        a['installed_cost'] += _safe(g.get_installed_cost, 0.) * 1e6  # MM$ -> $
        a['material_cost'] += _safe(lambda: _group_material_cost(g, counted), 0.)  # USD/h
        a['heating_duty'] += _safe(g.get_heating_duty, 0.)  # GJ/h
        a['cooling_duty'] += _safe(g.get_cooling_duty, 0.)  # GJ/h
        a['electricity_consumption'] += _safe(g.get_electricity_consumption, 0.)  # MW
        a['electricity_production'] += _safe(g.get_electricity_production, 0.)  # MW
        a['utility_cost'] += _utility_cost(g.units)
    return {f'{area}__{q}': v for area, qs in acc.items() for q, v in qs.items()}


def describe_chemical(chem, T_ferm=310.15) -> dict:
    """Thermophysical descriptors of a thermosteam Chemical (NaN if
    unavailable)."""
    d = {}
    d['MW'] = _safe(lambda: chem.MW)
    d['Tb'] = _safe(lambda: chem.Tb)
    d['Tm'] = _safe(lambda: chem.Tm)
    d['Tc'] = _safe(lambda: chem.Tc)
    d['Pc'] = _safe(lambda: chem.Pc)
    d['omega'] = _safe(lambda: chem.omega)
    d['Hvap_298'] = _safe(lambda: chem.Hvap(298.15))  # J/mol
    d['Psat_ferm'] = _safe(lambda: chem.Psat(T_ferm))  # Pa
    d['Hf'] = _safe(lambda: chem.Hf)  # J/mol
    d['LHV'] = _safe(lambda: chem.LHV)  # J/mol
    d['HHV'] = _safe(lambda: chem.HHV)  # J/mol
    d['Cn_l_298'] = _safe(lambda: chem.Cn('l', 298.15))
    d['rho_l_298'] = _safe(lambda: chem.rho('l', 298.15, 101325))
    atoms = {}
    try:
        atoms = dict(chem.atoms)
    except Exception:
        pass
    for el in ('C', 'H', 'O', 'N', 'S', 'Na', 'K', 'Ca'):
        d[f'n_{el}'] = float(atoms.get(el, 0))
    try:
        d['formula'] = chem.formula
    except Exception:
        d['formula'] = ''
    return d


_SUGAR_KEYS = {
    'glucose': ('Glucose', 'Dextrose', 'GlucoseMonohydrate'),
    'glucan': ('Glucan', 'Cellulose'),
    'sucrose': ('Sucrose',),
    'starch': ('Starch',),
    'xylose': ('Xylose',),
    'xylan': ('Xylan', 'Hemicellulose'),
    'arabinose': ('Arabinose',),
    'arabinan': ('Arabinan',),
    'galactose': ('Galactose',), 'galactan': ('Galactan',),
    'mannose': ('Mannose',), 'mannan': ('Mannan',),
    'fructose': ('Fructose',),
    'lignin': ('Lignin', 'SolubleLignin'),
    'ash': ('Ash',),
    'protein': ('Protein',),
    'extractives': ('Extract', 'Extractives'),
    'lipid': ('Lipid', 'Oil', 'TAG', 'OleicAcid'),
    'fiber': ('Fiber',),
    'water': ('Water', 'H2O'),
}


def describe_feedstock(stream) -> dict:
    """Composition descriptors of a feedstock stream (mass fractions, wet
    basis), price and flow."""
    d = {'fs_F_mass': _safe(lambda: stream.F_mass),
         'fs_price': _safe(lambda: stream.price)}
    F = max(_safe(lambda: stream.F_mass, 0.), 1e-12)
    try:
        IDs = stream.chemicals.IDs
        mass = stream.mass
    except Exception:
        return d
    used = set()
    for key, names in _SUGAR_KEYS.items():
        tot = 0.
        for i, ID in enumerate(IDs):
            if ID in names or any(ID.startswith(n) for n in names):
                if ID not in used:
                    tot += float(mass[i])
                    used.add(ID)
        d[f'fs_w_{key}'] = tot / F
    d['fs_w_other'] = max(0., 1. - sum(v for k, v in d.items() if k.startswith('fs_w_')))
    # glucose-equivalent fermentable sugar mass fraction (hydrolysis mass gain)
    gain = {'glucan': 180.156 / 162.141, 'starch': 180.156 / 162.141,
            'sucrose': 2 * 180.156 / 342.297, 'xylan': 150.13 / 132.115,
            'arabinan': 150.13 / 132.115, 'galactan': 180.156 / 162.141,
            'mannan': 180.156 / 162.141}
    d['fs_w_fermentable_sugar_eq'] = sum(
        d.get(f'fs_w_{k}', 0.) * gain.get(k, 1.)
        for k in ('glucose', 'glucan', 'sucrose', 'starch', 'xylose', 'xylan',
                  'arabinose', 'arabinan', 'galactose', 'galactan', 'mannose',
                  'mannan', 'fructose'))
    return d


def _versions():
    v = {'python': platform.python_version()}
    for mod in ('biosteam', 'thermosteam', 'numpy', 'scipy'):
        try:
            v[mod] = __import__(mod).__version__
        except Exception:
            v[mod] = None
    try:
        import biorefineries
        path = os.path.dirname(os.path.dirname(biorefineries.__file__))
        v['bioindustrial_park_commit'] = subprocess.check_output(
            ['git', '-C', path, 'rev-parse', 'HEAD'], text=True).strip()
        v['bioindustrial_park_branch'] = subprocess.check_output(
            ['git', '-C', path, 'rev-parse', '--abbrev-ref', 'HEAD'], text=True).strip()
    except Exception:
        pass
    return v


# %% TRY design

def try_design(dom, n_y, n_t, rng, jitter, max_titer):
    """Return list of (y, t_upper) yields and the per-yield titer grids."""
    y_lo, y_hi = dom['y_min'], dom['y_max']
    ys = np.linspace(y_lo, y_hi, n_y)
    if jitter and n_y > 2:
        step = (y_hi - y_lo) / (n_y - 1)
        ys[1:-1] += rng.uniform(-0.3, 0.3, n_y - 2) * step
    design = []
    for y in ys:
        t_hi = dom['t_max']
        mt = max_titer(y)
        if mt is not None and np.isfinite(mt):
            t_hi = min(t_hi, mt)
        t_hi = max(t_hi, dom['t_min'] * 1.05)
        ts = np.linspace(dom['t_min'], t_hi, n_t)
        if jitter and n_t > 2:
            step = (t_hi - dom['t_min']) / (n_t - 1)
            ts[1:-1] += rng.uniform(-0.3, 0.3, n_t - 2) * step
        design.append((float(y), [float(t) for t in ts]))
    return design


# %% Campaign

class _Alarm:
    def __init__(self, seconds):
        self.seconds = int(seconds)

    def __enter__(self):
        if self.seconds > 0:
            signal.signal(signal.SIGALRM, self._handler)
            signal.alarm(self.seconds)

    def _handler(self, signum, frame):
        raise PointTimeout(f'point exceeded {self.seconds} s')

    def __exit__(self, *args):
        if self.seconds > 0:
            signal.alarm(0)
        return False


def _jsonable(v):
    if isinstance(v, (np.floating, float)):
        v = float(v)
        return v if math.isfinite(v) else None
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.bool_,)):
        return bool(v)
    return v


def _write_row(fh, row):
    fh.write(json.dumps({k: _jsonable(v) for k, v in row.items()}) + '\n')
    fh.flush()


def _completed_samples(path):
    done = set()
    if os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                try:
                    rec = json.loads(line)
                except Exception:
                    continue
                if rec.get('record') == 'sample_end':
                    done.add(rec['sample_id'])
    return done


def _sample_plan(adapter, config, n_random, n_productivity, seed, prod_range):
    """Deterministic list of sample specs for a configuration."""
    params = adapter.parameters()
    rng = np.random.default_rng([seed, _stable_hash(config['config_id'])])
    plan = [{'kind': 'baseline', 'u': None, 'r_mult': 1.0}]
    lo, hi = np.log(prod_range[0]), np.log(prod_range[1])
    for k in range(n_productivity):
        # stratified over log-productivity, baseline parameters
        q = (k + rng.uniform(0.1, 0.9)) / max(n_productivity, 1)
        plan.append({'kind': 'productivity', 'u': None,
                     'r_mult': float(np.exp(lo + (hi - lo) * q))})
    U = latin_hypercube(n_random, len(params) + 1, rng)
    for k in range(n_random):
        plan.append({'kind': 'random', 'u': U[k, :-1],
                     'r_mult': float(np.exp(lo + (hi - lo) * U[k, -1]))})
    for i, p in enumerate(plan):
        p['index'] = i
        p['seed'] = int(rng.integers(0, 2**31 - 1))
    return params, plan


def _stable_hash(s):
    h = 0
    for c in str(s):
        h = (h * 131 + ord(c)) % 2**31
    return h


def run_campaign(adapter: BiorefineryAdapter, out_dir: str, *,
                 config_ids: Optional[list] = None,
                 n_random: int = 20, n_productivity: int = 4, seed: int = 1,
                 n_y: int = 5, n_t: int = 5,
                 productivity_range=(0.2, 5.0),
                 point_timeout: float = 600., sample_limit: Optional[int] = None,
                 tag: str = '', verbose: bool = True):
    """Run (or resume) a data-generation campaign.

    Writes ``{out_dir}/{adapter.name}{tag}.jsonl`` containing records:
    'config' (descriptors), 'sample_start', 'point', 'sample_end'.
    """
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f'{adapter.name}{tag}.jsonl')
    done = _completed_samples(path)
    configs = adapter.configurations()
    if config_ids:
        configs = [c for c in configs if c['config_id'] in config_ids]
    versions = adapter.versions()
    n_run = 0
    with open(path, 'a') as fh:
        for config in configs:
            if sample_limit is not None and n_run >= sample_limit:
                return path
            cid = config['config_id']
            if verbose:
                print(f'\n=== {adapter.name}: loading {cid}', flush=True)
            t0 = time.time()
            adapter.load(config)
            dom = adapter.domain()
            params, plan = _sample_plan(adapter, config, n_random, n_productivity,
                                        seed, productivity_range)
            pending = [p for p in plan if f'{adapter.name}|{cid}|s{seed}|{p["index"]}' not in done]
            if not pending:
                continue
            desc = adapter.descriptors(config)
            _write_row(fh, {'record': 'config', 'biorefinery': adapter.name,
                            'config_id': cid, **config, **dom, **desc,
                            'parameters': json.dumps([(p.name, p.baseline, p.description) for p in params]),
                            'versions': json.dumps(versions),
                            'load_seconds': time.time() - t0})
            for p in pending:
                if sample_limit is not None and n_run >= sample_limit:
                    return path
                sid = f'{adapter.name}|{cid}|s{seed}|{p["index"]}'
                _run_sample(adapter, fh, sid, config, dom, params, p,
                            n_y, n_t, point_timeout, verbose)
                n_run += 1
    return path


def _run_sample(adapter, fh, sid, config, dom, params, plan_item,
                n_y, n_t, point_timeout, verbose):
    rng = np.random.default_rng(plan_item['seed'])
    # parameter values
    values = {}
    for i, prm in enumerate(params):
        if plan_item['u'] is None:
            values[prm.name] = prm.baseline
        else:
            values[prm.name] = prm.inverse_cdf(float(plan_item['u'][i]))
    r = dom['baseline_productivity'] * plan_item['r_mult']
    base = {'record': 'point', 'biorefinery': adapter.name, 'sample_id': sid,
            'sample_kind': plan_item['kind'], 'config_id': config['config_id'],
            'r_mult': plan_item['r_mult']}
    t_sample = time.time()
    _write_row(fh, {'record': 'sample_start', 'sample_id': sid, **values,
                    'r': r, 'sample_kind': plan_item['kind']})
    param_error = ''
    try:
        for prm in params:
            prm.setter(values[prm.name])
    except Exception as e:  # parameter setting failed: record and skip
        param_error = f'{type(e).__name__}: {e}'[:300]
    n_ok = 0
    design = try_design(dom, n_y, n_t, rng, plan_item['kind'] != 'baseline',
                        adapter.max_titer) if not param_error else []
    point_id = 0
    for y, ts in design:
        for t in ts:
            row = dict(base, point_id=point_id, y=y, t=t, r=r)
            point_id += 1
            t0 = time.time()
            try:
                with _Alarm(point_timeout):
                    adapter.set_try(y, t, r)
                    adapter.simulate()
                    res = adapter.collect()
                    res.update(collect_tea_metrics(adapter.tea,
                                                   getattr(adapter, 'product_stream', None)))
                    res.update(collect_area_metrics(adapter.unit_groups, adapter.area_map))
                mpsp = res.get('MPSP', float('nan'))
                if mpsp is None or not np.isfinite(mpsp) or mpsp <= 0:
                    raise RuntimeError(f'invalid MPSP {mpsp}')
                row.update(res)
                row['status'] = 'ok'
                n_ok += 1
            except InfeasibleSpecification as e:
                row['status'] = 'infeasible'
                row['error'] = str(e)[:300]
            except Exception as e:
                row['status'] = 'failed'
                row['error'] = f'{type(e).__name__}: {e}'[:300]
                if verbose:
                    traceback.print_exc(limit=2)
                try:
                    with _Alarm(point_timeout):
                        adapter.reset()
                        for prm in params:
                            prm.setter(values[prm.name])
                except Exception:
                    pass
            row['sim_seconds'] = time.time() - t0
            _write_row(fh, row)
    _write_row(fh, {'record': 'sample_end', 'sample_id': sid, 'n_ok': n_ok,
                    'param_error': param_error,
                    'seconds': time.time() - t_sample})
    if verbose:
        print(f'{sid}: {n_ok} ok points in {time.time() - t_sample:.0f} s', flush=True)
    # restore baseline parameters for the next sample
    try:
        for prm in params:
            prm.setter(prm.baseline)
    except Exception:
        pass
