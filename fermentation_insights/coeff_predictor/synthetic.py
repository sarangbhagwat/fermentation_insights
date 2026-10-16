# -*- coding: utf-8 -*-
"""
Synthetic raw-data generator with known ground truth.

Fake biorefineries are implemented as :class:`BiorefineryAdapter` subclasses
with fake TEA / UnitGroup objects and are run through the *real*
``data_generation.harness.run_campaign``, so the JSON-lines output has exactly
the harness format (config / sample_start / point / sample_end records,
TEA metrics via ``collect_tea_metrics``, area metrics via
``collect_area_metrics``).

Ground-truth structure of a sample (parameters constant within a sample)::

    F_sugar   = F_feedstock * w_sugar * eta                [kg/h]
    F_ferm    = y * F_sugar;   F_final = g(t) * F_ferm      [kg/h]
    V_broth   = F_ferm / t  [m3/h];   V_fermenters = 1.2 F_ferm / r  [m3]
    installed costs (USD) linear in F_feedstock, V_fermenters (~1/r),
        V_broth (~1/t), F_ferm, unconverted sugar, duties
    electricity, heating/cooling duties, natural gas (purchased only when the
        heating demand exceeds the residue fuel supply -> non-linear),
        base/acid, enzymes/nutrients, other utilities, coproduct sales
    FCI = 1.6 * (I + 0.175 * I_ISBL)   (site development/piping/warehouse
          factors on ISBL areas only -> FCI/I varies across TY points)
    R = sum(opex) + phi*FCI + F (+ optional pseudo-noise of relative size kappa)
    MPSP = R / (F_final * h)

``phi`` (capital charge) is a function of IRR/tax, ``F`` the fixed cost.
:func:`ground_truth` recomputes the true channels for every recorded point and
fits them with :func:`physics.fit_ols_many` (primed basis), so a correct
decomposition of the JSONL columns must recover them (exactly for
``kappa = 0``).

Harness v2 objects: the fake TEA carries ``labor_cost`` (and
``labor_burden``, as ``CellulosicEthanolTEA``; SYNC mimics
``ConventionalEthanolTEA``: no burden attribute, multiplier 1.6, no OSBL
units) and ``OSBL_units`` (wastewater, CHP, facilities); natural gas is an
inlet stream registered as a 'Natural gas' utility of the boiler -- unpriced
(stream utility, charged at the stand-in ``biosteam.stream_prices``; SYNA,
SYNC) or priced (material pool; SYNB) -- and SYNA also feeds unpriced methane
to its dryer (material pool, cost 0, like TAL); ``feedstock_stream`` gives the
per-point ``fs_*`` descriptors (price = the sampled feedstock price).
Generation runs with a minimal stand-in ``biosteam`` module
(:func:`biosteam_stub`), so results do not depend on whether (or which)
BioSTEAM is installed. ``generate(..., harness_version=1)`` strips the v2
columns to emulate legacy files.

CLI::

    python -m fermentation_insights.coeff_predictor.synthetic --out <dir> \
        [--n-random 6] [--n-productivity 2] [--kappa 0]
"""
from __future__ import annotations

import argparse
import contextlib
import json
import os
import sys
import types

import numpy as np
import pandas as pd

from .data_generation.harness import (BiorefineryAdapter, InfeasibleSpecification,
                                      Parameter, run_campaign)
from . import physics

__all__ = ('SYNTHETIC_BIOREFINERIES', 'SyntheticAdapter', 'point_model',
           'generate', 'ground_truth', 'biosteam_stub', 'strip_v2', 'V2_POINT_KEYS')

#: natural-gas utility name of the fake boiler inlet
NG_NAME = 'Natural gas'
#: per-point keys added by harness v2 (``strip_v2`` removes them, plus fs_*)
V2_POINT_KEYS = physics.NG_POOL_COLUMNS + ('labor_cost', 'labor_burden')


def _stub_module():
    m = types.ModuleType('biosteam')
    m.__version__ = 'synthetic-stub'
    m.stream_prices = {}
    m.utils = types.SimpleNamespace(
        feeds_from_units=lambda units: [s for u in units for s in getattr(u, 'ins', ())],
        filter_out_missing_streams=lambda streams: None,
        feeds=lambda streams: [s for s in streams if getattr(s, 'source', None) is None])
    return m


_BST = _stub_module()


@contextlib.contextmanager
def biosteam_stub():
    """Temporarily make ``import biosteam`` return a minimal stand-in that
    provides what the harness collectors use (``stream_prices``,
    ``utils.feeds_from_units/filter_out_missing_streams/feeds``)."""
    old = sys.modules.get('biosteam')
    sys.modules['biosteam'] = _BST
    try:
        yield _BST
    finally:
        if old is None:
            sys.modules.pop('biosteam', None)
        else:
            sys.modules['biosteam'] = old


class _Stream:
    """Hashable fake stream (feed: ``source`` None)."""
    source = None

    def __init__(self, F_mass=0., price=0., cost=None, CH4=0., IDs=(), mass=()):
        self.F_mass, self.price = float(F_mass), float(price)
        self._cost = cost
        self.imass = {'CH4': float(CH4)}
        self.chemicals = _Chemicals(tuple(IDs) + (('CH4',) if CH4 else ()))
        self.mass = np.asarray(mass, float)

    @property
    def cost(self):
        return self.F_mass * self.price if self._cost is None else self._cost


class _Chemicals(tuple):
    @property
    def IDs(self):
        return tuple(self)


class _Unit:
    """Hashable fake unit."""
    def __init__(self):
        self.utility_cost = 0.
        self.installed_cost = 0.
        self.ins = []
        self._inlet_utility_indices = {}

# %% Definitions

_BASE_CONFIG = dict(n_acid_groups=0, pKa1=float('nan'), upgrading=False,
                    separation_variant='', aerobic=False,
                    fermentation_mode='batch', T_fermentation=305.15,
                    product_water_solubility=1000., neutralizing_agent='')

SYNTHETIC_BIOREFINERIES = {
    'SYNA': dict(
        feedstock='glucose', feedstock_type='pure_sugar', w_sugar=0.95, eta=0.98,
        w_lignin=0.0, capacity=20000., fs_price=0.40, ng_pool='utility',
        labor_burden=0.9, osbl=True, ng_unpriced=0.05,
        fp='SynAcid', pp='SynAcid', fp_MW=90.08, pp_MW=90.08, stoich=1.0,
        theoretical_yield=1.0, y_min=0.3, inflection_yield=0.9, t_min=20., t_max=150.,
        S_max=600., baseline_productivity=0.8, g0=0.85, t_g=0.0, steam_supply=0.0,
        configs=[dict(config_id='SYNA_neutral', neutralization=True, acidulation=True,
                      neutralizing_agent='Ca(OH)2', n_acid_groups=1, pKa1=3.86,
                      separation_principle='acidulation+evaporation', c_sep=2.0e5),
                 dict(config_id='SYNA_lowpH', neutralization=False, acidulation=False,
                      n_acid_groups=1, pKa1=3.86,
                      separation_principle='evaporation', c_sep=2.6e5)],
        drop_columns=(),
    ),
    'SYNB': dict(
        feedstock='cornstover', feedstock_type='lignocellulose', w_sugar=0.60, eta=0.85,
        w_lignin=0.18, capacity=100000., fs_price=0.09, ng_pool='material',
        labor_burden=0.9, osbl=True, ng_unpriced=0.,
        fp='SynAlcohol', pp='SynFuel', fp_MW=74.12, pp_MW=56.1, stoich=0.757,
        theoretical_yield=0.41, y_min=0.12, inflection_yield=0.38, t_min=10., t_max=120.,
        S_max=300., baseline_productivity=0.5, g0=0.70, t_g=1.0, steam_supply=1.2e-3,
        configs=[dict(config_id='SYNB_distill', neutralization=False, acidulation=False,
                      upgrading=True, separation_principle='distillation', c_sep=1.5e5),
                 dict(config_id='SYNB_extract', neutralization=False, acidulation=False,
                      upgrading=True, separation_variant='extraction',
                      separation_principle='extraction+distillation', c_sep=1.1e5)],
        drop_columns=(),
    ),
    'SYNC': dict(
        feedstock='sugarcane', feedstock_type='sugar_juice', w_sugar=0.14, eta=0.95,
        w_lignin=0.05, capacity=300000., fs_price=0.035, ng_pool='utility',
        labor_burden=None, labor_multiplier=1.6, osbl=False, ng_unpriced=0.,
        fp='SynDiol', pp='SynDiol', fp_MW=90.12, pp_MW=90.12, stoich=1.0,
        theoretical_yield=0.5, y_min=0.15, inflection_yield=0.45, t_min=15., t_max=130.,
        S_max=450., baseline_productivity=1.2, g0=0.80, t_g=0.5, steam_supply=0.15e-3,
        configs=[dict(config_id='SYNC_aerobic', neutralization=False, acidulation=False,
                      aerobic=True, separation_principle='evaporation+distillation',
                      c_sep=1.8e5)],
        # columns this fake adapter does not provide (loader must cope)
        drop_columns=('acid_cost', 'GWP', 'FEC', 'tau_fermentation'),
    ),
}

_PARAMS = (
    # name, baseline multiplier/value, (lo, hi) relative to baseline or absolute
    ('p_feedstock_price', 'fs_price', (0.75, 1.25), 'rel'),
    ('p_electricity_price', 0.07, (0.03, 0.11), 'abs'),
    ('p_natural_gas_price', 0.25, (0.15, 0.40), 'abs'),
    ('p_IRR', 0.10, (0.05, 0.15), 'abs'),
    ('p_income_tax', 0.21, (0.15, 0.35), 'abs'),
    ('p_operating_days', 330., (300., 350.), 'abs'),
    ('p_base_price', 0.20, (0.15, 0.30), 'abs'),
    ('p_acid_price', 0.10, (0.07, 0.14), 'abs'),
    ('bp_capex_scale', 1.0, (0.75, 1.3), 'abs'),
    ('bp_labor_scale', 1.0, (0.7, 1.4), 'abs'),
)


def capital_charge(IRR, tax, life=20.):
    """Synthetic effective capital charge per USD of FCI [1/yr] (CRF x
    TCI/FCI x tax gross-up + FCI-proportional fixed costs)."""
    crf = IRR / (1 - (1 + IRR) ** -life)
    return float(crf * 1.05 / (1 - 0.5 * tax) + 0.03)


def fixed_cost(bio, params):
    return float(2.0e6 * (bio['capacity'] / 1e5) ** 0.25 * params['bp_labor_scale'])


def point_model(bio, cfg, params, y, t, r, kappa=0.):
    """All quantities at one TY point (dict). Raises InfeasibleSpecification
    above the sugar-concentration limit."""
    if t > y * bio['S_max'] * (1 + 1e-9):
        raise InfeasibleSpecification(f'titer {t:.1f} needs sugar > {bio["S_max"]} g/L')
    h = params['p_operating_days'] * 24.
    cs = params['bp_capex_scale']
    lig = bio['feedstock_type'] == 'lignocellulose'
    F_fs = bio['capacity']
    S = F_fs * bio['w_sugar'] * bio['eta']
    F_ferm = y * S
    g = bio['g0'] * t / (t + bio['t_g'])
    F_final = g * F_ferm
    V_broth = F_ferm / t
    V_ferm = 1.2 * F_ferm / r
    water = V_broth * 1000. * 0.95
    # duties GJ/h
    HD = 2.5e-3 * 0.25 * water + 0.4e-3 * F_fs + 1.0e-3 * F_ferm
    CD = 0.6 * HD + (8e-3 if cfg.get('aerobic') else 1e-3) * F_ferm
    supply = bio['steam_supply'] * F_fs * (1. + 0.5 * (1 - y / bio['theoretical_yield']))
    ng_flow = max(0., HD - supply) / (0.050 * 0.85)  # kg/h (LHV 50 MJ/kg, 85 %)
    # electricity MW
    e_cons = 1e-3 * (0.02 * F_fs + 0.6 * V_broth + 0.004 * V_ferm
                     + (0.5 if cfg.get('aerobic') else 0.05) * F_ferm + 0.01 * CD * 1e3)
    e_prod = 1e-3 * (0.9 * bio['w_lignin'] * F_fs) if bio['w_lignin'] > 0 else 0.
    # materials USD/h
    p_fs = params['p_feedstock_price']
    fs_cost = F_fs * p_fs
    base_cost = params['p_base_price'] * 0.45 * F_ferm if cfg.get('neutralization') else 0.
    acid_cost = params['p_acid_price'] * 0.6 * F_ferm if cfg.get('acidulation') else 0.
    other_mat = (0.012 * F_fs if lig else 0.) + 0.004 * S + (0.05 * F_ferm if cfg.get('upgrading') else 0.)
    ng_cost = ng_flow * params['p_natural_gas_price']
    el_cost = (e_cons - e_prod) * 1000. * params['p_electricity_price']
    other_util = 0.45 * CD + 0.02 * HD + 0.002 * F_fs   # cooling water, chilled water, ash
    cop_sales = (0.02 * 1.7 * F_ferm if cfg.get('acidulation') else 0.) \
        + (0.003 * F_fs if bio['feedstock_type'] == 'sugar_juice' else 0.)  # USD/h
    # installed costs USD per area
    I_area = {
        'feedstock': 150. * F_fs,
        'conversion': (600. if lig else 60. if bio['feedstock_type'] == 'sugar_juice' else 0.) * F_fs,
        'fermentation': 3000. * V_ferm + 2.0e4 * V_broth,
        'separation': cfg['c_sep'] * V_broth + 600. * F_ferm,
        'upgrading': 1500. * F_ferm if cfg.get('upgrading') else 0.,
        'wastewater': 4.0e4 * V_broth + 900. * (S - F_ferm),
        'chp': 400. * F_fs if lig else 3.0e5 * HD,
        'facilities': 1.0e5 * CD + 5.0e6,
    }
    I_area = {k: v * cs for k, v in I_area.items()}
    I = sum(I_area.values())
    I_ISBL = sum(v for k, v in I_area.items() if k not in ('wastewater', 'chp', 'facilities'))
    FCI = 1.6 * (I + 0.175 * I_ISBL)
    phi = capital_charge(params['p_IRR'], params['p_income_tax'])
    F = fixed_cost(bio, params)
    E = dict(feedstock=fs_cost * h, electricity=el_cost * h, natural_gas=ng_cost * h,
             neutralization=(base_cost + acid_cost) * h, materials_other=other_mat * h,
             utilities_other=other_util * h, coproducts=-cop_sales * h)
    R_lin = sum(E.values()) + phi * FCI + F
    # optional deterministic pseudo-noise (numerical/solver tolerance-like,
    # relative amplitude kappa); not representable by the channel model
    eps = kappa * R_lin * np.sin(1299.0 * y + 78.233 * t + 3.7 * r)
    R = R_lin + eps
    MPSP = R / (F_final * h)
    ng_in_mat = bio['ng_pool'] == 'material'
    return dict(
        h=h, MPSP=MPSP, recovery=g, F_final_product=F_final, F_ferm_product=F_ferm,
        F_sugar_fed=S, F_feedstock=F_fs, feedstock_price=p_fs, V_broth=V_broth,
        F_water_broth=water, V_fermenters=V_ferm, tau_fermentation=t / r,
        base_cost=base_cost, acid_cost=acid_cost, heating_duty=HD, cooling_duty=CD,
        electricity_consumption=e_cons, electricity_production=e_prod,
        natural_gas_flow=ng_flow, electricity_price=params['p_electricity_price'],
        natural_gas_price=params['p_natural_gas_price'],
        material_cost_h=fs_cost + base_cost + acid_cost + other_mat + (ng_cost if ng_in_mat else 0.),
        utility_cost_h=el_cost + other_util + (0. if ng_in_mat else ng_cost),
        cop_sales_h=cop_sales, I_area=I_area, I=I, FCI=FCI, phi=phi, F=F, E=E, eps=eps,
        el_cost=el_cost, other_util=other_util, ng_cost=ng_cost, fs_cost=fs_cost,
        other_mat=other_mat, ng_in_mat=ng_in_mat, I_ISBL=I_ISBL,
        # unpriced methane to a dryer (material pool, cost 0; not a cost channel)
        ng_unpriced_flow=bio.get('ng_unpriced', 0.) * F_ferm,
    )


#: OSBL areas of the synthetic TEAs (``bio['osbl']``)
OSBL_AREAS = ('wastewater', 'chp', 'facilities')


def area_fci_weights(res, osbl=True):
    """True per-area capital weights: FCI contributions ``I_a + 0.175*ISBL_a``
    (synthetic FCI = 1.6*(I + 0.175*ISBL)) when the OSBL split is recorded,
    else installed costs (what the pipeline can know)."""
    return {a: v * (1. if a in OSBL_AREAS else 1.175) if osbl else v
            for a, v in res['I_area'].items()}


# %% Fake biosteam-like objects

class _Group:
    def __init__(self, name):
        self.name = name
        self.I = self.mat = self.hd = self.cd = self.ec = self.ep = 0.
        self.unit = _Unit()
        self.units = [self.unit]
        self.feed = _Stream(cost=0.)     # lumped purchased feeds of the group
        self.unit.ins = [self.feed]

    def sync(self):
        """Push I / mat into the fake unit (installed cost, feed cost net of
        the separately modelled feeds such as natural gas)."""
        self.unit.installed_cost = self.I
        self.feed._cost = self.mat - sum(s.cost for s in self.unit.ins if s is not self.feed)

    def get_installed_cost(self): return self.I / 1e6   # MM$
    def get_material_cost(self): return self.mat        # USD/h
    def get_heating_duty(self): return self.hd
    def get_cooling_duty(self): return self.cd
    def get_electricity_consumption(self): return self.ec
    def get_electricity_production(self): return self.ep


#: fake UnitGroup names -> canonical areas (one ignored group, one unmapped)
_GROUPS = {'feedstock_handling': 'feedstock', 'pretreatment': 'conversion',
           'fermentation': 'fermentation', 'separation': 'separation',
           'upgrading': 'upgrading', 'wastewater': 'wastewater',
           'boiler_turbogenerator': 'chp', 'cooling_tower': 'facilities',
           'duplicate_view': None}
_UNMAPPED = 'storage'   # not in area_map -> harness puts it in 'facilities'


class SyntheticAdapter(BiorefineryAdapter):
    """Fake adapter for one synthetic biorefinery."""

    def __init__(self, name, kappa=0., fail_every=0, drop_columns=None):
        self.name = name
        self.bio = SYNTHETIC_BIOREFINERIES[name]
        self.kappa = kappa
        self.fail_every = fail_every
        self.drop_columns = tuple(self.bio['drop_columns'] if drop_columns is None else drop_columns)
        self.area_map = {k: v for k, v in _GROUPS.items()}
        self.values = {}
        self._n = 0

    def configurations(self):
        out = []
        for c in self.bio['configs']:
            d = dict(_BASE_CONFIG)
            d.update(feedstock=self.bio['feedstock'], feedstock_type=self.bio['feedstock_type'],
                     fermentation_product=self.bio['fp'], final_product=self.bio['pp'],
                     ferm_product_chemical=self.bio['fp'], final_product_chemical=self.bio['pp'],
                     final_product_stoich_mass_ratio=self.bio['stoich'])
            d.update({k: v for k, v in c.items() if k not in ('acidulation', 'c_sep')})
            out.append(d)
        return out

    def load(self, config):
        self.config = config
        self.cfg = next(dict(c) for c in self.bio['configs'] if c['config_id'] == config['config_id'])
        self.cfg.setdefault('aerobic', False)
        self.unit_groups = [_Group(n) for n in list(_GROUPS) + [_UNMAPPED]]
        G = {g.name: g for g in self.unit_groups}
        self.tea = types.SimpleNamespace()
        if self.bio.get('osbl'):
            self.tea.OSBL_units = [G[n].unit for n, a in list(_GROUPS.items()) + [(_UNMAPPED, 'facilities')]
                                   if a in OSBL_AREAS]
        self.product_stream = types.SimpleNamespace(price=0., F_mass=0.)
        # natural gas: inlet 1 of the boiler, registered as a utility
        self.ng_stream = _Stream(CH4=1.)
        boiler = G['boiler_turbogenerator'].unit
        boiler.ins.append(self.ng_stream)
        boiler._inlet_utility_indices = {NG_NAME: 1}
        # unpriced methane to the dryer (separation), not registered: material pool
        self.dryer_gas = _Stream(CH4=1.)
        if self.bio.get('ng_unpriced'):
            G['separation'].unit.ins.append(self.dryer_gas)
        units = [g.unit for g in self.unit_groups]
        self.system = types.SimpleNamespace(
            cost_units=units, units=units, feeds=[s for u in units for s in u.ins])
        self.values = {p.name: p.baseline for p in self.parameters()}

    @property
    def feedstock_stream(self):
        """Fake feedstock stream (composition of the synthetic biorefinery,
        price = the sampled feedstock price)."""
        b = self.bio
        F = b['capacity']
        w = (b['w_sugar'], b['w_lignin'], 1. - b['w_sugar'] - b['w_lignin'])
        return _Stream(F_mass=F, price=self.values.get('p_feedstock_price', b['fs_price']),
                       IDs=('Glucose', 'Lignin', 'Water'), mass=[F * x for x in w])

    def domain(self):
        b = self.bio
        return dict(theoretical_yield=b['theoretical_yield'], inflection_yield=b['inflection_yield'],
                    y_min=b['y_min'], y_max=b['inflection_yield'], t_min=b['t_min'],
                    t_max=b['t_max'], baseline_productivity=b['baseline_productivity'],
                    baseline_yield=0.5 * (b['y_min'] + b['inflection_yield']),
                    baseline_titer=0.5 * (b['t_min'] + b['t_max']))

    def max_titer(self, y):
        return None  # deliberately None so some points are infeasible

    def parameters(self):
        out = []
        for name, base, (lo, hi), kind in _PARAMS:
            b = self.bio[base] if isinstance(base, str) else base
            if kind == 'rel':
                lo, hi = lo * b, hi * b
            out.append(Parameter(name=name, baseline=float(b),
                                 inverse_cdf=lambda u, lo=lo, hi=hi: float(lo + (hi - lo) * u),
                                 setter=lambda v, n=name: self.values.__setitem__(n, float(v)),
                                 description=f'synthetic {name}'))
        return out

    def set_try(self, y, t, r):
        self.try_ = (y, t, r)
        if t > y * self.bio['S_max'] * (1 + 1e-9):
            raise InfeasibleSpecification('sugar concentration limit')

    def simulate(self):
        self._n += 1
        if self.fail_every and self._n % self.fail_every == 0:
            raise RuntimeError('synthetic convergence failure')
        y, t, r = self.try_
        self.res = res = point_model(self.bio, self.cfg, self.values, y, t, r, self.kappa)
        h = res['h']
        purity = 0.98
        self.product_stream.price = res['MPSP'] * purity
        self.product_stream.F_mass = res['F_final_product'] / purity
        own = self.product_stream.price * self.product_stream.F_mass * h
        I = res['I']
        tea = self.tea
        tea.installed_equipment_cost = I
        tea.FCI = res['FCI']
        tea.TCI = 1.05 * res['FCI']
        tea.material_cost = res['material_cost_h'] * h
        tea.utility_cost = res['utility_cost_h'] * h
        tea.sales = res['cop_sales_h'] * h + own
        tea.FOC = res['F'] + 0.03 * tea.FCI
        tea.VOC = tea.material_cost + tea.utility_cost - res['cop_sales_h'] * h
        tea.AOC = tea.FOC + tea.VOC
        tea.operating_hours = h
        tea.IRR = self.values['p_IRR']
        tea.income_tax = self.values['p_income_tax']
        # unit groups
        G = {g.name: g for g in self.unit_groups}
        for g in self.unit_groups:
            g.I = g.mat = g.hd = g.cd = g.ec = g.ep = 0.
            g.units[0].utility_cost = 0.
        Ia = res['I_area']
        for gname, area in _GROUPS.items():
            if area is not None:
                G[gname].I = Ia[area]
        G['cooling_tower'].I -= 1.0e6 * self.values['bp_capex_scale']  # part in 'storage'
        G[_UNMAPPED].I = 1.0e6 * self.values['bp_capex_scale']
        G['duplicate_view'].I = 0.3 * Ia['separation']                 # ignored (None)
        G['feedstock_handling'].mat = res['fs_cost']
        G['fermentation'].mat = res['base_cost'] + 0.004 * res['F_sugar_fed']
        G['separation'].mat = res['acid_cost']
        G['pretreatment'].mat = res['other_mat'] - 0.004 * res['F_sugar_fed'] \
            - (0.05 * res['F_ferm_product'] if self.cfg.get('upgrading') else 0.)
        G['upgrading'].mat = 0.05 * res['F_ferm_product'] if self.cfg.get('upgrading') else 0.
        if res['ng_in_mat']:
            G['boiler_turbogenerator'].mat = res['ng_cost']
        G['fermentation'].hd = 0.2 * res['heating_duty']
        G['separation'].hd = 0.8 * res['heating_duty']
        G['fermentation'].cd = 0.5 * res['cooling_duty']
        G['separation'].cd = 0.5 * res['cooling_duty']
        G['fermentation'].ec = 0.6 * res['electricity_consumption']
        G['separation'].ec = 0.4 * res['electricity_consumption']
        G['boiler_turbogenerator'].ep = res['electricity_production']
        p_el = self.values['p_electricity_price']
        G['fermentation'].units[0].utility_cost = 0.6 * res['electricity_consumption'] * 1e3 * p_el
        G['separation'].units[0].utility_cost = 0.4 * res['electricity_consumption'] * 1e3 * p_el
        G['boiler_turbogenerator'].units[0].utility_cost = (
            -res['electricity_production'] * 1e3 * p_el + (0. if res['ng_in_mat'] else res['ng_cost']))
        G['cooling_tower'].units[0].utility_cost = res['other_util']
        # harness v2 objects: labor (FOC = labor*(1 + burden) + 0.03*FCI)
        b = self.bio
        if b.get('labor_burden') is not None:
            tea.labor_cost = res['F'] / (1. + b['labor_burden'])
            tea.labor_burden = b['labor_burden']
        else:   # ConventionalEthanolTEA-like: no labor_burden attribute
            tea.labor_cost = res['F'] / b.get('labor_multiplier', 1.)
        # natural gas: stream utility (unpriced stream charged at
        # biosteam.stream_prices) or priced feed (material pool)
        p_ng = self.values['p_natural_gas_price']
        self.ng_stream.F_mass = res['natural_gas_flow']
        self.ng_stream.imass['CH4'] = res['natural_gas_flow']
        self.ng_stream.price = p_ng if res['ng_in_mat'] else 0.
        _BST.stream_prices[NG_NAME] = p_ng
        self.dryer_gas.F_mass = self.dryer_gas.imass['CH4'] = res['ng_unpriced_flow']
        for g in self.unit_groups:
            g.sync()

    def collect(self):
        keys = ('MPSP', 'recovery', 'F_final_product', 'F_ferm_product', 'F_sugar_fed',
                'F_feedstock', 'feedstock_price', 'V_broth', 'F_water_broth', 'V_fermenters',
                'tau_fermentation', 'base_cost', 'acid_cost', 'heating_duty', 'cooling_duty',
                'electricity_consumption', 'electricity_production', 'natural_gas_flow',
                'electricity_price', 'natural_gas_price')
        out = {k: self.res[k] for k in keys}
        out['GWP'] = float('nan')
        out['FEC'] = float('nan')
        out['x_true_eps'] = self.res['eps']
        for k in self.drop_columns:
            out.pop(k, None)
        return out

    def descriptors(self, config):
        b = self.bio
        return {'fp_MW': b['fp_MW'], 'pp_MW': b['pp_MW'],
                'fs_F_mass': b['capacity'], 'fs_price': b['fs_price'],
                'fs_w_fermentable_sugar_eq': b['w_sugar'], 'fs_w_lignin': b['w_lignin']}

    def versions(self):
        return {'synthetic': True, 'numpy': np.__version__}


# %% Generation and ground truth

def generate(out_dir, names=None, n_random=6, n_productivity=2, n_y=5, n_t=5, seed=1,
             kappa=0., fail_every=23, sample_limit=None, tag='', harness_version=2):
    """Write ``{out_dir}/{name}/{name}{tag}.jsonl`` for each synthetic
    biorefinery via the real harness (with the :func:`biosteam_stub`).
    ``harness_version=1`` strips the v2 columns (:func:`strip_v2`). Returns
    the list of paths."""
    paths = []
    with biosteam_stub():
        for name in names or SYNTHETIC_BIOREFINERIES:
            ad = SyntheticAdapter(name, kappa=kappa, fail_every=fail_every)
            d = os.path.join(out_dir, name)
            paths.append(run_campaign(ad, d, n_random=n_random, n_productivity=n_productivity,
                                      seed=seed, n_y=n_y, n_t=n_t, point_timeout=0,
                                      sample_limit=sample_limit, tag=tag, verbose=False))
    if harness_version == 1:
        for p in paths:
            strip_v2(p)
    return paths


def _is_v2_key(k):
    return k in V2_POINT_KEYS or k.endswith('__installed_cost_osbl') or k.startswith('fs_')


def strip_v2(path):
    """Rewrite a JSONL file as a legacy (v1) harness file: drop the v2 point
    columns (natural gas per pool, labor, OSBL installed cost, per-point
    ``fs_*``) and ``harness_version`` / ``fs_*`` of ``sample_start``."""
    with open(path) as fh:
        lines = fh.readlines()
    out = []
    for line in lines:
        try:
            rec = json.loads(line)
        except Exception:
            out.append(line)
            continue
        if rec.get('record') in ('point', 'sample_start'):
            rec = {k: v for k, v in rec.items() if k != 'harness_version' and not _is_v2_key(k)}
        out.append(json.dumps(rec) + '\n')
    with open(path, 'w') as fh:
        fh.writelines(out)
    return path


def ground_truth(paths, kappa=0.):
    """True phi, F and channel/intensity coefficients (primed basis, same OLS
    and the same feasible points / g as the pipeline) for every complete
    sample in the given synthetic JSONL files."""
    rows = []
    for path in paths:
        starts, pts, done = {}, {}, set()
        with open(path) as fh:
            for line in fh:
                try:
                    rec = json.loads(line)
                except Exception:
                    continue
                kind = rec.get('record')
                if kind == 'sample_start':
                    starts[rec['sample_id']] = rec
                    pts[rec['sample_id']] = []
                elif kind == 'point' and rec.get('status') == 'ok':
                    pts.setdefault(rec['sample_id'], []).append(rec)
                elif kind == 'sample_end':
                    done.add(rec['sample_id'])
        for sid in done:
            P = pts.get(sid, [])
            if not P:
                continue
            name = P[0]['biorefinery']
            bio = SYNTHETIC_BIOREFINERIES[name]
            cfg = next(dict(c) for c in bio['configs'] if c['config_id'] == P[0]['config_id'])
            params = {k: v for k, v in starts[sid].items() if k.startswith(('p_', 'bp_'))}
            y = np.array([p['y'] for p in P])
            t = np.array([p['t'] for p in P])
            res = [point_model(bio, cfg, params, p['y'], p['t'], p['r'], kappa) for p in P]
            g = physics.median_recovery([r['recovery'] for r in res])
            denom = np.array([r['F_final_product'] * r['h'] for r in res])
            m = {}
            for ch in physics.OPEX_CHANNELS:
                m[ch] = np.array([r['E'][ch] for r in res]) / denom
            phi = res[0]['phi']
            # capital by FCI contribution when the OSBL split is recorded (v2)
            osbl = bool(bio.get('osbl')) and any(k.endswith('__installed_cost_osbl') for k in P[0])
            for ch, areas in physics.CAPEX_CHANNEL_AREAS.items():
                ws = [area_fci_weights(r, osbl) for r in res]
                m[ch] = np.array([phi * r['FCI'] / sum(w.values()) * sum(w[a] for a in areas)
                                  for r, w in zip(res, ws)]) / denom
            m['fixed'] = np.array([r['F'] for r in res]) / denom
            m['closure'] = np.array([r['eps'] for r in res]) / denom
            q = {'fs': np.array([r['F_feedstock'] / r['F_final_product'] for r in res]),
                 'el': np.array([(r['electricity_consumption'] - r['electricity_production'])
                                 * 1e3 / r['F_final_product'] for r in res]),
                 'ng': np.array([r['natural_gas_flow'] / r['F_final_product'] for r in res])}
            cols = [m[ch] for ch in physics.CHANNELS] + [q[k] for k in physics.INTENSITY_KEYS]
            C, _ = physics.fit_ols_many(y, t, np.column_stack(cols))
            C = g * C
            row = {'sample_id': sid, 'biorefinery': name, 'phi_true': phi,
                   'F_fixed_true': res[0]['F'], 'g_true': g,
                   'w_sugar': bio['w_sugar'], 'eta': bio['eta']}
            for j, ch in enumerate(physics.CHANNELS):
                for i, k in enumerate(physics.COEFF_KEYS):
                    row[f'ch__{ch}__{k}'] = C[i, j]
            for j, qk in enumerate(physics.INTENSITY_KEYS):
                for i, k in enumerate(physics.COEFF_KEYS):
                    row[f'int__{qk}__{k}'] = C[i, len(physics.CHANNELS) + j]
            rows.append(row)
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--out', required=True)
    ap.add_argument('--names', nargs='*', default=None)
    ap.add_argument('--n-random', type=int, default=6)
    ap.add_argument('--n-productivity', type=int, default=2)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--kappa', type=float, default=0.)
    ap.add_argument('--harness-version', type=int, default=2, choices=(1, 2))
    a = ap.parse_args(argv)
    paths = generate(a.out, a.names, a.n_random, a.n_productivity, seed=a.seed, kappa=a.kappa,
                     harness_version=a.harness_version)
    gt = ground_truth(paths, a.kappa)
    gt.to_csv(os.path.join(a.out, 'ground_truth.csv.gz'), index=False)
    print(f'wrote {len(paths)} files, {len(gt)} complete samples -> {a.out}')


if __name__ == '__main__':
    main()
