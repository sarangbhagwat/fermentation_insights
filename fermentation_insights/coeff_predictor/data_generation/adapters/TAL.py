# -*- coding: utf-8 -*-
"""
Data-generation adapter for the TAL biorefineries of Bioindustrial-Park
(``biorefineries.TAL``): triacetic acid lactone (TAL) by aerobic batch
fermentation (Y. lipolytica, not neutralized; NaOH only controls the pH by
neutralizing byproduct acids), recovered by the "solubility-exploit" process,
and optionally upgraded catalytically to potassium sorbate (KSA).

Configurations (as in the fermentation_insights publication):
``TAL_<feedstock>`` (final product TAL) and ``TAL_SA_<feedstock>`` (final
product potassium sorbate, MPSP in USD/kg pure KSA) for feedstock in
{glucose, corn, sugarcane, cornstover}.

Definitions follow the biorefinery's TRY-analysis scripts
(``TAL/analyses/fermentation/TRY_analysis_{TAL,SA}_FGI.py``):

* parameters: ``parameter-distributions_{TAL|Sorbate}_A_FGI_<feedstock>.xlsx``;
* MPSP = ``tea.solve_price(product) / purity`` [USD/kg pure TAL or KSA];
* recovery = ``product.imass[TAL|KSA] / R302.outs[1].imass['TAL']`` (mass
  basis; includes the TAL -> KSA mass change for KSA);
* yields are converted to the spec's fraction of theoretical
  (0.4667 g TAL/g glucose); TRY domain: 5.3% of theoretical to the inflection
  yield (1 - cell mass - citrate conversion = 57.24% of theoretical) and
  10-100 g/L; baseline productivity 0.12 g/L/h;
* each point: ``spec.load_specifications(y, t, r)`` followed by 3 system
  simulations (as in the TRY scripts), plus further simulations until the
  MPSP changes by < 0.1%.

Requires the fi-TAL branch of Bioindustrial-Park (biosteam 2.54 port).

Each configuration builds module-level systems on a shared flowsheet, so only
ONE configuration can be loaded per Python process; run configurations in
separate processes (``--configs <config_id>``; see scripts/run_TAL.sh).
"""
from __future__ import annotations

import importlib
import os
import time
import warnings

import numpy as np

from ..harness import BiorefineryAdapter, InfeasibleSpecification, Parameter

__all__ = ('Adapter',)

_FEEDSTOCKS = {
    # feedstock: (feedstock_type, feedstock stream ID)
    'glucose': ('pure_sugar', 'glucose_feedstock'),
    'corn': ('starch', 'corn'),
    'sugarcane': ('sugar_juice', 'sugarcane'),
    'cornstover': ('lignocellulose', 'cornstover'),
}

# Canonical names for parameters of the FGI parameter-distribution files
_CANONICAL = {
    'Plant annual operating days': 'p_operating_days',
    'Feedstock unit price': 'p_feedstock_price',
    'Natural gas unit price': 'p_natural_gas_price',
    'Electricity unit price': 'p_electricity_price',
    'Feedstock capacity': 'p_feedstock_capacity',
    'Federal corporate tax rate': 'p_income_tax',
    'Internal rate of return': 'p_IRR',
    'Fermentation Y. lipolytica yield': 'p_cell_mass_yield',
    'Boiler efficiency': 'p_boiler_efficiency',
    'Turbogenerator efficiency': 'p_turbogenerator_efficiency',
}

# Parameters of the distribution files that have no effect on the model and
# are therefore excluded: a dummy, and 'Fermentation acetate loading', whose
# load statement sets `R302.acetate_loading`, an attribute the fermenter never
# reads (the acetate spiking target is `R302.acetate_target_loading`).
_EXCLUDED = ('Blank parameter', 'Fermentation acetate loading')

_AREA_MAP = {
    'feedstock acquisition': 'feedstock',
    'feedstock juicing': 'conversion',
    'feedstock pretreatment': 'conversion',
    'feedstock pretreatment and saccharification': 'conversion',
    'pretreatment': 'conversion',
    'corn processing': 'conversion',
    'fermentation': 'fermentation',
    'separation': 'separation',
    'upgrading': 'upgrading',
    'wastewater': 'wastewater',
    'storage & other facilities': 'facilities',
    'boiler & turbogenerator': 'chp',
    'cooling utility facilities': 'facilities',
    'other facilities': 'facilities',
    'heat exchanger network': 'facilities',
    'unassigned units': 'facilities',
    # metric-only pseudo-groups of the TAL module (no units)
    'natural gas (for steam generation)': None,
    'natural gas (for product drying)': None,
    'chilled brine': None,
    'fixed operating cost': None,
    'electricity consumption': None,
    'heating duty': None,
    'excess electricity': None,
}

# Storage/feeding units of fermentation nutrients (TAL groups them with
# 'storage & other facilities'); moved to the fermentation area so that the
# CSL, acetate, and DAP purchase costs are attributed to fermentation.
_FERMENTATION_NUTRIENT_UNITS = ('T601', 'T601_P', 'T607', 'T607_P', 'T609', 'T609_P')

_THEORETICAL_YIELD_FRACTION_MIN = 0.053  # of theoretical, TRY_analysis_*_FGI.py
_T_MIN, _T_MAX = 10., 100.  # g/L, TRY_analysis_*_FGI.py

# Maximum feasible titer [g/L] vs. absolute yield y [g/g]: t_max = k*y + c,
# calibrated with this adapter by bisection of the sugar-concentration
# feasibility frontier at baseline parameters (y = 5.3, 13 and 25% of
# theoretical; frontier is linear in y); max_titer stays 0.5% below it.
# (For cornstover the frontier at y_min is below t_min = 10 g/L.)
_MAX_TITER_LINE = {
    'TAL_glucose': (1250., -1.9), 'TAL_SA_glucose': (1250., -1.9),
    'TAL_corn': (656., -0.1), 'TAL_SA_corn': (656., -0.1),
    'TAL_sugarcane': (678., -1.4), 'TAL_SA_sugarcane': (614., -0.6),
    'TAL_cornstover': (421., -1.1), 'TAL_SA_cornstover': (421., -1.1),
}

_loaded_module = None  # module-level systems: one configuration per process


def _is_infeasible_error(e):
    msg = str(e).lower()
    name = type(e).__name__
    return (name == 'InfeasibleRegion' or 'sugar concentration' in msg
            or 'infeasible' in msg)


def _is_bracket_error(e):
    msg = str(e).lower()
    return ('opposite sign' in msg or 'same sign' in msg or 'bracket' in msg)


class Adapter(BiorefineryAdapter):
    name = 'TAL'
    area_map = _AREA_MAP

    #: Minimum and maximum number of system simulations per TRY point
    #: (the TRY scripts use 3) and relative MPSP tolerance between the last
    #: two simulations.
    min_simulations = 3
    max_simulations = 8
    MPSP_rtol = 1e-3

    def __init__(self):
        self.config = None
        self._parameters = None
        self._n_simulations = 0
        self._MPSP_rel_change = float('nan')

    # %% Configurations

    def configurations(self):
        configs = []
        for product in ('TAL', 'TAL_SA'):
            for feedstock, (ftype, _) in _FEEDSTOCKS.items():
                upgrading = product == 'TAL_SA'
                configs.append(dict(
                    config_id=f'{product}_{feedstock}',
                    feedstock=feedstock,
                    feedstock_type=ftype,
                    fermentation_product='TAL',
                    final_product='KSA' if upgrading else 'TAL',
                    ferm_product_chemical='TAL',
                    final_product_chemical='PotassiumSorbate' if upgrading else 'TAL',
                    neutralization=False,  # TAL is not neutralized
                    neutralizing_agent='',
                    pH_control_base='NaOH',  # pH 6.5 by neutralizing byproduct acids (citrate, acetate, phosphate)
                    n_acid_groups=0,  # lactone; acidic enol (vinylogous acid), no carboxylic group
                    pKa1=4.7,
                    final_product_stoich_mass_ratio=(150.2169 / 126.1100) if upgrading else 1.,
                    upgrading=upgrading,
                    separation_variant='',
                    separation_principle=(
                        'solubility exploit: thermal decarboxylation of dissolved TAL, '
                        'evaporation, cooling crystallization, centrifugation, drying'
                        + ('; upgrading: hydrogenation (Ni-SiO2), dehydration (Amberlyst-70), '
                           'ring-opening hydrolysis with KOH, KSA crystallization from IPA, drying'
                           if upgrading else '')),
                    aerobic=True,
                    fermentation_mode='batch',
                    T_fermentation=28. + 273.15,
                    product_water_solubility=float('nan'),  # set in load (TAL solubility model)
                ))
        return configs

    # %% Loading

    def load(self, config):
        global _loaded_module
        warnings.filterwarnings('ignore')
        feedstock = config['feedstock']
        upgrading = config['final_product'] == 'KSA'
        modname = ('biorefineries.TAL.models.models_'
                   f'{"SA" if upgrading else "TAL"}_solubility_exploit_{feedstock}')
        if _loaded_module is not None and _loaded_module != modname:
            raise RuntimeError(
                f'{_loaded_module} is already loaded in this process; TAL configurations '
                'use module-level systems on a shared flowsheet and must be run in '
                'separate processes (use --configs <config_id>)')
        import biosteam as bst
        from biorefineries import TAL
        models = importlib.import_module(modname)
        _loaded_module = modname
        self.config = config
        self.models = models
        self.bst = bst
        f = models.TAL_sys.flowsheet
        self.u, self.s = u, s = f.unit, f.stream
        self.system = models.TAL_sys
        self.tea = models.TAL_tea
        self.lca = models.TAL_lca
        self.spec = models.spec
        self.model = models.TAL_model
        self.R302 = u.R302
        self.broth = u.R302.outs[1]
        self.product_ID = 'KSA' if upgrading else 'TAL'
        self.product_stream = s.KSA_product if upgrading else s.TAL_product
        self.feedstock_stream = models.feedstock
        self.BT = u.BT701
        self._original_converge_method = getattr(self.system, 'converge_method', None)
        chems = self.product_stream.chemicals
        self.theoretical_yield = 2 * chems.TAL.MW / (3 * chems.Glucose.MW)  # g TAL / g glucose (= xylose basis)
        self._glucose_eq_sucrose = 2 * chems.Glucose.MW / chems.Sucrose.MW if 'Sucrose' in chems else 0.

        # Parameter distributions used by the TRY analysis of the publication
        fn = os.path.join(os.path.dirname(TAL.__file__), 'analyses', 'full', 'parameter_distributions',
                          f'parameter-distributions_{"Sorbate" if upgrading else "TAL"}_A_FGI_{feedstock}.xlsx')
        self.parameter_distributions_file = fn
        self.model.parameters = ()
        self.model.load_parameter_distributions(fn, models.namespace_dict)
        self._parameters = None

        # Unit groups for the harness (copies; TAL's own groups are untouched)
        self.unit_groups = self._make_unit_groups(models.unit_groups)

        # Product water solubility from TAL's fitted solubility model [g/L water, 25 C]
        try:
            from biorefineries.TAL.models.solubility.fit_TAL_solubility_in_water_one_parameter_van_laar_activity import \
                get_TAL_solubility_in_water_gpL
            config['product_water_solubility'] = float(get_TAL_solubility_in_water_gpL(298.15))
        except Exception:
            config['product_water_solubility'] = 8.6

        # Baseline state as in the TRY scripts: baseline parameter values, baseline TRY
        for p in self.model.parameters:
            p.setter(p.baseline)
        self.baseline_spec = dict(spec_1=self.spec.baseline_yield, spec_2=self.spec.baseline_titer,
                                  spec_3=self.spec.baseline_productivity)
        self.spec.load_specifications(**self.baseline_spec)
        self._simulate_converged()
        self.baseline_MPSP = self._MPSP()

    def _make_unit_groups(self, groups):
        UnitGroup = self.bst.UnitGroup
        moved = [self.u[i] for i in _FERMENTATION_NUTRIENT_UNITS if i in self.u]
        moved_set = set(moved)
        new = []
        for g in groups:
            units = [i for i in g.units if i not in moved_set]
            if g.name == 'fermentation':
                units = units + [i for i in moved if i not in units]
            new.append(UnitGroup(g.name, units=units))
        covered = set(sum([list(g.units) for g in new], []))
        missing = [i for i in self.system.units if i not in covered]
        # e.g., the supernatant recycle splitter S403 of the KSA systems
        separation = [i for i in missing if i.ID[:2] == 'S4' or i.ID[:2] == 'F4' or i.ID[:2] == 'H4']
        for g in new:
            if g.name == 'separation' and separation:
                g.units.extend(separation)
        missing = [i for i in missing if i not in separation]
        if missing:
            new.append(UnitGroup('unassigned units', units=missing))
        self.unassigned_units = [i.ID for i in missing]
        self.units_added_to_separation = [i.ID for i in separation]
        for g in new:
            if g.name not in _AREA_MAP:
                warnings.warn(f'TAL adapter: unit group {g.name!r} not in area_map (-> facilities)')
        return new

    # %% Domain

    def domain(self):
        R302 = self.R302
        theo = self.theoretical_yield
        inflection = theo * (1. - R302.regular_microbe_conversion - R302.regular_citric_acid_conversion)
        return dict(
            theoretical_yield=theo,
            inflection_yield=inflection,
            y_min=_THEORETICAL_YIELD_FRACTION_MIN * theo,
            y_max=inflection,
            t_min=_T_MIN, t_max=_T_MAX,
            baseline_productivity=self.spec.baseline_productivity,
            baseline_yield=self.spec.baseline_yield * theo,
            baseline_titer=self.spec.baseline_titer,
        )

    def max_titer(self, y):
        k, c = _MAX_TITER_LINE[self.config['config_id']]
        return 0.995 * (k * y + c)

    # %% Parameters

    def parameters(self):
        if self._parameters is None:
            params = []
            for p in self.model.parameters:
                if p.name in _EXCLUDED:
                    continue
                prm = Parameter.from_biosteam(p, name=_CANONICAL.get(p.name))
                prm.description = f'{p.element}: {p.name} [{p.units}]'
                params.append(prm)
            self._parameters = params
        return self._parameters

    # %% Simulation

    def set_try(self, y, t, r):
        y_frac = y / self.theoretical_yield
        if not (0. < y_frac < 1.):
            raise InfeasibleSpecification(f'yield {y:.4g} g/g outside (0, theoretical)')
        self._target = dict(spec_1=y_frac, spec_2=t, spec_3=r)
        try:
            self.spec.load_specifications(**self._target)
        except Exception as e:
            if _is_infeasible_error(e):
                raise InfeasibleSpecification(f'{type(e).__name__}: {e}')
            if _is_bracket_error(e):
                raise InfeasibleSpecification(f'titer not attainable: {type(e).__name__}: {e}')
            raise

    def _MPSP(self):
        product = self.product_stream
        purity = product.imass[self.product_ID] / product.F_mass
        return self.tea.solve_price(product) / purity

    def _simulate_converged(self):
        system = self.system
        last = None
        rel = float('nan')
        n = 0
        while True:
            system.simulate()
            n += 1
            if n >= self.min_simulations - 1:
                MPSP = self._MPSP()
                if last is not None:
                    rel = abs(MPSP - last) / abs(MPSP)
                    if (n >= self.min_simulations and rel < self.MPSP_rtol) or n >= self.max_simulations:
                        break
                last = MPSP
        self._n_simulations = n
        self._MPSP_rel_change = rel

    def _run_evaporators(self):
        from biosteam.units import MultiEffectEvaporator
        for unit in self.system.units:
            if isinstance(unit, MultiEffectEvaporator):
                unit._V_first_effect = None
                try: unit._run()
                except Exception: pass

    def _bugfix_barrage(self):
        """Reset/reload and solver switching as in the biorefinery's model
        specification (models_*.run_bugfix_barrage); the original
        convergence method is restored afterwards."""
        system, spec = self.system, self.spec
        method = self._original_converge_method
        target = getattr(self, '_target', self.baseline_spec)
        try:
            system.reset_cache(); system.empty_recycles(); self._run_evaporators()
            spec.load_specifications(**self.baseline_spec)
            system.simulate()
            spec.load_specifications(**target)
            self._simulate_converged()
            return
        except Exception as e:
            if _is_infeasible_error(e): raise
            last_error = e
        for solver in ('fixedpoint', 'aitken'):
            try:
                system.reset_cache(); system.empty_recycles(); self._run_evaporators()
                system.converge_method = solver
                spec.load_specifications(**target)
                self._simulate_converged()
                return
            except Exception as e:
                if _is_infeasible_error(e): raise
                last_error = e
            finally:
                if method is not None: system.converge_method = method
        raise last_error

    def simulate(self):
        try:
            self._simulate_converged()
        except Exception as e:
            if _is_infeasible_error(e):
                raise InfeasibleSpecification(f'{type(e).__name__}: {e}')
            try:
                self._bugfix_barrage()
            except Exception as e2:
                if _is_infeasible_error(e2):
                    raise InfeasibleSpecification(f'{type(e2).__name__}: {e2}')
                raise

    def reset(self):
        system = self.system
        if self._original_converge_method is not None:
            system.converge_method = self._original_converge_method
        system.reset_cache()
        system.empty_recycles()
        self._run_evaporators()
        self.spec.load_specifications(**self.baseline_spec)
        for i in range(2): system.simulate()

    # %% Collection

    def collect(self):
        u, s = self.u, self.s
        product, ID = self.product_stream, self.product_ID
        broth, R302, BT = self.broth, self.R302, self.BT
        system, tea = self.system, self.tea
        bst = self.bst
        hours = tea.operating_hours

        purity = product.imass[ID] / product.F_mass
        price = tea.solve_price(product)
        product.price = price  # as in the biorefinery's get_MSP
        MPSP = price / purity
        F_final = product.imass[ID]
        F_TAL_broth = broth.imass['TAL']

        # Fermentable sugars (glucose-eq.) entering the fermentation area
        sugar_feed = u.S302.ins[0]
        F_sugar = sugar_feed.imass['Glucose'] + sugar_feed.imass['Xylose']
        if self._glucose_eq_sucrose: F_sugar += self._glucose_eq_sucrose * sugar_feed.imass['Sucrose']

        # Fermenter volume
        D = R302.design_results
        V_fermenters = D.get('Reactor volume', float('nan')) * D.get('Number of reactors', float('nan'))

        # Utilities of process units (facilities such as boiler, cooling
        # tower, chilled water packages, and the HXN are excluded)
        process_units = [i for i in system.units if not isinstance(i, bst.Facility)]
        heating = cooling = 0.
        for unit in process_units:
            for hu in unit.heat_utilities:
                if hu.duty > 0: heating += hu.duty
                else: cooling -= hu.duty
        HXN = u.HXN1001 if 'HXN1001' in u else None
        HXN_heating = HXN_cooling = 0.
        if HXN is not None:
            for hu in HXN.heat_utilities:
                if hu.duty > 0: HXN_heating += hu.duty
                else: HXN_cooling -= hu.duty

        # Natural gas: boiler + direct-fired dryers
        NG_BT = BT.natural_gas.F_mass
        NG_total = sum([i.F_mass for i in system.feeds if 'natural_gas' in i.ID])
        base = R302.ins[6]

        def cost(ID):
            return s[ID].cost if ID in s else 0.

        try: GWP = self.lca.GWP
        except Exception: GWP = float('nan')
        try: FEC = self.lca.FEC
        except Exception: FEC = float('nan')
        try: HXN_error = HXN.energy_balance_percent_error
        except Exception: HXN_error = float('nan')

        theo = self.theoretical_yield
        out = dict(
            MPSP=MPSP,
            recovery=F_final / F_TAL_broth,
            GWP=GWP, FEC=FEC,
            F_final_product=F_final,
            F_ferm_product=F_TAL_broth,
            F_sugar_fed=F_sugar,
            F_feedstock=self.feedstock_stream.F_mass,
            feedstock_price=self.feedstock_stream.price,
            V_broth=broth.F_vol,
            F_water_broth=broth.imass['Water'],
            V_fermenters=V_fermenters,
            tau_fermentation=R302.tau,
            base_cost=base.cost,  # NaOH for fermentation pH control (no product neutralization)
            acid_cost=0.,  # no acidulation
            heating_duty=heating / 1e6,
            cooling_duty=cooling / 1e6,
            electricity_consumption=system.power_utility.consumption / 1e3,
            electricity_production=system.power_utility.production / 1e3,
            natural_gas_flow=NG_total,
            electricity_price=bst.PowerUtility.price,
            natural_gas_price=BT.natural_gas_price,
            # adapter-specific
            x_purity=purity,
            x_MPSP_impure=price,
            x_n_simulations=self._n_simulations,
            x_MPSP_rel_change=self._MPSP_rel_change,
            x_yield_fraction_theoretical=self.spec.spec_1,
            x_inflection_yield=theo * (1. - R302.regular_microbe_conversion - R302.regular_citric_acid_conversion),
            x_titer_measured=F_TAL_broth / broth.F_vol,
            x_yield_measured=F_TAL_broth / F_sugar if F_sugar else float('nan'),
            x_F_acetate_fed=(R302.ins[3].imass['AceticAcid'] + sugar_feed.imass['AceticAcid']),
            x_F_cell_mass=broth.imass['FermMicrobe'],
            x_F_citrate=broth.imass['CitricAcid'] + broth.imass['SodiumCitrate'],
            x_CO2_vent=R302.outs[0].imass['CO2'],
            x_F_air=R302.ins[5].F_mass,
            x_F_base_pH_control=base.F_mass,
            x_CSL_cost=cost('CSL_fresh'),
            x_DAP_cost=cost('DAP_fresh'),
            x_acetate_cost=cost('Acetate_fresh'),
            x_decarboxylation_base_cost=cost('base_decarboxylation_fresh'),
            x_N_fermenters=D.get('Number of reactors', float('nan')),
            x_F301_evaporated=(u.F301.outs[1].F_mass if 'F301' in u else 0.),
            x_dilution_water=(u.M304.ins[1].F_mass if 'M304' in u and len(u.M304.ins) > 1 else float('nan')),
            x_natural_gas_BT=NG_BT,
            x_natural_gas_dryers=NG_total - NG_BT,
            x_HXN_heating_duty=HXN_heating / 1e6,
            x_HXN_cooling_duty=HXN_cooling / 1e6,
            x_HXN_energy_balance_error=HXN_error,
            x_annual_production=F_final * hours / 1e3,  # metric ton/y
        )
        return out
