# -*- coding: utf-8 -*-
"""
Adapter for the 3-hydroxypropionic acid (3-HP) -> acrylic acid biorefineries
of Bioindustrial-Park (``biorefineries.HP``).

Configurations (16, as in the fermentation_insights publication)::

    {glucose, corn, sugarcane, cornstover} x
    {'HP'                  : low-pH fermentation (no neutralization), improved separations,
     'HP_neutral'          : fermentation with lime neutralization, improved separations,
     'HP_hexanol'          : low-pH fermentation, hexanol-extraction separation,
     'HP_neutral_hexanol'  : neutralization, hexanol-extraction separation}

Models and specifications follow the publication's scripts
(``MPSP_TRY_fit/uncertainties_coeffs/HP*_uncertainties_coeffs.py``) and the TRY
scripts in ``HP/analyses/fermentation``:

* systems are loaded from ``biorefineries.HP.models.<feedstock>.models_*``
  (``*_improved_separations`` or ``*_hexanol``);
* ``spec.reactor.neutralization`` is set per configuration;
* parameter distributions: ``parameter-distributions_<feedstock>_Acrylic_300L_FGI.xlsx``
  (the publication's "fixed gross input" mode: feedstock capacity is an
  uncertain parameter and is held constant across TRY; no
  ``set_production_capacity``). These files are current on Bioindustrial-Park
  master and (checked for cornstover) identical to the non-FGI '300L' files
  except that yield, titer, productivity and desired annual production are
  replaced by feedstock capacity;
* TRY: ``spec.load_specifications(spec_1=y/theoretical_yield, spec_2=t,
  spec_3=r)`` then ``system.simulate()``, with the scripts' bugfix barrage;
* MPSP = ``tea.solve_price(AcrylicAcid) / purity`` [USD/kg pure AA];
* recovery = ``get_product_recover_FGI`` = kg AA / kg 3-HP-eq. in broth.

IMPORTANT: the HP systems are module-level singletons registered in the main
flowsheet. Different feedstocks/separation variants cannot be loaded in the
same process; run each configuration in its own process (``--configs``).
``load`` raises if asked to switch to a different system module.
HP and HP_neutral of the same feedstock/separation share a module and *can*
run in one process (neutralization is toggled and the baseline re-simulated),
but the launcher still runs one configuration per process.
"""
from __future__ import annotations

import importlib
import inspect
import math
import os

import numpy as np

from ..harness import BiorefineryAdapter, InfeasibleSpecification, Parameter

__all__ = ('Adapter',)

_MW_HP = 90.07794
_MW_AA = 72.06266
_MW_glucose = 180.15588
_MW_sucrose = 342.29648

_FEEDSTOCKS = {
    # feedstock: (models subpackage, improved-separations module, hexanol module, feedstock_type)
    'glucose': ('glucose', 'models_glucose_improved_separations', 'models_glucose_hexanol', 'pure_sugar'),
    'corn': ('corn', 'models_corn_improved_separations', 'models_corn_hexanol', 'starch'),
    'sugarcane': ('sugarcane', 'models_sc_improved_separations', 'models_sc_hexanol', 'sugar_juice'),
    'cornstover': ('cornstover', 'models_cs_improved_separations', 'models_cs_hexanol', 'lignocellulose'),
}

_VARIANTS = {
    # variant: (neutralization, hexanol)
    'HP': (False, False),
    'HP_neutral': (True, False),
    'HP_hexanol': (False, True),
    'HP_neutral_hexanol': (True, True),
}

#: Error substrings that the HP TRY scripts treat as infeasible TRY combinations.
_INFEASIBLE_KEYS = ('sugar concentration', 'permeate moisture', 'opposite sign')

#: Canonical names for the uncertain parameters of the 300L_FGI distribution files.
_CANONICAL = {
    'Plant annual operating days': 'p_operating_days',
    'Feedstock unit price': 'p_feedstock_price',
    'Natural gas unit price': 'p_natural_gas_price',
    'Electricity unit price': 'p_electricity_price',
    'Lime unit price': 'p_base_price',
    'Sulfuric acid unit price': 'p_acid_price',
    'Feedstock capacity': 'p_feedstock_capacity',
    'Federal corporate tax rate': 'p_income_tax',
    'Internal rate of return': 'p_IRR',
    'Fermentation I. orientalis yield': 'p_cell_mass_yield',
    'Boiler efficiency': 'p_boiler_efficiency',
    'Turbogenerator efficiency': 'p_turbogenerator_efficiency',
}
_EXCLUDED_PARAMETERS = ('Blank parameter', 'Fermentation HP yield',
                        'Fermentation HP titer', 'Fermentation HP productivity',
                        'Desired annual production')

_AREA_MAP = {
    'feedstock acquisition': 'feedstock',
    'feedstock juicing': 'conversion',  # glucose dissolution / corn liquefaction & saccharification / cane juicing
    'feedstock pretreatment and saccharification': 'conversion',
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
    # metric-only groups without units (would contribute 0 anyway)
    'natural gas (for steam generation)': None,
    'natural gas (for product drying)': None,
    'chilled brine': None,
    'fixed operating cost': None,
    'electricity consumption': None,
    'heating duty': None,
    'excess electricity': None,
}


def _is_infeasible(e):
    s = str(e).lower()
    return any(k in s for k in _INFEASIBLE_KEYS)


class Adapter(BiorefineryAdapter):
    name = 'HP'
    area_map = _AREA_MAP

    #: Safety factor applied to the computed maximum attainable titer.
    max_titer_safety_factor = 0.99
    #: Maximum tolerated relative error between achieved and specified titer.
    titer_rtol = 0.05
    #: Distribution file mode (publication: '300L_FGI').
    distribution_mode = '300L_FGI'

    def __init__(self):
        self._module_name = None
        self.models = None
        self.config = None

    # %% Configurations

    def configurations(self):
        configs = []
        for variant, (neutral, hexanol) in _VARIANTS.items():
            for feedstock, (_, _, _, ftype) in _FEEDSTOCKS.items():
                configs.append(dict(
                    config_id=f'{variant}_{feedstock}',
                    feedstock=feedstock,
                    feedstock_type=ftype,
                    fermentation_product='HP',
                    final_product='AcrylicAcid',
                    ferm_product_chemical='HP',
                    final_product_chemical='AcrylicAcid',
                    neutralization=neutral,
                    neutralizing_agent='Ca(OH)2' if neutral else '',
                    n_acid_groups=1,
                    pKa1=4.51,
                    final_product_stoich_mass_ratio=_MW_AA / _MW_HP,
                    upgrading=True,  # catalytic (TiO2) dehydration of 3-HP to acrylic acid
                    separation_variant='hexanol_extraction' if hexanol else '',
                    separation_principle=(
                        'acidulation+evaporation+LLE(hexanol)+distillation' if hexanol else
                        'acidulation+adsorption+ion_exchange+evaporation'
                    ) + '|dehydration+distillation',
                    aerobic=True,
                    fermentation_mode='batch',
                    T_fermentation=28. + 273.15,
                    product_water_solubility=1000.,  # 3-HP is miscible with water
                ))
        return configs

    # %% Loading

    def load(self, config):
        feedstock = config['feedstock']
        neutral, hexanol = _VARIANTS[config['config_id'][:-len(feedstock) - 1]]
        subpkg, mod_improved, mod_hexanol, _ = _FEEDSTOCKS[feedstock]
        mod_name = f'biorefineries.HP.models.{subpkg}.{mod_hexanol if hexanol else mod_improved}'
        if self._module_name is not None and self._module_name != mod_name:
            raise RuntimeError(
                f'HP adapter: {self._module_name} is already loaded in this process; '
                f'{mod_name} (config {config["config_id"]}) must run in its own '
                f'process (use --configs with one configuration per process).')
        first_load = self._module_name is None
        models = importlib.import_module(mod_name)
        self._module_name = mod_name
        self.models = models
        self.config = config
        self.hexanol = hexanol
        import biosteam as bst
        self.bst = bst
        from biorefineries.HP import units as HP_units
        self.system = system = models.HP_sys
        self.tea = models.HP_tea
        self.lca = models.HP_lca
        self.spec = spec = models.spec
        self.u = u = models.u
        self.s = models.s
        self.feedstock_stream = models.feedstock
        self.product_stream = models.product_stream
        self.R302 = R302 = spec.reactor
        self.R303 = spec.seed_train
        self.broth = R302.outs[1]
        self.theoretical_yield = float(models.theoretical_max_g_HP_per_g_glucose)
        self._default_converge_method = (system.converge_method if first_load
                                         else self._default_converge_method)
        # Neutralization (publication: spec.reactor.neutralization = 'neutral' in tag)
        R302.neutralization = neutral

        # Units referenced in collect
        self.BT = u.BT701
        self.HXN = u.HXN1001
        self.S302 = u.S302
        self.F301 = u.F301
        self.M304 = u.M304
        self.acidulation_reactors = [i for i in system.units
                                     if isinstance(i, HP_units.AcidulationReactor)]
        self.separation_evaporators = [
            i for i in models.unit_groups_dict['separation'].units
            if isinstance(i, bst.MultiEffectEvaporator)]
        self.evaporators = [i for i in system.units if isinstance(i, bst.MultiEffectEvaporator)]
        self._acid_feeds = self._find_acidulation_acid_feeds()
        self._gypsum = getattr(models.s, 'gypsum', None)
        self.utility_facilities = [i for i in system.units
                                   if isinstance(i, (bst.BoilerTurbogenerator, bst.CoolingTower,
                                                     bst.ChilledWaterPackage))]

        # Unit groups (HP's own + any unassigned units so areas sum to the TEA total)
        groups = [g for g in models.unit_groups]
        assigned = set()
        for g in groups:
            assigned.update(g.units)
        leftover = [i for i in system.units if i not in assigned]
        if leftover:
            groups.append(bst.UnitGroup('unassigned units', units=leftover))
        self.unit_groups = groups

        # Uncertain parameters (300L_FGI distributions, as in the publication)
        self._load_parameters()
        for p in self._bst_parameters:  # baseline values (as model.metrics_at_baseline)
            p.setter(p.baseline)

        # Baseline TRY and simulation
        self.baseline_spec = dict(spec_1=spec.baseline_yield, spec_2=spec.baseline_titer,
                                  spec_3=spec.baseline_productivity)
        self._simulate_baseline()
        self._compute_max_titer_curve()
        self._simulate_baseline()
        self.baseline_MPSP = self._MPSP()

    def _simulate_baseline(self):
        spec = self.spec
        b = self.baseline_spec
        self._requested = (b['spec_1'] * self.theoretical_yield, b['spec_2'], b['spec_3'])
        self.system.converge_method = self._default_converge_method
        spec.load_specifications(**b)
        try:
            self.system.simulate()
        except Exception as e:
            self._bugfix_barrage(e)

    def _find_acidulation_acid_feeds(self):
        """Purchased sulfuric acid feeds upstream of acidulation but not upstream
        of fermentation (i.e., excluding pretreatment acid)."""
        feeds = []
        ferm_upstream = set(self.R302.get_upstream_units()) | {self.R302}
        for R in self.acidulation_reactors:
            for unit in set(R.get_upstream_units()) | {R}:
                if unit in ferm_upstream:
                    continue
                for s in unit.ins:
                    if (s.source is None and s.price and 'H2SO4' in s.chemicals
                            and s not in feeds):
                        feeds.append(s)
        return feeds

    def _load_parameters(self):
        import biorefineries.HP as HP
        models = self.models
        model = models.HP_model
        product_folder = 'acrylic_acid_product'
        path = os.path.join(os.path.dirname(HP.__file__), 'analyses', 'full',
                            'parameter_distributions', product_folder,
                            f'parameter-distributions_{self.config["feedstock"]}_Acrylic_'
                            f'{self.distribution_mode}.xlsx')
        self.parameter_distributions_file = path
        model.parameters = ()
        model.load_parameter_distributions(path, models.namespace_dict)
        excluded = [i.lower() for i in _EXCLUDED_PARAMETERS]
        self._bst_parameters = [p for p in model.get_parameters()
                                if p.name.lower() not in excluded]
        params = []
        canonical = {k.lower(): v for k, v in _CANONICAL.items()}
        for p in self._bst_parameters:
            name = canonical.get(p.name.lower())  # biosteam lower-cases parameter names
            params.append(Parameter.from_biosteam(p, name=name))  # default name: bp_<slug>
        names = [p.name for p in params]
        assert len(set(names)) == len(names), names
        self._parameters = params

    # %% Domain

    def domain(self):
        theo = self.theoretical_yield
        R302 = self.R302
        inflection = (1. - R302.regular_biomass_conversion
                      - R302.regular_succinic_acid_conversion
                      - R302.regular_xylitol_conversion) * theo
        spec = self.spec
        return dict(
            theoretical_yield=theo,
            inflection_yield=inflection,
            y_min=0.05 * theo,
            y_max=inflection,
            t_min=10.,
            t_max=200.,
            baseline_productivity=float(spec.baseline_productivity),
            baseline_yield=float(spec.baseline_yield) * theo,
            baseline_titer=float(spec.baseline_titer),
        )

    def _compute_max_titer_curve(self):
        """Titer at the maximum attainable sugar concentration (pre-fermentation
        evaporator F301 at the vapor fraction used by the system's
        ``load_titer_with_glucose`` and no dilution water in M304) as a function
        of yield, evaluated by running the fermentation area only."""
        # The system module (where load_titer_with_glucose etc. are defined)
        sysmod = inspect.getmodule(self.spec.load_spec_2)
        load_titer = sysmod.load_titer_with_glucose
        default_V = inspect.signature(load_titer).parameters['set_F301_V'].default
        F301_lb = 0. if default_V is None else default_V
        R302 = self.R302
        dom = self.domain()
        ys = np.linspace(dom['y_min'], dom['y_max'], 5)
        ts = []
        spec = self.spec
        try:
            for y in ys:
                spec.load_spec_1(y / self.theoretical_yield)
                R302.titer_to_load = 1.
                sysmod.F301_titer_obj_fn(F301_lb)
                ts.append(sysmod.M304_titer_obj_fn(0.) + 1.)
            self._max_titer_ys = ys
            self._max_titer_ts = np.array(ts)
        except Exception as e:  # fall back to no bound
            print(f'HP adapter: could not compute max titer curve ({e})', flush=True)
            self._max_titer_ys = None

    def max_titer(self, y):
        if getattr(self, '_max_titer_ys', None) is None:
            return None
        ys, ts = self._max_titer_ys, self._max_titer_ts
        # linear in yield (titer ~ y x max sugar concentration); extrapolate linearly
        slope = (ts[-1] - ts[0]) / (ys[-1] - ys[0])
        if y < ys[0]:
            t = ts[0] + slope * (y - ys[0])
        elif y > ys[-1]:
            t = ts[-1] + slope * (y - ys[-1])
        else:
            t = float(np.interp(y, ys, ts))
        return self.max_titer_safety_factor * t

    # %% Parameters

    def parameters(self):
        return self._parameters

    # %% Simulation

    def set_try(self, y, t, r):
        spec_1 = y / self.theoretical_yield
        if not (0. < spec_1 < 1.):
            raise InfeasibleSpecification(f'yield {y} outside (0, theoretical)')
        self._requested = (y, t, r)
        self._load_error = None
        self._n_barrage = 0
        self.system.converge_method = self._default_converge_method
        try:
            self.spec.load_specifications(spec_1=spec_1, spec_2=t, spec_3=r)
        except Exception as e:
            if _is_infeasible(e):
                raise InfeasibleSpecification(str(e)[:200])
            self._load_error = e

    def simulate(self):
        system = self.system
        if self._load_error is not None:
            self._bugfix_barrage(self._load_error)
            return
        try:
            system.simulate()
        except Exception as e:
            str_e = str(e).lower()
            if 'sugar concentration' in str_e or 'permeate moisture' in str_e:
                raise InfeasibleSpecification(str(e)[:200])
            self._bugfix_barrage(e)

    # Bugfix barrage as in HP_uncertainties_coeffs.py / HP_hexanol_uncertainties_coeffs.py
    # (without set_production_capacity, i.e., fixed feedstock capacity).
    def _reset_evaporators(self):
        for unit in self.evaporators:
            unit._V_first_effect = None
        for unit in self.evaporators:
            try:
                unit._run()
            except Exception:
                pass

    def _reset_and_reload(self):
        system, spec = self.system, self.spec
        requested = getattr(self, '_requested', None)
        if requested is None:
            spec_1, spec_2, spec_3 = spec.spec_1, spec.spec_2, spec.spec_3
        else:  # spec attributes may be stale if load_specifications failed midway
            y, spec_2, spec_3 = requested
            spec_1 = y / self.theoretical_yield
        system.reset_cache()
        system.empty_recycles()
        if self.hexanol:
            self._reset_evaporators()
        spec.load_specifications(**self.baseline_spec)
        system.simulate()
        spec.load_specifications(spec_1=spec_1, spec_2=spec_2, spec_3=spec_3)
        system.simulate()

    def _reset_and_switch_solver(self, method):
        system = self.system
        system.reset_cache()
        system.empty_recycles()
        system.converge_method = method
        if self.hexanol:
            self._reset_evaporators()
        system.simulate()

    def _bugfix_barrage(self, error):
        attempts = []
        if self.hexanol:
            def evaporators_then_simulate():
                for unit in self.evaporators:
                    unit._V_first_effect = None
                self._reset_evaporators()
                self.system.simulate()
            attempts.append(evaporators_then_simulate)
        attempts += [self._reset_and_reload,
                     lambda: self._reset_and_switch_solver('fixedpoint'),
                     lambda: self._reset_and_switch_solver('aitken')]
        last = error
        for f in attempts:
            self._n_barrage = getattr(self, '_n_barrage', 0) + 1
            try:
                f()
                return
            except Exception as e:
                last = e
                str_e = str(e).lower()
                if 'sugar concentration' in str_e or 'permeate moisture' in str_e:
                    raise InfeasibleSpecification(str(e)[:200])
        raise last

    def reset(self):
        system = self.system
        system.converge_method = self._default_converge_method
        try:
            system.reset_cache()
            system.empty_recycles()
            if self.hexanol:
                self._reset_evaporators()
        except Exception:
            pass
        self._simulate_baseline()

    # %% Collection

    def _purity(self):
        AA = self.product_stream
        return AA.imass['AcrylicAcid'] / AA.F_mass

    def _MPSP(self):
        return self.tea.solve_price(self.product_stream) / self._purity()

    def _ferm_product_mass(self):
        broth = self.broth
        # 'CalciumLactate' is HP's surrogate for calcium 3-hydroxypropionate (2 HP-eq.)
        return broth.imass['HP'] + broth.imol['CalciumLactate'] * 2. * _MW_HP

    def collect(self):
        bst = self.bst
        AA = self.product_stream
        broth = self.broth
        R302, R303 = self.R302, self.R303
        F_AA = AA.imass['AcrylicAcid']
        if not F_AA > 0.:
            raise RuntimeError('no acrylic acid produced')
        purity = self._purity()
        # set the (impure) product price to the solved price, as HP's get_AA_MPSP does,
        # so TEA sales/product_stream_sales are evaluated at the MPSP
        AA.price = self.tea.solve_price(AA)
        MPSP = AA.price / purity
        F_fp = self._ferm_product_mass()
        sugar_feed = self.S302.ins[0]
        F_sugar = (sugar_feed.imass['Glucose'] + sugar_feed.imass['Xylose']
                   + sugar_feed.imass['Sucrose'] * 2 * _MW_glucose / _MW_sucrose)
        V_broth = broth.F_vol
        titer = F_fp / V_broth
        y_req, t_req, r_req = self._requested
        titer_err = (titer - t_req) / t_req
        if abs(titer_err) > self.titer_rtol:
            raise RuntimeError(f'titer {titer:.2f} g/L deviates from specification {t_req:.2f} g/L')
        units = self.system.units
        process_units = [i for i in units if i not in self.utility_facilities]
        heat_utilities = sum([i.heat_utilities for i in process_units], [])
        utils = bst.process_tools.utils
        power = [i.power_utility for i in units]
        lime = R302.ins[5]
        D = R302.design_results
        N = D.get('Number of reactors', float('nan'))
        V_reactor = D.get('Reactor volume', float('nan'))
        D3 = R303.design_results
        BT = self.BT
        natural_gas = getattr(BT, 'natural_gas', None)
        try:
            GWP = self.lca.GWP
            FEC = self.lca.FEC
        except Exception:
            GWP = FEC = float('nan')
        acid_cost = sum([s.cost for s in self._acid_feeds])
        out = dict(
            MPSP=MPSP,
            recovery=F_AA / F_fp,
            GWP=GWP, FEC=FEC,
            F_final_product=F_AA,
            F_ferm_product=F_fp,
            F_sugar_fed=F_sugar,
            F_feedstock=self.feedstock_stream.F_mass,
            feedstock_price=self.feedstock_stream.price,
            V_broth=V_broth,
            F_water_broth=broth.imass['Water'],
            V_fermenters=N * V_reactor,
            tau_fermentation=R302.tau,
            base_cost=lime.cost,
            acid_cost=acid_cost,
            heating_duty=utils.get_heating_duty(heat_utilities),
            cooling_duty=utils.get_cooling_duty(heat_utilities),
            electricity_consumption=utils.get_electricity_consumption(power),
            electricity_production=utils.get_electricity_production(power),
            natural_gas_flow=natural_gas.F_mass if natural_gas is not None else 0.,
            electricity_price=bst.PowerUtility.price,
            natural_gas_price=BT.natural_gas_price,
            # adapter-specific
            x_yield_measured=F_fp / F_sugar if F_sugar else float('nan'),
            x_titer_measured=titer,
            x_titer_rel_error=titer_err,
            x_AA_purity=purity,
            x_AA_price=AA.price,
            x_F_sugar_ferm_feed_conc=(sugar_feed.imass['Glucose', 'Xylose', 'Sucrose'].sum()
                                      / sugar_feed.F_vol),
            x_F_other_sugars_fed=sugar_feed.imass['Arabinose'] if 'Arabinose' in sugar_feed.chemicals else 0.,
            x_F301_V=self.F301.V,
            x_F301_water_evaporated=self.F301.outs[1].imass['Water'],
            x_dilution_water=self.M304.ins[1].F_mass,
            x_evaporator_water_removed=sum([i.outs[1].imass['Water'] for i in self.separation_evaporators]),
            x_cell_mass=broth.imass['FermMicrobe'],
            x_xylitol_broth=broth.imass['Xylitol'],
            x_succinic_acid_broth=broth.imass['SuccinicAcid'],
            x_residual_sugar_broth=broth.imass['Glucose', 'Xylose', 'Sucrose'].sum(),
            x_CO2_fermentation_vent=R302.outs[0].imass['CO2'] + R303.outs[0].imass['CO2'],
            x_air_fermentation=R302.ins[8].F_mass,
            x_CO2_fresh_fermentation=R302.ins[6].F_mass,
            x_lime_fermentation=lime.F_mass,
            x_gypsum=self._gypsum.F_mass if self._gypsum is not None else float('nan'),
            x_acid_acidulation=sum([s.imass['H2SO4'] for s in self._acid_feeds]),
            x_N_fermenters=N,
            x_V_fermenter=V_reactor,
            x_V_seed_train=D3.get('Reactor volume', float('nan')) * D3.get('Number of reactors', float('nan')),
            x_fermenters_installed_cost=R302.installed_cost,
            x_HXN_energy_balance_percent_error=getattr(self.HXN, 'energy_balance_percent_error', float('nan')),
            x_n_bugfix_attempts=getattr(self, '_n_barrage', 0),
            x_cell_mass_conversion=R302.glucose_to_biomass_rxn.X,
        )
        return out
