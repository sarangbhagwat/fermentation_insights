# -*- coding: utf-8 -*-
"""
Adapter for the succinic acid biorefineries of Bioindustrial-Park
(``biorefineries.succinic``; Tran, Mishra, Bhagwat et al.), ported to
BioSTEAM 2.54 on the ``fi-succinic`` branch.

Configurations (as in the fermentation_insights publication):
feedstock in {glucose, corn, sugarcane, cornstover} x
{'succinic' (low-pH fermentation with partial pH control; no neutralization),
 'succinic_neutral' (full neutralization with lime)}.

Notes
-----
* Each feedstock lives in its own module-level system
  (``biorefineries.succinic.models_<feedstock>``) that also changes global
  BioSTEAM settings at import, so only ONE feedstock can be loaded per
  process. Switching between 'succinic' and 'succinic_neutral' of the same
  feedstock within a process is supported (``spec.neutralization``). The
  launcher runs every configuration in its own process (``--configs``).
* Uncertain parameters come from the current default (non-FGI) distribution
  files ``analyses/parameter_distributions/parameter-distributions_pilot-scale_batch_<feedstock>.xlsx``
  (commit "make non-FGI param distribs default"); their baselines are loaded
  with ``model.metrics_at_baseline()`` exactly as in the TRY scripts.
* Yield convention of the biorefinery: fraction of the theoretical maximum
  (2*MW_SA/MW_glucose = 1.311 g/g, same mass yield from xylose); the harness
  passes absolute g/g which are divided by 1.311 in ``set_try``.
* Neutral configurations (adapter_version 3): ``R302.base_neutralizes_product
  = True``, so the lime (1 mol Ca(OH)2/mol SA + 0.5 mol/mol acetic acid)
  forms calcium succinate (and acetate) in the broth, which the acidulation
  reactor R401 converts back to succinic acid with H2SO4 (93 wt%,
  ``sulfuric_acid_acidulation``, priced) before the gypsum filter S405
  removes the CaSO4 (stream ``gypsum``, priced at the module's gypsum price,
  0 at baseline: neither sold nor charged). Adapter v2 set only
  ``spec.neutralization`` (as the publication's TRY script), which with the
  module default ``base_neutralizes_product=False`` discards the purchased
  lime and bypasses R401 / R401_P / S405 (DESIGN.md defect D15). The low-pH
  configurations keep the module default (``base_neutralizes_product=False``:
  the partial pH-control lime is bought, but no salt forms and nothing is
  acidulated; adapter_version 2).
* MPSP = tea.solve_price(product) / purity (USD/kg pure succinic acid) and
  recovery = SA in product / SA(acid-eq.) in fermentation broth, as in
  ``analyses/TRY_analysis_*.py``. Points whose heat exchanger network energy
  balance error exceeds ``spec.tolerable_HXN_energy_balance_percent_error``
  are failed (the TRY scripts return NaN for them).
"""
from __future__ import annotations

import importlib
import os
import re

import numpy as np

from ..harness import BiorefineryAdapter, InfeasibleSpecification, Parameter

__all__ = ('Adapter',)


class _TiterBracketError(RuntimeError):
    """The biorefinery's load_titer could not bracket the titer specification."""

FEEDSTOCKS = ('glucose', 'corn', 'sugarcane', 'cornstover')
FEEDSTOCK_TYPES = {'glucose': 'pure_sugar', 'corn': 'starch',
                   'sugarcane': 'sugar_juice', 'cornstover': 'lignocellulose'}

#: Approximate maximum attainable "effective sugar concentration" [g sugar-eq / L broth]:
#: max feasible titer ~= y [g/g] * S_MAX. Estimated from the feasible-region boundary of
#: the publication's 50x50 TRY grids (t_max/y = 1268-1289, 651-672, 606-632, 381-405 g/L)
#: plus ~1 grid step of margin.
S_MAX = {'glucose': 1300., 'corn': 680., 'sugarcane': 640., 'cornstover': 410.}

#: Publication TRY domain (fractions of theoretical yield; g/L)
YIELD_FRACTION_MIN = 0.038
T_MIN, T_MAX = 10., 120.

#: Parameters set by the harness (excluded) and dummy parameters.
_EXCLUDED_PARAMETERS = ('blank parameter', 'succinic acid yield', 'succinic acid titer',
                        'succinic acid productivity')

_CANONICAL = {
    'plant annual operating days': 'p_operating_days',
    'feedstock unit price': 'p_feedstock_price',
    'natural gas unit price': 'p_natural_gas_price',
    'electricity unit price': 'p_electricity_price',
    'lime unit price': 'p_base_price',
    'feedstock capacity': 'p_feedstock_capacity',
    'a. succinogenes yield': 'p_cell_mass_yield',
    'boiler efficiency': 'p_boiler_efficiency',
}

_SUGAR_MASS_FACTORS = {'Glucose': 1., 'Xylose': 1.,  # SA mass yield from xylose = from glucose (1.311 g/g)
                       'Sucrose': 2 * 180.156 / 342.297}  # glucose-eq.

_UTILITY_PRODUCERS = ('BoilerTurbogenerator', 'CoolingTower', 'ChilledWaterPackage')
_FACILITY_TYPES = ('HeatExchangerNetwork', 'ChilledWaterPackage', 'CoolingTower',
                   'ProcessWaterCenter', 'FireWaterTank', 'CIPpackage',
                   'AirDistributionPackage')

_AREA_BY_NUMBER = {0: 'facilities', 100: 'feedstock', 200: 'conversion', 300: 'fermentation',
                   400: 'separation', 500: 'wastewater', 600: 'facilities', 700: 'chp',
                   800: 'facilities', 900: 'facilities', 1000: 'facilities'}


def _unit_area_number(unit):
    m = re.search(r'(\d+)', str(unit.ID))
    if not m:
        return None
    n = int(m.group(1))
    return (n // 100) * 100 if n >= 100 else 0


class Adapter(BiorefineryAdapter):
    name = 'succinic'
    area_map = {a: a for a in ('feedstock', 'conversion', 'fermentation', 'separation',
                               'upgrading', 'wastewater', 'chp', 'facilities')}

    def __init__(self):
        self._feedstock_loaded = None
        self.config = None

    # %% Configurations
    def configurations(self):
        configs = []
        for neutral in (False, True):
            for feed in FEEDSTOCKS:
                configs.append(dict(
                    config_id=f"succinic{'_neutral' if neutral else ''}_{feed}",
                    feedstock=feed,
                    feedstock_type=FEEDSTOCK_TYPES[feed],
                    fermentation_product='SuccinicAcid',
                    final_product='SuccinicAcid',
                    ferm_product_chemical='SuccinicAcid',
                    final_product_chemical='SuccinicAcid',
                    # 3 (neutral only): lime neutralizes the product, acidulation + gypsum
                    # filtration run (D15); 2: neutralizing_agent of the low-pH case (M1)
                    adapter_version=3 if neutral else 2,
                    neutralization=neutral,
                    # base dosed in fermentation (base_cost = its cost); the low-pH case
                    # still doses lime for partial pH control (0.42 mol-OH-eq/mol SA
                    # baseline) but does not convert the product (no acidulation
                    # downstream), hence neutralization=False there (v1: agent '')
                    neutralizing_agent='Ca(OH)2',
                    partial_pH_control_base='' if neutral else 'Ca(OH)2',
                    n_acid_groups=2,
                    pKa1=4.21, pKa2=5.64,
                    final_product_stoich_mass_ratio=1.0,
                    upgrading=False,
                    separation_variant='',
                    separation_principle=(
                        ('acidulation(H2SO4)+gypsum_filtration+' if neutral else '')
                        + 'evaporation+crystallization'),
                    aerobic=True,  # model sparges air (3.57 mol/mol SA) and fixes CO2
                    co2_fixation=True,
                    fermentation_mode='batch',
                    T_fermentation=303.15,
                    product_water_solubility=83.,  # g/L at 25 C (literature)
                ))
        return configs

    # %% Loading
    def load(self, config):
        feed = config['feedstock']
        if self._feedstock_loaded is not None and self._feedstock_loaded != feed:
            raise RuntimeError(
                'succinic adapter: only one feedstock can be loaded per process '
                f'(loaded {self._feedstock_loaded!r}, requested {feed!r}); run each '
                'configuration in its own process with --configs')
        if self._feedstock_loaded is None:
            self._import_models(feed)
        import biosteam as bst
        self.config = config
        spec = self.spec
        neutral = bool(config['neutralization'])
        spec.neutralization = neutral
        # neutral: the lime forms calcium succinate, acidulated by R401 (H2SO4) and
        # gypsum-filtered by S405 (their bypass specifications read this flag);
        # low pH: module default (partial pH-control lime, no salt, R401/S405 bypassed)
        self.R302.base_neutralizes_product = neutral
        # baseline parameter values from the default distribution file, as in the TRY scripts
        self.model.metrics_at_baseline()
        self.system.converge_method = self._converge_method
        self.unit_groups = self._make_unit_groups(bst)
        self.T_fermentation = float(self.R302.T)

    def _import_models(self, feed):
        import biosteam as bst
        models = importlib.import_module(f'biorefineries.succinic.models_{feed}')
        self.models = models
        self.model = model = models.model
        self.spec = models.spec
        self.tea = models.succinic_tea
        self.lca = models.succinic_LCA
        self.system = model._system
        self.product_stream = models.product_stream
        self.feedstock_stream = models.feedstock
        self.R302 = models.R302
        self.R303 = models.R303
        self.BT = models.BT
        s, u = models.s, models.u
        self.base_stream = s.base_fermentation
        self.acid_stream = s.sulfuric_acid_acidulation
        self.gypsum_stream = s.gypsum
        self.CO2_purchased_stream = s.CO2_fermentation
        self.natural_gas_drying = s.natural_gas_drying
        self.M304 = u.M304
        self.F301 = u.F301
        self.theoretical_yield = float(models.theoretical_max_g_succinic_acid_per_g_glucose)
        self._pre_fermenter_units = models.pre_fermenter_units_path
        self._converge_method = self.system.converge_method
        path = os.path.join(os.path.dirname(models.__file__), 'analyses',
                            'parameter_distributions',
                            f'parameter-distributions_pilot-scale_batch_{feed}.xlsx')
        self.parameter_distributions_file = path
        model.parameters = ()
        model.load_parameter_distributions(path, models.namespace_dict)
        model.exception_hook = 'raise'
        self._feedstock_loaded = feed
        self._utility_producers = [i for i in self.system.units
                                   if type(i).__name__ in _UTILITY_PRODUCERS]
        self._process_units = [i for i in self.system.units if i not in self._utility_producers]

    def _make_unit_groups(self, bst):
        # units upstream of the sugar-concentration evaporator F301 (pretreatment, saccharification,
        # juicing, liquefaction; also the feed flash of the glucose system)
        upstream_of_sugar_prep = set(self.F301.get_upstream_units())
        assign = {}
        for unit in self.system.units:
            tname = type(unit).__name__
            if tname == 'BoilerTurbogenerator':
                area = 'chp'
            elif tname in _FACILITY_TYPES:
                area = 'facilities'
            else:
                n = _unit_area_number(unit)
                area = _AREA_BY_NUMBER.get(n, 'facilities')
                # saccharification/liquefaction units carrying 300-series IDs (corn, corn stover)
                if area == 'fermentation' and unit in upstream_of_sugar_prep:
                    area = 'conversion'
            assign.setdefault(area, []).append(unit)
        return [bst.UnitGroup(area, units) for area, units in assign.items()]

    # %% Domain
    def domain(self):
        spec = self.spec
        theo = self.theoretical_yield
        inflection = (1. - spec.baseline_cell_mass_yield - spec.baseline_pyruvic_acid_yield) * theo
        return dict(
            theoretical_yield=theo,
            inflection_yield=inflection,
            y_min=YIELD_FRACTION_MIN * theo,
            y_max=inflection,
            t_min=T_MIN, t_max=T_MAX,
            baseline_productivity=float(spec.baseline_productivity),
            baseline_yield=float(spec.baseline_yield) * theo,
            baseline_titer=float(spec.baseline_titer),
        )

    def max_titer(self, y):
        return float(y) * S_MAX[self.config['feedstock']]

    # %% Parameters
    def parameters(self):
        out = []
        for p in self.model.parameters:
            key = p.name.lower().strip()
            if key in _EXCLUDED_PARAMETERS:
                continue
            prm = Parameter.from_biosteam(p, name=_CANONICAL.get(key))
            if key == 'a. succinogenes yield':
                prm.setter = self._cell_mass_yield_setter(p.setter)
                prm.description += (' (adapter: also sets spec.baseline_cell_mass_yield so that '
                                    'spec.load_yield does not reset it)')
            out.append(prm)
        return out

    def _cell_mass_yield_setter(self, setter):
        spec = self.spec
        def set_cell_mass_yield(x):
            # spec.load_yield resets R302 cell mass conversions to spec.baseline_cell_mass_yield,
            # which would make this parameter ineffective; keep both consistent.
            spec.baseline_cell_mass_yield = x
            setter(x)
        return set_cell_mass_yield

    # %% Simulation
    def set_try(self, y, t, r):
        yf = float(y) / self.theoretical_yield
        if not 0. < yf < 1.:
            raise InfeasibleSpecification(f'yield {y:.4g} g/g outside (0, theoretical)')
        self._target = (yf, float(t), float(r))

    def _load_specs(self, yf, t, r):
        """spec.load_specifications, translating the biorefinery's infeasibility signals:
        InfeasibleRegion('sugar concentration') -> InfeasibleSpecification (as in the TRY
        scripts, no retry); failure to bracket the titer specification in the biorefinery's
        load_titer function (titer unattainable by evaporation/dilution within its bounds)
        -> _TiterBracketError (confirmed after a reset in the bugfix barrage)."""
        try:
            self.spec.load_specifications(spec_1=yf, spec_2=t, spec_3=r)
        except Exception as e:
            msg = str(e).lower()
            if 'sugar concentration' in msg:
                raise InfeasibleSpecification(f'sugar concentration infeasible: {e}')
            if 'opposite signs' in msg:
                raise _TiterBracketError(
                    f'titer {t:.3g} g/L not attainable (titer bracketing failed; closest '
                    f'attainable {float(self.R302.effluent_titer):.3g} g/L)') from e
            raise

    def _load_and_simulate(self):
        yf, t, r = self._target
        for i in self._pre_fermenter_units:
            i.simulate()
        self._load_specs(yf, t, r)
        self.system.simulate()

    def simulate(self):
        # The TRY scripts simulate each point twice: once in spec.evaluate_across_specs and
        # once more in the MPSP metric (model.specification); do the same.
        for _ in range(2):
            try:
                self._load_and_simulate()
            except InfeasibleSpecification:
                raise
            except Exception:
                self._bugfix_barrage()
        self._check_point()

    def _bugfix_barrage(self):
        """As in biorefineries.succinic._process_specification.evaluate_across_specs:
        reset & reload from baseline, then fixedpoint, then aitken."""
        system, spec = self.system, self.spec
        yf, t, r = self._target
        try:
            try:
                system.reset_cache()
                system.empty_recycles()
                self._load_specs(spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity)
                system.simulate()
                self._load_specs(yf, t, r)
                system.simulate()
            except _TiterBracketError as e:
                # unattainable titer confirmed from a fresh baseline state
                raise InfeasibleSpecification(str(e))
            except InfeasibleSpecification:
                raise
            except Exception:
                try:
                    self._reset_and_switch_solver('fixedpoint')
                except Exception:
                    self._reset_and_switch_solver('aitken')
        except _TiterBracketError as e:
            raise InfeasibleSpecification(str(e))
        finally:
            system.converge_method = self._converge_method

    def _reset_and_switch_solver(self, solver_ID):
        system, spec = self.system, self.spec
        yf, t, r = self._target
        system.reset_cache()
        system.empty_recycles()
        system.converge_method = solver_ID
        self._load_specs(yf, t, r)
        system.simulate()

    def _check_point(self):
        spec = self.spec
        yf, t, r = self._target
        titer = float(self.R302.effluent_titer)
        if not np.isfinite(titer) or abs(titer - t) > 0.02 * t:
            raise InfeasibleSpecification(f'titer target {t:.3g} g/L not attainable (got {titer:.3g})')
        err = float(spec.HXN.energy_balance_percent_error)
        if abs(err) >= spec.tolerable_HXN_energy_balance_percent_error:
            raise RuntimeError(f'HXN energy balance error {err:.3g}% >= '
                               f'{spec.tolerable_HXN_energy_balance_percent_error}% (TRY scripts discard such points)')

    def reset(self):
        system, spec = self.system, self.spec
        system.converge_method = self._converge_method
        system.reset_cache()
        system.empty_recycles()
        spec.load_specifications(spec_1=spec.baseline_yield, spec_2=spec.baseline_titer,
                                 spec_3=spec.baseline_productivity)
        system.simulate()

    # %% Collection
    def _sugar_fed(self):
        tot = 0.
        for s in (self.R302.ins[0], self.R303.ins[0]):
            for ID, f in _SUGAR_MASS_FACTORS.items():
                tot += f * s.imass[ID]
        return tot

    def collect(self):
        import biosteam as bst
        from biosteam.process_tools import utils as ptu
        from biorefineries.succinic.units import compute_succinic_acid_mass
        product, tea, spec, R302 = self.product_stream, self.tea, self.spec, self.R302
        broth = R302.outs[0]
        F_SA = product.imass['SuccinicAcid']
        purity = F_SA / product.F_mass
        product.price = price = tea.solve_price(product)
        MPSP = price / purity
        F_ferm = compute_succinic_acid_mass(broth)
        F_sugar = self._sugar_fed()
        heat_utilities = [hu for u in self._process_units for hu in u.heat_utilities]
        pu = self.system.power_utility
        NG_BT = self.BT.natural_gas.F_mass
        NG_dry = self.natural_gas_drying.F_mass
        tau = float(R302.tau)
        yf, t, r = self._target
        out = dict(
            MPSP=MPSP,
            recovery=F_SA / F_ferm if F_ferm else float('nan'),
            F_final_product=F_SA,
            F_ferm_product=F_ferm,
            F_sugar_fed=F_sugar,
            F_feedstock=self.feedstock_stream.F_mass,
            feedstock_price=self.feedstock_stream.price,
            V_broth=broth.F_vol,
            F_water_broth=broth.imass['Water'],
            V_fermenters=broth.F_vol * (R302.tau_batch_turnaround + tau) / R302.V_wf,
            tau_fermentation=tau,
            base_cost=self.base_stream.cost,
            acid_cost=self.acid_stream.cost,
            heating_duty=ptu.get_heating_duty(heat_utilities),
            cooling_duty=ptu.get_cooling_duty(heat_utilities),
            electricity_consumption=pu.consumption / 1e3,
            electricity_production=pu.production / 1e3,
            natural_gas_flow=NG_BT + NG_dry,
            electricity_price=bst.PowerUtility.price,
            natural_gas_price=self.BT.natural_gas_price,
        )
        try:
            out['GWP'] = float(self.lca.GWP)
            out['FEC'] = float(self.lca.FEC)
        except Exception:
            pass
        theo = self.theoretical_yield
        x = dict(
            x_yield_fraction_of_theoretical=yf,
            x_yield_measured=F_ferm / F_sugar if F_sugar else float('nan'),
            x_titer_achieved=float(R302.effluent_titer),
            x_inflection_yield=(1. - spec.baseline_cell_mass_yield - spec.baseline_pyruvic_acid_yield) * theo,
            x_cell_mass_yield=float(spec.baseline_cell_mass_yield),
            x_product_purity=purity,
            x_product_stream_price=price,
            x_cell_mass=broth.imass['FermMicrobe'],
            x_base_flow=self.base_stream.F_mass,
            x_acid_flow=self.acid_stream.F_mass,
            # adapter v3: salt formed in the broth and gypsum filtered off (0 in the low-pH case)
            x_calcium_succinate_broth=broth.imass['CalciumSuccinate'],
            x_gypsum_flow=self.gypsum_stream.F_mass,
            x_CO2_purchased=self.CO2_purchased_stream.F_mass,
            x_CO2_recycled=R302.ins[5].F_mass,
            x_air_fermentation=R302.ins[8].F_mass,
            x_natural_gas_BT=NG_BT,
            x_natural_gas_drying=NG_dry,
            x_sugar_evaporator_V=float(getattr(self.F301, 'V', float('nan'))),
            x_sugar_evaporator_water_removed=self.F301.outs[1].F_mass,
            x_dilution_water=self.M304.ins[1].F_mass,
            x_sugar_conc_fermenter_feed=sum(f * R302.ins[0].imass[ID] for ID, f in _SUGAR_MASS_FACTORS.items()) / R302.ins[0].F_vol,
            x_HXN_energy_balance_error=float(spec.HXN.energy_balance_percent_error),
            x_operating_days=float(tea.operating_days),
            x_BT_steam_heating_duty=ptu.get_heating_duty(self.BT.heat_utilities),
        )
        out.update(x)
        return out
