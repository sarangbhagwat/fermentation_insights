# -*- coding: utf-8 -*-
"""
Adapter for the oxalic acid biorefinery (``biorefineries.oxalic``,
Bioindustrial-Park; Lu et al., Nat. Commun. 2025, submitted): sugarcane
juicing -> aerobic batch fermentation by *I. orientalis* with partial lime
(Ca(OH)2) neutralization (0.91 mol Ca(OH)2 / mol oxalic acid, calcium oxalate
in broth) -> cell-mass centrifugation -> acidulation with H2SO4 -> gypsum
filtration. The final product is the clarified oxalic acid broth (crude
aqueous solution, ~2-5 wt% oxalic acid); MPSP is per kg *pure* oxalic acid.

Faithful to ``analyses/fermentation/TRY_analysis_sugarcane_oxalic_broth.py``:

* the model is loaded with the default parameter distributions
  (``parameter-distributions_sugarcane_oxalic-broth_100mL.xlsx``) and
  evaluated at their baseline, then production capacity is set at baseline
  TRY (feedstock flow adjusted to the desired annual production); during the
  TRY sweep the **feedstock capacity is fixed**;
* yield: ``spec.spec_1`` is the fraction of the theoretical yield
  (2 mol oxalic acid / mol glucose = 0.9995 g/g); titer: 5-100 g/L; titer is
  set by evaporation (F301, V <= 0.8) or dilution (M304); sugar concentration
  is capped at 800 g/L;
* MPSP = ``tea.solve_price(product) / purity``.

Only the sugarcane system exists in the module (``systems/sugarcane``); the
glucose/corn/cornstover parameter-distribution files have no system, and the
"glucose" in the published TRY results filenames is a hard-coded string, so no
other feedstock is offered. The published analyses use neutralization only.
"""
from __future__ import annotations

import os
import shutil
import tempfile
import warnings

import numpy as np

from ..harness import (BiorefineryAdapter, InfeasibleSpecification, Parameter,
                       describe_chemical, describe_feedstock)

__all__ = ('Adapter',)

_PARAM_FILE = ('analyses', 'full', 'parameter_distributions', 'oxalic_broth_product',
               'parameter-distributions_sugarcane_oxalic-broth_100mL.xlsx')

#: Parameters set by the harness (TRY) or meaningless (blank)
_SKIP = ('blank parameter', 'fermentation hp yield', 'fermentation hp titer',
         'fermentation hp productivity')

#: Canonical names for the module's parameters (lower-case; biosteam
#: normalizes capitalization of parameter names)
_CANONICAL = {
    'plant annual operating days': 'p_operating_days',
    'feedstock unit price': 'p_feedstock_price',
    'electricity unit price': 'p_electricity_price',
    'lime unit price': 'p_base_price',
    'sulfuric acid unit price': 'p_acid_price',
    'federal corporate tax rate': 'p_income_tax',
    'internal rate of return': 'p_IRR',
    'fermentation i. orientalis yield': 'p_cell_mass_yield',
    'boiler efficiency': 'p_boiler_efficiency',
    'turbogenerator efficiency': 'p_turbogenerator_efficiency',
}

#: UnitGroup names (system_sc_broth.unit_groups) -> canonical areas
_AREA_MAP = {
    'feedstock acquisition': 'feedstock',
    'feedstock juicing': 'conversion',
    'fermentation': 'fermentation',
    'separation': 'separation',
    'upgrading': 'upgrading',              # empty in this biorefinery
    'wastewater': 'wastewater',
    'storage & other facilities': 'facilities',
    'boiler & turbogenerator': 'chp',
    'cooling utility facilities': 'facilities',
    'heat exchanger network': 'facilities',
    # metric-only groups without units
    'natural gas (for steam generation)': None,
    'fixed operating cost': None,
    'electricity consumption': None,
    'heating duty': None,
    'excess electricity': None,
}

_MW_OA = 90.03488   # OxalicAcid (thermosteam)
_SUCROSE_TO_GLUCOSE = 2 * 180.15588 / 342.29648


class Adapter(BiorefineryAdapter):
    name = 'oxalic'
    area_map = _AREA_MAP

    def __init__(self):
        self._mod = None
        self._loaded_config = None

    # ------------------------------------------------------------------ configs
    def configurations(self):
        return [dict(
            config_id='oxalic_neutral_sugarcane',
            feedstock='sugarcane',
            feedstock_type='sugar_juice',
            fermentation_product='OxalicAcid',
            final_product='OxalicAcid',
            ferm_product_chemical='OxalicAcid',
            final_product_chemical='OxalicAcid',
            product_form='clarified_broth',   # crude aqueous solution, not crystals
            neutralization=True,
            neutralizing_agent='Ca(OH)2',
            neutralization_degree=0.91,       # mol Ca(OH)2 / mol oxalic acid (baseline)
            neutralized_salt='CalciumOxalate',
            n_acid_groups=2,
            pKa1=1.25, pKa2=4.14,
            final_product_stoich_mass_ratio=1.0,
            upgrading=False,
            separation_variant='',
            separation_principle='centrifugation+acidulation+gypsum_filtration',
            aerobic=True,
            fermentation_mode='batch',
            T_fermentation=303.15,
            product_water_solubility=143.,    # g/L at 25 C (oxalic acid)
            salt_water_solubility=0.0067,     # g/L calcium oxalate
            organism='I. orientalis',
        )]

    # --------------------------------------------------------------------- load
    def load(self, config):
        if self._mod is None:
            cwd = os.getcwd()
            tmp = tempfile.mkdtemp(prefix='oxalic_load_')
            os.chdir(tmp)  # the system module writes diagrams/figures to cwd
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    from biorefineries import oxalic
                    from biorefineries.oxalic.models.sugarcane import models_sc_broth as mod
            finally:
                os.chdir(cwd)
                shutil.rmtree(tmp, ignore_errors=True)
            self._mod = mod
            self._oxalic_path = os.path.dirname(oxalic.__file__)
        mod = self._mod
        import biosteam as bst
        self._bst = bst
        sysmod = __import__('biorefineries.oxalic.systems.sugarcane.system_sc_broth',
                            fromlist=['x'])
        self._sysmod = sysmod
        u, s = mod.u, mod.s
        self.system = mod.oxalic_sys
        self.tea = mod.oxalic_tea
        self.lca = mod.oxalic_lca
        self.spec = spec = mod.spec
        self.model = model = mod.model
        self.R302, self.R303, self.F301, self.M304 = u.R302, u.R303, u.F301, u.M304
        self.S302, self.BT = u.S302, u.BT701
        self.feedstock_stream = mod.feedstock
        self.product_stream = mod.product_stream
        self.broth = self.R302.outs[1]
        self.lime = s.lime_fresh
        self.acid = s.sulfuric_acid_fresh2
        self.gypsum = s.gypsum
        self.theoretical_yield = float(mod.theoretical_max_g_oxalic_per_g_glucose)

        # --- as in the TRY-analysis script
        spec.reactor.neutralization = bool(config['neutralization'])
        model.parameters = ()
        model.load_parameter_distributions(os.path.join(self._oxalic_path, *_PARAM_FILE))
        model.exception_hook = 'warn'
        model.metrics_at_baseline()
        sysmod.spec_set_production_capacity(spec.desired_annual_production)
        self._mpsp()

        # Unit groups for the harness: same units as the biorefinery's groups,
        # but material costs counted once (extend_feed_ends=False; the default
        # traces feeds through U101 and storage tanks, counting the feedstock in
        # both 'feedstock acquisition' and 'feedstock juicing' and stored
        # chemicals both in their consumer group and in storage).
        self.unit_groups = [
            bst.UnitGroup(g.name, g.units, extend_feed_ends=False,
                          filter_savings=(g.name != 'heat exchanger network'))
            for g in mod.unit_groups]

        self._baseline = dict(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
        self._F_feed_base = float(self.feedstock_stream.F_mass)
        self._biomass_base = float(self.R302.regular_biomass_conversion)
        self._converge_method = self.system.converge_method
        self._probe_max_titer()
        self._loaded_config = config

    def _mpsp(self):
        P = self.product_stream
        purity = P.imass['OxalicAcid'] / P.F_mass
        price = self.tea.solve_price(P)
        P.price = price
        return price / purity, purity

    def _probe_max_titer(self):
        """Maximum titer at a few yields: no dilution and evaporation (F301)
        up to the first of the two limits of ``load_titer_with_glucose`` --
        sugar concentration after evaporation = 800 g/L or F301.V = 0.8."""
        import flexsolve as flx
        sysmod, spec = self._sysmod, self.spec
        tis = spec.titer_inhibitor_specification
        Smax = tis.max_sugar_concentration
        ys, ts = [], []

        def dS(V):
            sysmod.F301_titer_obj_fn(V)
            return tis.calculate_sugar_concentration() - Smax
        try:
            for yf in (0.05, 0.08, 0.12, 0.17, 0.25, 0.4, 0.6, 0.9):
                spec.load_spec_1(yf)
                self.M304.water_to_sugar_mol_ratio = 0.
                self.R302.titer_to_load = 1.
                V = 0.8 if dS(0.8) <= 0. else flx.IQ_interpolation(dS, 0., 0.8, ytol=0.5, maxiter=100)
                dS(V)
                ys.append(yf * self.theoretical_yield)
                ts.append(float(self.R302.effluent_titer))
            self._max_titer_pts = (np.array(ys), np.array(ts))
        except Exception as e:  # fall back to the published TRY grid limit
            warnings.warn(f'max titer probe failed ({e}); using t_max ~ 520 y')
            self._max_titer_pts = None
        finally:
            self._reload_baseline()

    def _reload_baseline(self):
        spec = self.spec
        spec.load_specifications(**self._baseline)
        self.system.simulate()

    # ------------------------------------------------------------------- domain
    def domain(self):
        ty = self.theoretical_yield
        # Cell mass yield (regular_biomass_conversion, fraction of sugar) is an
        # uncertain parameter, U(0.12, 0.18), baseline 0.15. Above
        # y = (1 - biomass) * ty the cell-mass yield is decreased
        # ('sequential' policy). y_max uses the upper bound of the parameter so
        # that no sample crosses its inflection; the baseline inflection is
        # reported as inflection_yield and the per-sample one as
        # x_inflection_yield.
        biomass_hi = 0.18
        for p in self.model.parameters:
            if p.name.lower() == 'fermentation i. orientalis yield':
                biomass_hi = float(np.ravel(p.distribution.upper)[0])
        return dict(
            theoretical_yield=ty,
            inflection_yield=(1. - self._biomass_base) * ty,
            y_min=0.05 * ty,
            y_max=(1. - biomass_hi) * ty,
            t_min=5., t_max=100.,
            baseline_productivity=float(self._baseline['spec_3']),
            baseline_yield=float(self._baseline['spec_1']) * ty,
            baseline_titer=float(self._baseline['spec_2']),
        )

    def max_titer(self, y):
        pts = self._max_titer_pts
        if pts is None:
            return 520. * y / self.theoretical_yield
        ys, ts = pts
        # t_max(y) is slightly concave, so linear interpolation is conservative;
        # 0.995 keeps the grid's top point clear of the 800 g/L check.
        return float(np.interp(y, ys, ts) * 0.995)

    # -------------------------------------------------------------- descriptors
    def descriptors(self, config):
        # Same as BiorefineryAdapter.descriptors, which fails on biosteam 2.54
        # (`System` has no `thermo` attribute); chemicals from a stream instead.
        out = {}
        chems = self.product_stream.chemicals
        for key, prefix in (('ferm_product_chemical', 'fp_'),
                            ('final_product_chemical', 'pp_')):
            ID = config.get(key)
            if ID and ID in chems:
                out.update({prefix + k: v for k, v in describe_chemical(chems[ID]).items()})
        out.update(describe_feedstock(self.feedstock_stream))
        return out

    # --------------------------------------------------------------- parameters
    def parameters(self):
        out = []
        for p in self.model.parameters:
            key = p.name.lower()
            if key in _SKIP:
                continue
            if key == 'desired annual production':
                out.append(self._feedstock_capacity_parameter(p))
                continue
            out.append(Parameter.from_biosteam(p, name=_CANONICAL.get(key)))
        return out

    def _feedstock_capacity_parameter(self, p):
        # In the TRY analysis the feedstock capacity is fixed at the value that
        # gives the desired annual production at baseline TRY; production is
        # proportional to feedstock flow at fixed TRY, so the sampled desired
        # production maps linearly onto the wet feedstock flow.
        dist, x0, F0 = p.distribution, float(p.baseline), self._F_feed_base
        feed = self.feedstock_stream

        def setter(F):
            feed.F_mass = F
        return Parameter(
            name='p_feedstock_capacity', baseline=F0,
            inverse_cdf=lambda u, d=dist: float(d.inv(u)) / x0 * F0,
            setter=setter,
            description=(f'{p.element_name}: {p.name} [{p.units}] mapped to wet '
                         f'sugarcane flow [kg/h] = F_base*x/x_base '
                         f'(F_base={F0:.0f} kg/h, x_base={x0:.0f} t/y)'))

    # --------------------------------------------------------------- simulation
    def set_try(self, y, t, r):
        spec, sysmod = self.spec, self._sysmod
        yf = y / self.theoretical_yield
        if not (0. < yf < 1.):
            raise InfeasibleSpecification(f'yield {y:.3f} g/g outside (0, theoretical)')
        mt = self.max_titer(y)
        if t > mt * 1.03:  # clearly beyond the sugar-concentration/evaporation limit
            raise InfeasibleSpecification(f'titer {t:.1f} g/L > attainable {mt:.1f} g/L')
        try:
            spec.load_specifications(spec_1=yf, spec_2=t, spec_3=r)
        except Exception as e:
            msg = str(e).lower()
            if 'sugar concentration' in msg or type(e).__name__ == 'InfeasibleRegion':
                raise InfeasibleSpecification('sugar concentration > 800 g/L') from e
            if 'opposite sign' in msg:
                # IQ bracket failed: check whether the titer is unattainable at
                # the evaporation (V=0.8) or dilution limits.
                self.M304.water_to_sugar_mol_ratio = 0.
                if sysmod.F301_titer_obj_fn(0.8) < 0.:
                    raise InfeasibleSpecification('titer not attainable by evaporation (F301.V<=0.8)') from e
                sysmod.F301_titer_obj_fn(0.)
                if sysmod.M304_titer_obj_fn(40000.) > 0.:
                    raise InfeasibleSpecification('titer not attainable by dilution') from e
            raise

    def simulate(self):
        system = self.system
        try:
            system.simulate()
            return
        except Exception as e:
            if 'sugar concentration' in str(e).lower():
                raise InfeasibleSpecification('sugar concentration > 800 g/L') from e
            first = e
        # bugfix barrage (as in the TRY script, but keeping feedstock flow fixed)
        target = dict(spec_1=self.spec.spec_1, spec_2=self.spec.spec_2, spec_3=self.spec.spec_3)
        F = self.feedstock_stream.F_mass
        attempts = (None, 'fixedpoint', 'aitken')
        try:
            for method in attempts:
                try:
                    system.reset_cache()
                    system.empty_recycles()
                    if method is None:
                        self.spec.load_specifications(**self._baseline)
                        system.simulate()
                    else:
                        system.converge_method = method
                    self.feedstock_stream.F_mass = F
                    self.spec.load_specifications(**target)
                    system.simulate()
                    return
                except Exception as e:
                    last = e
            raise RuntimeError(f'bugfix barrage failed: {type(first).__name__}: {first}; '
                               f'last: {type(last).__name__}: {last}')
        finally:
            system.converge_method = self._converge_method
            self.feedstock_stream.F_mass = F

    def reset(self):
        F = self.feedstock_stream.F_mass
        system = self.system
        system.converge_method = self._converge_method
        system.reset_cache()
        system.empty_recycles()
        self.feedstock_stream.F_mass = F
        self._reload_baseline()

    # ------------------------------------------------------------------ collect
    def collect(self):
        bst = self._bst
        mpsp, purity = self._mpsp()
        P, broth, R302 = self.product_stream, self.broth, self.R302
        F_final = float(P.imass['OxalicAcid'])
        F_ferm = float((broth.imol['OxalicAcid'] + broth.imol['CalciumOxalate']) * _MW_OA)
        sugar = self.S302.ins[0]
        F_sugar = float(sugar.imass['Glucose'] + sugar.imass['Xylose']
                        + sugar.imass['Sucrose'] * _SUCROSE_TO_GLUCOSE)
        V_broth = float(broth.F_vol)
        D302 = R302.design_results
        D303 = self.R303.design_results
        excluded = (bst.BoilerTurbogenerator, bst.CoolingTower, bst.ChilledWaterPackage)
        process_units = [i for i in self.system.units if not isinstance(i, excluded)]
        hus = [hu for i in process_units for hu in i.heat_utilities]
        heating = sum(hu.duty for hu in hus if hu.flow * hu.duty > 0) / 1e6
        cooling = -sum(hu.duty for hu in hus if hu.flow * hu.duty < 0) / 1e6
        pu = self.system.power_utility
        BT = self.BT
        lca = self.lca
        tis = self.spec.titer_inhibitor_specification
        out = dict(
            MPSP=mpsp,
            recovery=F_final / F_ferm if F_ferm else float('nan'),
            GWP=float(lca.GWP), FEC=float(lca.FEC),
            F_final_product=F_final,
            F_ferm_product=F_ferm,
            F_sugar_fed=F_sugar,
            F_feedstock=float(self.feedstock_stream.F_mass),
            feedstock_price=float(self.feedstock_stream.price),
            V_broth=V_broth,
            F_water_broth=float(broth.imass['Water']),
            V_fermenters=float(D302.get('Reactor volume', np.nan) * D302.get('Number of reactors', np.nan)),
            tau_fermentation=float(R302.tau),
            base_cost=float(self.lime.cost),
            acid_cost=float(self.acid.cost),
            heating_duty=heating, cooling_duty=cooling,
            electricity_consumption=pu.consumption / 1e3,
            electricity_production=pu.production / 1e3,
            natural_gas_flow=float(BT.natural_gas.F_mass),
            electricity_price=float(bst.PowerUtility.price),
            natural_gas_price=float(BT.natural_gas_price),
            # extras
            x_purity=purity,
            x_MPSP_per_kg_broth=float(P.price),
            x_product_flow_total=float(P.F_mass),
            x_annual_production_pure=F_final * self.tea.operating_hours / 1e3,  # t/y
            x_yield_measured=F_ferm / F_sugar if F_sugar else float('nan'),
            x_titer_measured=F_ferm / V_broth if V_broth else float('nan'),
            x_titer_module=float(R302.effluent_titer),
            x_yield_fraction_theoretical=float(self.spec.spec_1),
            x_inflection_yield=(1. - float(R302.regular_biomass_conversion)) * self.theoretical_yield,
            x_cell_mass_yield_param=float(R302.regular_biomass_conversion),
            x_cell_mass=float(broth.imass['FermMicrobe']),
            x_acetic_acid_broth=float(broth.imass['AceticAcid']),
            x_calcium_oxalate_broth=float(broth.imass['CalciumOxalate']),
            x_evaporator_V=float(self.F301.V),
            x_evaporator_water_removed=float(self.F301.outs[1].F_mass),
            x_dilution_water=float(self.M304.ins[1].F_mass),
            x_sugar_conc_after_evaporation=float(tis.calculate_sugar_concentration()),
            x_lime_flow=float(self.lime.F_mass),
            x_H2SO4_flow=float(self.acid.imass['H2SO4']),
            x_gypsum_flow=float(self.gypsum.F_mass),
            x_gypsum_cost=float(self.gypsum.cost),
            x_CO2_fermentation_vent=float(R302.outs[0].imass['CO2']),
            x_air_flow=float(R302.ins[8].F_mass),
            x_fresh_CO2=float(R302.ins[6].F_mass),
            x_CSL_flow=float(R302.ins[2].F_mass),
            x_n_fermenters=float(D302.get('Number of reactors', np.nan)),
            x_V_fermenter_each=float(D302.get('Reactor volume', np.nan)),
            x_V_seed_train=float(D303.get('Reactor volume', np.nan) * D303.get('Number of reactors', np.nan)),
            x_wastewater_flow=float(self._sysmod.u.M501.outs[0].F_mass),
            x_HXN_qbal_error=float(getattr(self.spec.HXN, 'energy_balance_percent_error', np.nan)),
            x_GWP_no_ec=float(lca.GWP - lca.net_electricity_GWP),
            x_FEC_no_ec=float(lca.FEC - lca.net_electricity_FEC),
            x_BT_natural_gas_cost=float(BT.natural_gas_price * BT.natural_gas.F_mass),
        )
        return out
