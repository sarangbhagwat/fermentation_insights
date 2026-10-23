# -*- coding: utf-8 -*-
"""
Adapter for the corn -> isobutanol (+ ethanol) biorefinery
(``biorefineries.isobutanol``, Bioindustrial-Park; fed-batch S. cerevisiae
fermentation simulated with nskinetics).

TRY-imposed fermentation
------------------------
The biorefinery has no yield/titer/productivity-specified mode: its fermentor
V406 (nskinetics ``FermentationSaccharomycesEthanolIsobutanol``) integrates a
kinetic model, and the "TRY-informed" campaign is a kinetic-parameter
optimization (Bayesian optimization warm-started from TRY-objective studies),
not a TRY specification. This adapter therefore replaces V406's kinetic run
by a stoichiometric one that keeps every other process assumption:

* Reference: the module's ``opt_IRR`` scenario (best-IRR trial of the
  metabolic_split_12d profitability study; reproduces the module's pinned
  MPSPs). Its kinetic run is simulated once at load and calibrates the
  byproduct split per g glucose added: residual glucose, cell mass (Yeast),
  acetic acid, and the fraction of sugar to "other" sinks (cell mass +
  acetate + respiration), plus the O2 demand per kg cell mass.
* At yield ``y`` (g isobutanol / g glucose added; theoretical 0.4114), cell
  mass, acetate and residual glucose stay at their reference yields and the
  remaining sugar goes to ethanol (S. cerevisiae's competing fermentative
  sink, exactly what the kinetic model trades off). Ethanol vanishes at the
  inflection yield ``Y_IBO*(1 - f_residual - f_other)`` (~0.36 g/g); above it
  cell mass would have to fall (regime change), so ``y_max`` = inflection.
  CO2 closes the carbon balance and the vent/aeration are handled as in
  nskinetics (``_finalize_effluent``; air = O2 demand x safety factor).
* Titer ``t`` (g isobutanol produced / L broth leaving V406) is imposed by
  the module's own feeding machinery: total glucose loading S (g/L water) is
  solved so that the broth titer equals ``t``; S <= the reference initial
  charge concentration (175.9 g/L) -> batch at S (S301 split 1); above ->
  fed-batch with the initial charge at 175.9 g/L and 600 g/L spikes (split
  from the volume balance). The feed/spike concentrations are set with the
  module's concentration actuators (F301/M301, F302/M302). The module's
  fed-batch working-volume cap (final/initial <= 20) is enforced
  (InfeasibleSpecification).
* Productivity ``r`` sets the fermentation time ``tau = t/r``.

Products: isobutanol is the fermentation and final product; ethanol (fuel
grade), DDGS and corn oil are sold at their default prices as coproduct
credits. MPSP = the module's purity-adjusted isobutanol price solved at the
TEA's IRR (``system.solve_TEA_at_IRR``), USD/kg pure isobutanol.
``F_ferm_product`` is isobutanol leaving V406 in broth + vent (the vent
isobutanol is captured by the scrubber and recycled, as in the module's
recovery metric), so recovery = isobutanol in product / F_ferm_product.
"""
from __future__ import annotations

import ctypes
import glob
import math
import os
import sys

import numpy as np

from ..harness import BiorefineryAdapter, InfeasibleSpecification, Parameter

__all__ = ('Adapter',)

MW_GLUCOSE = 180.156
MW_IBO = 74.122
MW_ETOH = 46.069
Y_IBO_THEO = MW_IBO / MW_GLUCOSE           # 0.4114 g/g (C6H12O6 -> C4H10O + 2CO2 + H2O)
Y_ETOH_THEO = 2 * MW_ETOH / MW_GLUCOSE     # 0.5114 g/g
RHO_ETOH, RHO_IBO = 780., 795.             # kg/m3 liquid at ~32 C (broth-volume estimate only)

#: nskinetics commit matching BIP master's isobutanol kinetics (the
#: requirements.txt pin d8591e5 lacks k_17); install with pip --no-deps.
NSKINETICS_COMMIT = '1ea45768cc94d0bc0d4f02c0f2966ab5a40ed08a'

REFERENCE_SCENARIO = 'opt_IRR'
#: Values of the opt_IRR kinetic run measured with biosteam 2.54.0 /
#: thermosteam 0.54.2 / nskinetics 1ea4576 (used only to warn on drift).
REFERENCE_EXPECTED = dict(y=0.23082, f_residual=0.00819, f_other=0.11500)
#: Module MPSP pin of the reference scenario (scenarios.SCENARIOS).
REFERENCE_MPSP_PIN = 1.1635401445494427


def _preload_libpython():
    """libroadrunner's pip wheel links libpython3.X.so.1.0 but conda's python
    is static and the env's lib/ is not on its runpath; preload it globally
    so `import tellurium` works without LD_LIBRARY_PATH."""
    try:
        import roadrunner  # noqa: F401
        return
    except ImportError:
        pass
    v = f'{sys.version_info.major}.{sys.version_info.minor}'
    for path in glob.glob(os.path.join(sys.prefix, 'lib', f'libpython{v}.so*')):
        try:
            ctypes.CDLL(path, mode=ctypes.RTLD_GLOBAL)
            return
        except OSError:
            continue


class Adapter(BiorefineryAdapter):
    name = 'isobutanol'
    area_map = {
        'feedstock acquisition': 'feedstock',
        'feedstock saccharification': 'conversion',
        'sugar solution preparation': 'fermentation',  # S301 split, F301/F302 evaporators, dilution
        'fermentation': 'fermentation',
        'alcohol recovery': 'separation',
        'ethanol purification': 'separation',
        'isobutanol purification': 'separation',
        'storage and handling': 'facilities',
        'DDGS recovery': 'separation',        # stillage solids -> DDGS (coproduct recovery)
        'wastewater treatment': 'wastewater',
        'heat exchanger network': 'facilities',
        'boiler': 'chp',
        'turbogenerator': 'chp',
        'cooling utility facilities': 'facilities',
        'other facilities': 'facilities',
        # empty cost pseudo-groups (metrics only; no units)
        'natural gas (for steam generation)': None,
        'natural gas (for product drying)': None,
        'fixed operating cost': None,
        'excess electricity': None,
    }

    # TRY domain
    T_MIN = 5.          # g/L broth
    T_MAX = 80.         # g/L; module's highest kinetic isobutanol titer is 78 g/L
    Y_MIN = 0.05        # g/g (12 % of theoretical)

    def __init__(self):
        self._loaded = False
        self._spec = None
        self._S_guess = None

    # ------------------------------------------------------------------ configs
    def configurations(self):
        return [dict(
            config_id='IBO_EtOH_corn',
            feedstock='corn',
            feedstock_type='starch',
            fermentation_product='Isobutanol',
            final_product='Isobutanol',
            ferm_product_chemical='Isobutanol',
            final_product_chemical='Isobutanol',
            neutralization=False,
            neutralizing_agent='',
            n_acid_groups=0,
            pKa1=float('nan'),
            final_product_stoich_mass_ratio=1.0,
            upgrading=False,
            separation_variant='',
            separation_principle='heteroazeotropic distillation+decantation',
            aerobic=True,                 # aerated growth stage, then anaerobic production
            fermentation_mode='fed-batch',
            T_fermentation=305.15,
            product_water_solubility=85.,  # g/L at 25 C (8.5 wt%)
            coproducts='Ethanol,DDGS,crude_oil',
            reference_kinetic_scenario=REFERENCE_SCENARIO,
        )]

    # --------------------------------------------------------------------- load
    def load(self, config):
        if self._loaded:
            self.reset()
            return
        _preload_libpython()
        import biosteam as bst
        from biorefineries import isobutanol
        isobutanol.load(separation_processes=('IBO_EtOH',))
        from biorefineries.isobutanol import scenarios, system as S
        self.bst = bst
        self.S = S
        f = self.f = S.f
        self.V406 = V406 = f.V406
        self.system = S.corn_EtOH_IBO_sys
        self.tea = S.corn_EtOH_IBO_sys_tea
        self.unit_groups = list(S.unit_groups)
        self.feedstock_stream = f.corn
        self.product_stream = f.isobutanol
        self.ethanol_stream = f.ethanol
        self.fbs = S.fbs_spec
        # reference kinetic run (sets the scenario's kinetics, feeding and
        # enzyme burden, simulates it with the module's load_simulate)
        bundle = scenarios.load_scenario(REFERENCE_SCENARIO)
        self._model = bundle['model']
        self._reference_MPSP = S.solve_TEA_at_IRR(
            stream_IDs=('isobutanol',))['purity_adjusted_MPSPs']['isobutanol']
        self._calibrate()
        # switch to TRY mode: no kinetics, no enzyme burden
        S.set_active_burden(None)
        self._kinetic_run = V406._run
        V406._run = self._run_try
        self.set_try(self.ref['y'], self.ref['t'], self.ref['r'])
        self.simulate()
        self._baseline_MPSP = self._solve_MPSP()
        print(f'[isobutanol] reference kinetic MPSP {self._reference_MPSP:.5f} '
              f'(module pin {REFERENCE_MPSP_PIN:.5f}); TRY-imposed baseline MPSP '
              f'{self._baseline_MPSP:.5f} USD/kg', flush=True)
        self._snapshot = S.snapshot_flowsheet_state()
        self._snapshot_S = self._S_guess
        self._loaded = True

    def _calibrate(self):
        V406 = self.V406
        vent, eff = V406.outs
        feeds = [i for i in V406.ins if i is not V406.ins[3]]
        G = sum(i.imass['Glucose'] for i in feeds)  # corn mash: all sugar as glucose
        prod = lambda ID: eff.imass[ID] + vent.imass[ID]
        y = prod('Isobutanol') / G
        e = prod('Ethanol') / G
        f_res = eff.imass['Glucose'] / G
        yeast_in = sum(i.imass['Yeast'] for i in feeds)
        ref = dict(
            G=G, y=y, e=e, f_residual=f_res,
            x=(eff.imass['Yeast'] - yeast_in) / G,          # g cell mass / g glucose
            acetate=prod('AceticAcid') / G,
            f_other=1. - y / Y_IBO_THEO - e / Y_ETOH_THEO - f_res,
            O2_per_yeast=V406.stoich_O2_flow_req / eff.imass['Yeast'],  # kmol O2 / kg
            t=prod('Isobutanol') / eff.F_vol,                # g/L broth
            tau=V406.tau,
            c0=self.fbs.target_conc, cs=self.fbs.spike_conc,
            vol_cap=self.S.parameters['max_fed_batch_volume_ratio'],
        )
        W = sum(i.ivol['Water'] for i in feeds)
        ref['S'] = G / W                                       # g glucose / L water fed
        ref['k_solids'] = (eff.F_vol - eff.ivol['Water']) / G   # m3 non-water broth / kg glucose
        # ... of which not ethanol/isobutanol (solids, residual sugar, yeast)
        ref['k_nonalcohol'] = ref['k_solids'] - e / RHO_ETOH - y / RHO_IBO
        ref['r'] = ref['t'] / ref['tau']
        ref['y_inflection'] = Y_IBO_THEO * (1. - ref['f_residual'] - ref['f_other'])
        R = ref['vol_cap']
        ref['S_max'] = (ref['c0'] + (R - 1.) * ref['cs']) / R
        self.ref = ref
        self._S_guess = ref['S']
        for k, v in REFERENCE_EXPECTED.items():
            if abs(ref[k] - v) > 0.02 * max(abs(v), 0.05):
                print(f'[isobutanol] WARNING reference {k} = {ref[k]:.5f} '
                      f'(expected {v:.5f})', flush=True)

    # ------------------------------------------------------------------- domain
    def domain(self):
        ref = self.ref
        return dict(theoretical_yield=Y_IBO_THEO,
                    inflection_yield=ref['y_inflection'],
                    y_min=self.Y_MIN, y_max=ref['y_inflection'],
                    t_min=self.T_MIN, t_max=self.T_MAX,
                    baseline_productivity=ref['r'],
                    baseline_yield=ref['y'], baseline_titer=ref['t'])

    def max_titer(self, y):
        """Broth titer reached at the largest glucose loading the module's
        fed-batch envelope allows (initial charge 175.9 g/L, 600 g/L spikes,
        final/initial volume <= 20), with the non-water broth volume per g
        glucose (reference corn solids + yield-dependent alcohols); 1 %
        below the bound so the top grid point stays feasible."""
        # t [kg/m3] = y / (1/S [m3/kg water per kg glucose] + k(y) [m3 non-water / kg glucose])
        return 0.99 * y / (1. / self.ref['S_max'] + self._k_nonwater(y))

    def _ethanol_yield(self, y):
        ref = self.ref
        return max(0., Y_ETOH_THEO * (1. - y / Y_IBO_THEO - ref['f_residual'] - ref['f_other']))

    def _k_nonwater(self, y):
        """Non-water broth volume per kg glucose fed at yield y [m3/kg]."""
        return self.ref['k_nonalcohol'] + self._ethanol_yield(y) / RHO_ETOH + y / RHO_IBO

    # --------------------------------------------------------------- parameters
    _CANONICAL = {
        'Plant annual operating days': 'p_operating_days',
        'Feedstock unit price': 'p_feedstock_price',
        'Electricity unit price': 'p_electricity_price',
        'Federal corporate tax rate': 'p_income_tax',
        'Feedstock capacity': 'p_feedstock_capacity',
    }
    _EXCLUDE = ('Blank parameter',
                'Isobutanol selling unit price (sale)')  # solved for (MPSP)

    def parameters(self):
        out = []
        for p in self._model.parameters:
            if str(p.element_name).startswith('Kinetics') or p.name in self._EXCLUDE:
                continue
            out.append(Parameter.from_biosteam(p, name=self._CANONICAL.get(p.name)))
        return out

    # --------------------------------------------------------------- TRY specs
    def set_try(self, y, t, r):
        ref = self.ref
        if not (0. < y <= ref['y_inflection'] * (1 + 1e-9)):
            raise InfeasibleSpecification(
                f'yield {y:.4g} outside (0, inflection {ref["y_inflection"]:.4g}] g/g')
        inv_S = y / t - self._k_nonwater(y)   # 1/S needed [m3/kg]
        if inv_S * ref['S_max'] < 1. / 1.03:  # loading beyond S_max (3 % slack)
            raise InfeasibleSpecification(
                f'titer {t:.4g} g/L at yield {y:.4g} needs more glucose than the '
                f'fed-batch envelope allows')
        self._spec = dict(y=y, t=t, r=r, e=self._ethanol_yield(y))
        self.V406.tau = t / r

    def _run_try(self):
        """Stoichiometric replacement of V406's kinetic _run."""
        V406 = self.V406
        spec, ref = self._spec, self.ref
        vent, eff = V406.outs
        air = V406.ins[3]
        vent.empty()
        eff.empty()
        eff.mix_from([i for i in V406.ins if i is not air])
        for rxn in V406.pre_reactions:  # sucrose hydrolysis (no sucrose in corn mash)
            rxn.force_reaction(eff)
        G = eff.imass['Glucose']
        eff.imass['Glucose'] = ref['f_residual'] * G
        eff.imass['Isobutanol'] += spec['y'] * G
        eff.imass['Ethanol'] += spec['e'] * G
        eff.imass['Yeast'] += ref['x'] * G
        eff.imass['AceticAcid'] += ref['acetate'] * G
        V406._finalize_effluent(eff, vent, None)   # NH3 -> 0, CO2 by C balance, vent
        # aeration (nskinetics AerationSpec.set_air_stream, O2 demand per cell mass)
        air.empty()
        stoich_O2 = ref['O2_per_yeast'] * eff.imass['Yeast']
        V406.stoich_O2_flow_req = stoich_O2
        V406.stoich_air_flow_req = stoich_O2 / 0.21
        air.imol['O2'] = 0.21
        air.imol['N2'] = 0.79
        air.F_mol = stoich_O2 / 0.21 * V406.aeration.safety_factor
        vent.imol['O2'] += max(0., air.imol['O2'] - stoich_O2)
        vent.imol['N2'] += air.imol['N2']
        if V406.converge_air_supply:
            V406._converge_air_supply()

    def _ferm_titer(self):
        vent, eff = self.V406.outs
        return (eff.imass['Isobutanol'] + vent.imass['Isobutanol']) / eff.F_vol

    def _set_feeding(self, S):
        """Glucose loading S (g/L water, initial + spike feeds)."""
        ref, fbs, V406 = self.ref, self.fbs, self.V406
        c0, cs = ref['c0'], ref['cs']
        if S <= c0:
            sigma, c_init = 1., S
        else:
            sigma = (1. / S - 1. / cs) / (1. / c0 - 1. / cs)
            c_init = c0
            ratio = 1. + (1. - sigma) * c0 / (sigma * cs)
            if sigma <= 0. or ratio > ref['vol_cap']:
                raise InfeasibleSpecification(
                    f'fed-batch volume ratio {ratio:.3g} > {ref["vol_cap"]:g} '
                    f'(glucose loading {S:.4g} g/L)')
        fbs.splitter.split = sigma
        try:
            fbs._estimate_set_actuators(
                desired_conc=c_init, concentrator=fbs.feed_concentrator,
                diluter=fbs.feed_diluter,
                f_simulate_units_sequential=fbs._simulate_feed_units,
                terminal_stream=V406.ins[0])
            if sigma < 1.:
                fbs._estimate_set_actuators(
                    desired_conc=cs, concentrator=fbs.spike_concentrator,
                    diluter=fbs.spike_diluter,
                    f_simulate_units_sequential=fbs._simulate_spike_units,
                    terminal_stream=V406.ins[2])
            else:
                fbs._simulate_spike_units()
        except Exception as e:
            if type(e).__name__ == 'FeedingStrategyError':
                raise InfeasibleSpecification(f'feed concentration: {e}'[:250])
            raise

    def _impose_titer(self, rtol=2e-4, maxiter=12):
        t = self._spec['t']
        S = self._S_guess or self.ref['S']
        pts = []
        for i in range(maxiter):
            self._set_feeding(S)
            self.V406.simulate()
            tm = self._ferm_titer()
            pts.append((S, tm))
            if abs(tm / t - 1.) < rtol:
                break
            if len(pts) >= 2 and abs(pts[-1][1] - pts[-2][1]) > 1e-9:
                (S1, t1), (S2, t2) = pts[-2], pts[-1]
                S_new = S2 + (t - t2) * (S2 - S1) / (t2 - t1)
            else:
                S_new = S * t / tm
            S = float(np.clip(S_new, 0.3 * S, 3. * S))
            if S > 1.05 * self.ref['S_max']:
                raise InfeasibleSpecification(
                    f'glucose loading {S:.4g} g/L > {self.ref["S_max"]:.4g} g/L')
        else:
            raise RuntimeError(f'titer not reached ({tm:.4g} vs {t:.4g} g/L)')
        self._S_guess = S
        return S

    def _drift_state(self):
        return np.array([*self.S._simulation_drift_state()[:5],
                         self.product_stream.F_mass, self._ferm_titer()])

    def _load_simulate(self, n_sims=6, sim_rtol=1e-4, max_final_drift=1e-2):
        """TRY-mode analog of system.load_simulate: [impose titer at the
        fermentor; simulate system] until the tracked state stops moving."""
        prev = self._drift_state()
        floors = np.full(prev.shape, 1e-12)
        drift = np.inf
        for k in range(n_sims):
            self._impose_titer()
            self.system.simulate()
            curr = self._drift_state()
            drift = float(np.max(np.abs(curr - prev) / np.maximum(np.abs(prev), floors)))
            prev = curr
            if drift <= sim_rtol and k > 0:
                break
        if drift > max_final_drift:
            raise RuntimeError(f'TRY load_simulate not converged (drift {drift:.3g})')
        if abs(self._ferm_titer() / self._spec['t'] - 1.) > 5e-3:
            raise RuntimeError(f'titer drifted: {self._ferm_titer():.4g} vs {self._spec["t"]:.4g}')
        self._last_drift = drift

    def simulate(self):
        sys_ = self.system
        try:
            self._load_simulate()
            return
        except InfeasibleSpecification:
            raise
        except Exception as e:
            first = e
        # bugfix barrage (system.run_bugfix_barrage with the TRY loader)
        method = sys_.converge_method
        for solver in (None, 'fixedpoint', 'aitken'):
            try:
                sys_.reset_cache()
                sys_.empty_recycles()
                if solver:
                    sys_.converge_method = solver
                self._S_guess = self.ref['S']
                sys_.simulate()
                self._load_simulate()
                return
            except InfeasibleSpecification:
                raise
            except Exception:
                continue
            finally:
                sys_.converge_method = method
        raise first

    def reset(self):
        if getattr(self, '_snapshot', None) is not None:
            self.S.restore_flowsheet_state(self._snapshot)
            self._S_guess = self._snapshot_S

    def versions(self):
        v = super().versions()
        try:
            _preload_libpython()  # versions() runs before load()
            import nskinetics
            v['nskinetics'] = getattr(nskinetics, '__version__', None)
        except Exception:
            v['nskinetics'] = None
        v['nskinetics_commit_expected'] = NSKINETICS_COMMIT
        return v

    # ------------------------------------------------------------------ collect
    def _solve_MPSP(self):
        sol = self.S.solve_TEA_at_IRR(stream_IDs=('isobutanol',))
        self._last_IRR = sol['IRR']
        return sol['purity_adjusted_MPSPs']['isobutanol']

    def collect(self):
        f, V406, sys_, bst = self.f, self.V406, self.system, self.bst
        vent, eff = V406.outs
        MPSP = self._solve_MPSP()
        ibo, etoh = self.product_stream, self.ethanol_stream
        feeds = [i for i in V406.ins if i is not V406.ins[3]]
        G = sum(i.imass['Glucose'] for i in feeds)
        F_fp = eff.imass['Isobutanol'] + vent.imass['Isobutanol']
        F_final = ibo.imass['Isobutanol']
        BT = f.BT801
        units = list(sys_.units)
        facilities = set(sys_.facilities)
        heat = cool = 0.
        for u in units:  # process units + heat-exchanger network (net of integration)
            if u in facilities and not isinstance(u, bst.HeatExchangerNetwork):
                continue
            for hu in u.heat_utilities:  # biosteam UnitGroup sign convention
                if hu.flow * hu.duty > 0:
                    heat += hu.duty
                elif hu.flow * hu.duty < 0:
                    cool -= hu.duty
        el_cons = sum(u.power_utility.consumption for u in units) / 1e3
        el_prod = sum(u.power_utility.production for u in units) / 1e3
        ng_boiler = BT.natural_gas.F_mass if getattr(BT, 'natural_gas', None) is not None else 0.
        ng_other = 0.
        for u in units:
            if u is BT:
                continue
            ng = getattr(u, 'natural_gas', None)
            if ng is not None and hasattr(ng, 'F_mass'):
                ng_other += ng.F_mass
        design = V406.design_results
        hours = self.tea.operating_hours
        spec = self._spec
        out = dict(
            MPSP=MPSP,
            recovery=F_final / F_fp if F_fp else float('nan'),
            GWP=float('nan'), FEC=float('nan'),
            F_final_product=F_final,
            F_ferm_product=F_fp,
            F_sugar_fed=G,
            F_feedstock=f.corn.F_mass,
            feedstock_price=f.corn.price,
            V_broth=eff.F_vol,
            F_water_broth=eff.imass['Water'],
            V_fermenters=design['Reactor volume'] * design['Number of reactors'],
            tau_fermentation=V406.tau,
            base_cost=0., acid_cost=0.,
            heating_duty=heat / 1e6, cooling_duty=cool / 1e6,
            electricity_consumption=el_cons,
            electricity_production=el_prod,
            natural_gas_flow=ng_boiler + ng_other,
            electricity_price=bst.PowerUtility.price,
            natural_gas_price=getattr(BT, 'natural_gas_price', float('nan')),
            # adapter-specific
            x_IRR_default_prices=self._last_IRR,
            x_MPSP_reference_kinetic=self._reference_MPSP,
            x_purity_final_product=ibo.imass['Isobutanol'] / ibo.F_mass if ibo.F_mass else float('nan'),
            x_ethanol_yield=spec['e'],
            x_ethanol_titer=(eff.imass['Ethanol'] + vent.imass['Ethanol']) / eff.F_vol,
            x_F_ethanol_product=etoh.imass['Ethanol'],
            x_ethanol_price=etoh.price,
            x_ethanol_sales=etoh.price * etoh.F_mass * hours,
            x_F_DDGS=f.DDGS.F_mass, x_DDGS_price=f.DDGS.price,
            x_F_crude_oil=f.crude_oil.F_mass,
            x_cell_mass=eff.imass['Yeast'],
            x_CO2_vent=vent.imass['CO2'],
            x_F_ferm_product_vent=vent.imass['Isobutanol'],
            x_titer_broth_post_vent=eff.imass['Isobutanol'] / eff.F_vol,
            x_titer_water_basis=F_fp / eff.ivol['Water'],
            x_glucose_loading=self._S_guess,
            x_split_to_initial_feed=float(np.atleast_1d(self.fbs.splitter.split)[0]),
            x_F_air=V406.ins[3].F_mass,
            x_n_fermenters=design['Number of reactors'],
            x_natural_gas_boiler=ng_boiler,
            x_natural_gas_other=ng_other,
            x_ammonia_cost=f.ammonia.cost, x_lime_cost=f.lime.cost,
            x_sulfuric_acid_cost=f.sulfuric_acid.cost,
            x_sim_drift=self._last_drift,
        )
        return out
