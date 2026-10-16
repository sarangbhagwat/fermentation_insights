# -*- coding: utf-8 -*-
"""
Ethanol adapter: Bioindustrial-Park ethanol biorefineries on BioSTEAM 2.54.

Configurations
--------------
=====================  ==========================================================
config_id              module / design
=====================  ==========================================================
ethanol_cornstover     ``biorefineries.cellulosic.Biorefinery('corn stover ethanol')``
                       (= ``biorefineries.cornstover``; NREL 2011 dilute-acid
                       pretreatment, integrated saccharification and
                       co-fermentation by Z. mobilis, distillation + mol. sieves,
                       WWT, CHP).
ethanol_sugarcane      ``biorefineries.cane.Biorefinery('S1')`` after
                       ``YRCP2025()`` (latest assumptions; juicing, juice
                       evaporation, batch fermentation with yeast recycle,
                       distillation + mol. sieves, high-rate WWT, bagasse CHP).
ethanol_corn           ``biorefineries.corn.Biorefinery('conventional dry-grind')``
                       (liquefaction, SSF, distillation + mol. sieves, DDGS).
=====================  ==========================================================

Bioindustrial-Park changes required on BioSTEAM 2.54 (branch ``fi-ethanol``):
(i) ``cellulosic/units.py``: ``size_batch`` was called with positional
arguments in the wrong order (N_reactors passed as V_wf and V_wf as V_max),
which sized the cornstover SSCF vessels as ~3700 x 0.9 m3 tanks (fermentor
capital ~0); (ii) ``corn/units.py``: the SSF unit inherits a 1-inlet size from
``bst.BatchBioreactor`` but is given mash + yeast; ``_N_ins = 2`` was added.

How titer, rate and yield (TRY) are imposed
-------------------------------------------
Definitions are identical for all configurations:

* fermentable sugar fed ``S`` [kg/h, glucose-eq. mass] = sugars entering the
  fermenters that the module's organism can convert to ethanol:
  cornstover -- glucose + xylose after the in-vessel saccharification of the
  SSCF feed (sugars partially routed through the seed train are counted once);
  sugarcane -- glucose + sucrose x 2 MW_glc/MW_suc in the fermentor feed;
  corn -- glucose in the SSF feed (starch fully liquefied upstream).
* ``F_ferm_product`` = ethanol produced = ethanol leaving the fermenters
  (liquid broth + vent) minus ethanol entering them; measured yield
  ``F_ferm_product/F_sugar_fed`` equals ``y`` exactly.
* titer ``t`` = ethanol in the liquid broth leaving the main fermentor / broth
  volumetric flow (the module's own titer definition). Because ~1-2 % of the
  ethanol leaves with the CO2 vent (recovered by the vent scrubber),
  ``F_ferm_product/V_broth`` is 1-2 % above ``t``;
  ``x_F_ethanol_broth/V_broth`` equals ``t``.
* productivity: fermentation time tau = t / r (cornstover: co-fermentation
  time of the SSCF vessels; the 60 h saccharification time is kept).

Yield: theoretical yield 2 MW_EtOH / MW_glc = 0.511 g/g (3 xylose -> 5 EtOH
has the same mass yield). Byproduct/cell-mass reactions are kept at the
module's values.

* cornstover: all ethanol conversions (glucose & xylose, seed train &
  co-fermentation) are set to ``phi * X_max`` where ``X_max = 1 - sum of the
  competing (cell mass, glycerol, succinic, xylitol) conversions of the same
  sugar in the same reactor``; the 3 % contamination loss is kept. ``phi`` is
  solved exactly (quadratic, because the seed train consumes sugar) from the
  actual saccharified feed at run time by a specification on R303. The
  inflection yield is the yield at ``phi = 1``.
* sugarcane / corn: ethanol conversion ``X = y / 0.511``. In these modules cell
  growth only consumes the sugar left after ethanol formation (sequential
  reactions), so competing reactions never limit the ethanol conversion and
  the inflection yield equals the theoretical yield.

Domain: y_min = max(0.05 x theoretical, lowest yield at which 1.15 t_min is
attainable); y_max = 0.995 x inflection (cornstover; the inflection shifts
slightly with the parameter-dependent glucose/xylose ratio) or 0.99 x
inflection (sugarcane, corn: X = 1 means zero cell growth, a degenerate state in
which the sugarcane WWT sludge solver fails); t_min = 10 g/L; t_max =
min(150 g/L, max_titer(y_max)); ``max_titer(y)`` = 0.98 x the analytical titer
at the sugar-concentration limit below.

Titer: by the sugar concentration entering fermentation, adjusting water the
way each module does.

* cornstover: the saccharification dilution water (M301) is solved (Newton
  iterations on the whole system) for the target titer instead of the
  module's 20 % solids-loading spec. Infeasible if the water would have to be
  negative: the pretreated slurry is already ~22 % solids, so the maximum
  sugar concentration (~127 g/L broth) is only ~7 % above the module's and
  t_max ~ 58 g/L (module baseline 53 g/L).
* sugarcane: the module's own titer specification (juice multi-effect
  evaporation or dilution water), with an outer loop correcting its target
  (it is solved before the yeast recycle converges; error up to ~3 %).
  Infeasible when the fermentor feed would need > 350 g/L fermentable sugar
  (very-high-gravity limit; module baseline ~300 g/L; the module's own
  evaporation limit corresponds to > 1000 g/L) or the module cannot evaporate
  enough.
* corn: the total water entering the slurry tank (V307) is solved for the
  target titer (the recycled process water is set as its complement to the
  backwater recycle and other water inputs) instead of the module's fixed
  slurry solids content (0.311). Infeasible if the slurry total solids would
  exceed 35 wt% (very-high-gravity limit; module baseline 26 wt%) or no
  process water would remain.

Other choices: corn uses its TEA's design IRR (15 %, ``corn.tea``) instead of
the break-even IRR the module solves at load. Parameters: sugarcane -- the
cane model's parameters that apply to S1; cornstover -- as in
``cornstover/model.py`` (BioSTEAM default +-10 % prices, utility prices,
feedstock flow rate, income tax; saccharification conversion and time) plus
the boiler natural-gas price; corn -- BioSTEAM defaults (the module has no
model). Cost-item parameters (base cost, exponent, kW of every cost item) of
the BioSTEAM defaults are excluded. Area groups: one UnitGroup per canonical
area (TEA units plus zero-cost units, so every feed is attributed); corn uses
Lang-factor installed costs as its TEA does.
"""
from __future__ import annotations

import math
import warnings


from ..harness import BiorefineryAdapter, InfeasibleSpecification, Parameter
from ..schema import CANONICAL_AREAS

__all__ = ('Adapter',)

nan = float('nan')
KG_PER_GAL_ETHANOL = 3.78541 * 0.789

_COMMON = dict(
    fermentation_product='Ethanol',
    final_product='Ethanol',
    ferm_product_chemical='Ethanol',
    final_product_chemical='Ethanol',
    neutralization=False,
    neutralizing_agent='',
    n_acid_groups=0,
    pKa1=nan,
    final_product_stoich_mass_ratio=1.0,
    upgrading=False,
    separation_variant='',
    separation_principle='distillation+molecular_sieves',
    aerobic=False,
    product_water_solubility=1000.,
)

CONFIGURATIONS = (
    dict(config_id='ethanol_cornstover', feedstock='cornstover',
         feedstock_type='lignocellulose', fermentation_mode='batch',
         T_fermentation=305.15, organism='Z. mobilis (co-fermenting)',
         cell_recycle=False,
         module='biorefineries.cellulosic:corn stover ethanol', **_COMMON),
    dict(config_id='ethanol_sugarcane', feedstock='sugarcane',
         feedstock_type='sugar_juice', fermentation_mode='batch',
         T_fermentation=305.15, organism='S. cerevisiae',
         cell_recycle=True,
         module='biorefineries.cane:S1 (YRCP2025)', **_COMMON),
    dict(config_id='ethanol_corn', feedstock='corn',
         feedstock_type='starch', fermentation_mode='batch',
         T_fermentation=305.15, organism='S. cerevisiae',
         cell_recycle=False,
         module='biorefineries.corn:conventional dry-grind', **_COMMON),
)


# %% Helpers

def _ids(units):
    return {u.ID: u for u in units}


def _numeric_area(ID):
    digits = ''.join(c for c in ID if c.isdigit())
    return int(digits[0]) if digits else None


def _make_unit_groups(tea, system, assign, group_class=None):
    """One UnitGroup per canonical area holding the TEA units plus the
    zero-cost system units (mixers, splitters, WWT internals) so that every
    system feed is attributed to an area."""
    import biosteam as bst
    group_class = group_class or bst.UnitGroup
    groups = {a: [] for a in CANONICAL_AREAS}
    unassigned = []
    tea_units = set(tea.units)
    units = list(tea.units) + [u for u in system.units
                               if u not in tea_units and not u.purchase_cost]
    for u in units:
        a = assign(u)
        if a is None:
            unassigned.append(u.ID)
            a = 'facilities'
        groups[a].append(u)
    return [group_class(a, groups[a]) for a in CANONICAL_AREAS], unassigned


def _lang_factor_group_class(tea):
    """UnitGroup whose installed cost is ``tea.lang_factor`` x purchase cost,
    consistent with TEAs that use a Lang factor (corn)."""
    import biosteam as bst

    class LangFactorUnitGroup(bst.UnitGroup):
        __slots__ = ()

        def get_installed_cost(self):
            return tea.lang_factor * self.get_purchase_cost()
    return LangFactorUnitGroup


def _heat_duties(units):
    """Signed heating and cooling duty [GJ/h] of process units (utility
    producing facilities excluded; heat-exchanger-network savings included)."""
    import biosteam as bst
    exclude = (bst.BoilerTurbogenerator, bst.CoolingTower, bst.ChilledWaterPackage)
    heating_IDs = {a.ID for a in bst.HeatUtility.heating_agents}
    heating = cooling = 0.
    for u in units:
        if isinstance(u, exclude):
            continue
        for hu in u.heat_utilities:
            agent = hu.agent
            if agent is None:
                continue
            if agent.ID in heating_IDs:
                heating += hu.duty
            else:
                cooling -= hu.duty
    return heating / 1e6, cooling / 1e6


def _predict_water(impl, y, t, W):
    """Initial guess of the dilution water [kg/h] for target yield/titer from
    the last simulated state (broth volume is ~linear in water)."""
    try:
        b = impl.broth
        S = impl._sugar_fed() if hasattr(impl, '_sugar_fed') else sum(impl._sugars().values())
        produced = impl.ethanol_produced()
        if produced <= 0 or b.F_vol <= 0:
            return W
        E_broth = y * S * b.imass['Ethanol'] / produced
        rho = b.chemicals.Water.rho('l', b.T, b.P)
        return max(W + (E_broth / t - b.F_vol) * rho, 0.)
    except Exception:
        return W


def _pct_parameter(p, name, description=''):
    """Wrap a biosteam parameter defined in % as a fraction."""
    d = p.distribution
    return Parameter(
        name=name, baseline=float(p.baseline) / 100.,
        inverse_cdf=lambda u, d=d: float(d.inv(u)) / 100.,
        setter=lambda v, f=p.setter: f(100. * v),
        description=f'{p.element_name}: {p.name} [{p.units}] (stored as fraction) {description}'.strip(),
    )


def _default_parameters(system, feedstock, exclude_streams=()):
    """BioSTEAM's default uncertainty (``Model.load_default_parameters``,
    +-10 % triangular) for stream prices, electricity price, heat-utility
    prices, feedstock flow rate and income tax; cost-item parameters (base
    cost, exponent, kW of every cost item) are excluded."""
    import biosteam as bst
    from biosteam.evaluation.evaluation_tools import parameter as ep
    from biosteam.evaluation.evaluation_tools import triang
    model = bst.Model(system)
    ep.add_all_stream_price_params(model, triang, feedstock, True)
    ep.add_power_utility_price_param(model, triang)
    ep.add_heat_utility_price_params(model, triang)
    ep.add_feedstock_flow_rate_param(model, feedstock, triang)
    ep.add_basic_TEA_params(model, triang, False)
    params = [p for p in model.parameters
              if not any(p.element is i for i in exclude_streams)]
    return model, params


def _param_key(p):
    el = p.element
    el = getattr(el, 'ID', el)
    return (str(el), str(p.name))


# %% Configuration implementations

class _Impl:
    """Common machinery; subclasses implement module-specific parts."""
    theoretical_yield = None
    t_min = 10.
    t_max_cap = 150.

    def __init__(self, adapter):
        self.adapter = adapter
        self.infeasible_message = None

    # -- robust simulation --
    #: Whether emptying recycle streams is a valid last resort (not for corn,
    #: whose heat-integrated beer loop then has zero flow).
    allow_empty_recycles = True

    def snapshot_recycles(self):
        """Store the converged recycle streams (used to restart failed
        simulations from a good state)."""
        self._recycle_snapshot = [(i, i.copy()) for i in self.system.get_all_recycles()]

    def restore_recycles(self):
        for i, c in getattr(self, '_recycle_snapshot', ()):
            i.copy_like(c)

    def simulate_system(self):
        """Simulate with the usual fallbacks: reset cache, restart from the
        converged baseline recycles, and (if allowed) empty recycles."""
        sys = self.system
        steps = [None, 'reset_cache', 'restore']
        if self.allow_empty_recycles:
            steps.append('empty')
        last = None
        for step in steps:
            try:
                if step == 'reset_cache':
                    sys.reset_cache()
                elif step == 'restore':
                    self.recover()
                elif step == 'empty':
                    sys.empty_recycles()
                    sys.reset_cache()
                sys.simulate()
                return
            except InfeasibleSpecification:
                raise
            except Exception as e:
                msg = str(e)
                if 'evaporate' in msg:  # module's titer limit (sugarcane)
                    raise InfeasibleSpecification(f'sugar concentration unattainable: {msg[:150]}')
                last = e
        raise last

    def recover(self):
        """Bring the system to a state from which the next simulation can
        converge."""
        self.restore_recycles()
        self.system.reset_cache()

    def max_titer(self, y):
        """Approximate maximum titer [g/L] at yield y (sugar-concentration limit)."""
        return y * self.S_max

    def check_infeasible(self):
        if self.infeasible_message:
            msg = self.infeasible_message
            self.infeasible_message = None
            raise InfeasibleSpecification(msg)

    # -- common measurements --
    def ethanol_produced(self):
        return sum(s.imass['Ethanol'] for s in self.fermentor_outs) - \
               sum(s.imass['Ethanol'] for s in self.fermentor_ins)

    def titer(self):
        b = self.broth
        return b.imass['Ethanol'] / b.F_vol

    def mpsp(self):
        eth = self.ethanol
        price = self.tea.solve_price(eth)
        eth.price = price
        purity = eth.imass['Ethanol'] / eth.F_mass
        return price / purity, price, purity


class _Cornstover(_Impl):
    name = 'cornstover'

    def load(self):
        from biorefineries.cellulosic import Biorefinery
        br = Biorefinery('corn stover ethanol')
        self.br = br
        self.system = br.system
        self.tea = br.tea
        self.ethanol = br.ethanol
        self.feedstock = br.cornstover
        self.R302, self.R303, self.M301 = br.R302, br.R303, br.M301
        chems = self.R303.chemicals
        MW_E, MW_G, MW_X = chems.Ethanol.MW, chems.Glucose.MW, chems.Xylose.MW
        self.theoretical_yield = 2 * MW_E / MW_G
        self.m = {'Glucose': 2 * MW_E / MW_G, 'Xylose': 5 * MW_E / (3 * MW_X)}
        self.simulate_system()
        self.snapshot_recycles()
        # ---- module baseline (before any adapter change) ----
        self.module_baseline = self._snapshot()
        # ---- maximum (inflection) ethanol conversions ----
        seed = self.R302.reactions.X
        cof = self.R303.cofermentation.X
        loss = self.R303.loss.X
        self.a = {'Glucose': 1 - seed[1:4].sum(), 'Xylose': 1 - seed[5:9].sum()}
        self.b = {'Glucose': seed[1:4].sum(), 'Xylose': seed[5:9].sum()}
        self.c = {'Glucose': 1 - cof[1:4].sum(), 'Xylose': 1 - cof[5:9].sum()}
        self.L = {'Glucose': loss[0], 'Xylose': loss[1]}
        # ---- install TRY specifications ----
        from thermosteam.network import ProcessSpecification
        self.y_target = None
        self.t_target = None
        self.W_over_feed = self.M301.ins[2].imass['Water'] / self.feedstock.F_mass
        specs = self.M301._specifications
        assert specs[1].f.__name__ == 'update_moisture_content'
        specs[1] = ProcessSpecification(self._set_saccharification_water, (), ())
        self.R303._specifications.insert(0, ProcessSpecification(self._set_ethanol_conversions, (), ()))
        self.fermentor_ins = (self.R303.ins[0],)
        self.fermentor_outs = (self.R303.outs[0], self.R303.outs[1], self.R302.outs[0])
        self.broth = self.R303.outs[1]
        # baseline TRY (module values)
        mb = self.module_baseline
        self.baseline = dict(y=mb['y'], t=mb['t'], r=mb['t'] / self.R303.tau_cofermentation)
        S = mb['S']
        S_tot = sum(S.values())
        self.inflection_yield = sum(self.m[k] * S[k] * self._f(k, 1.) for k in S) / S_tot
        # maximum sugar concentration (no saccharification water)
        V_min = mb['V_broth'] - mb['W_sacc'] / 1000.
        self.S_max = S_tot / V_min  # kg/m3 = g/L
        self.unit_groups, self.unassigned = _make_unit_groups(self.tea, self.system, self._assign)

    def _snapshot(self):
        S = self._sugars()
        E = self.ethanol_produced_raw()
        b = self.R303.outs[1]
        mpsp, price, purity = self.mpsp_raw()
        return dict(S=S, y=E / sum(S.values()), t=b.imass['Ethanol'] / b.F_vol,
                    V_broth=b.F_vol, W_sacc=self.M301.ins[2].imass['Water'],
                    MPSP=mpsp, price=price, purity=purity,
                    MESP_gal=price * KG_PER_GAL_ETHANOL)

    def ethanol_produced_raw(self):
        R303, R302 = self.R303, self.R302
        return (R303.outs[0].imass['Ethanol'] + R303.outs[1].imass['Ethanol']
                + R302.outs[0].imass['Ethanol'] - R303.ins[0].imass['Ethanol'])

    def mpsp_raw(self):
        eth = self.ethanol
        price = self.tea.solve_price(eth)
        purity = eth.imass['Ethanol'] / eth.F_mass
        return price / purity, price, purity

    def _sugars(self):
        tmp = self.R303.ins[0].copy()
        self.R303.saccharification(tmp)
        return {'Glucose': float(tmp.imass['Glucose']), 'Xylose': float(tmp.imass['Xylose'])}

    def _f(self, k, phi):
        """Overall ethanol conversion of sugar k at fraction phi of the
        maximum conversions (seed train on a side draw + co-fermentation)."""
        q = self.R303.saccharification_split
        a, b, c, L = self.a[k], self.b[k], self.c[k], self.L[k]
        return q * phi * a + (1 - q * phi * a - q * b) * (1 - L) * phi * c

    def _set_ethanol_conversions(self):
        y = self.y_target
        if y is None:
            return
        S = self._sugars()
        S_tot = sum(S.values())
        q = self.R303.saccharification_split
        A = B = 0.
        for k, Sk in S.items():
            m, a, b, c, L = self.m[k], self.a[k], self.b[k], self.c[k], self.L[k]
            A += m * Sk * (q * a + (1 - q * b) * (1 - L) * c)
            B += m * Sk * q * a * (1 - L) * c
        C = y * S_tot
        disc = A * A - 4 * B * C
        if disc < 0:
            phi = 2.
        else:
            phi = (A - math.sqrt(disc)) / (2 * B)
        if phi > 1. + 1e-9:
            self.infeasible_message = (f'yield {y:.4f} above inflection yield '
                                       f'{sum(self.m[k]*S[k]*self._f(k, 1.) for k in S)/S_tot:.4f}')
            phi = 1.
        self.phi = phi
        seed = self.R302.reactions.X
        cof = self.R303.cofermentation.X
        seed[0] = phi * self.a['Glucose']
        seed[4] = phi * self.a['Xylose']
        cof[0] = phi * self.c['Glucose']
        cof[4] = phi * self.c['Xylose']

    def _set_saccharification_water(self):
        self.M301.ins[2].imass['Water'] = max(self.W_over_feed * self.feedstock.F_mass, 0.)

    def set_try(self, y, t, r):
        if y > self.theoretical_yield:
            raise InfeasibleSpecification('yield above theoretical')
        if t > 1.1 * y * self.S_max:
            raise InfeasibleSpecification(f'titer {t:.1f} g/L above max attainable ~{y*self.S_max:.1f} g/L')
        self.y_target = y
        self.t_target = t
        self.R303.tau_cofermentation = t / r
        F = self.feedstock.F_mass
        self.W_over_feed = _predict_water(self, y, t, self.W_over_feed * F) / F

    def simulate(self):
        t = self.t_target
        F = self.feedstock.F_mass
        for it in range(10):
            self.simulate_system()
            self.check_infeasible()
            if t is None:
                return
            t_cur = self.titer()
            if abs(t_cur - t) <= 1e-3 * t:
                return
            b = self.broth
            rho = b.chemicals.Water.rho('l', b.T, b.P)
            W = self.W_over_feed * F
            W_new = W + (1. / t - 1. / t_cur) * b.imass['Ethanol'] * rho
            if W_new < 0:
                if W <= 1e-6:
                    raise InfeasibleSpecification(
                        f'titer {t:.1f} g/L unattainable (max {t_cur:.1f} g/L without '
                        'saccharification water)')
                W_new = 0.
            self.W_over_feed = W_new / F
        if abs(self.titer() - t) > 5e-3 * t:
            raise RuntimeError(f'titer not converged: {self.titer():.2f} vs {t:.2f} g/L')

    def reset(self):
        bl = self.baseline
        self.W_over_feed = self.module_baseline['W_sacc'] / self.feedstock.F_mass
        self.infeasible_message = None
        self.set_try(bl['y'], bl['t'], bl['r'])
        self.recover()
        self.simulate()

    def measure(self):
        S = self._sugars()
        b = self.broth
        R303, R302 = self.R303, self.R302
        D = R303.design_results
        produced_cells = sum(s.imass['Z_mobilis'] for s in (R303.outs[1],)) - R303.ins[0].imass['Z_mobilis']
        return dict(
            F_sugar_fed=sum(S.values()),
            V_broth=b.F_vol,
            F_water_broth=b.imass['Water'],
            V_fermenters=D['Reactor volume'] * D['Number of reactors'],
            tau_fermentation=R303.tau_cofermentation,
            x_F_glucose_fed=S['Glucose'], x_F_xylose_fed=S['Xylose'],
            x_phi_max_conversion=getattr(self, 'phi', nan),
            x_glucose_to_ethanol_cofermentation=R303.cofermentation.X[0],
            x_xylose_to_ethanol_cofermentation=R303.cofermentation.X[4],
            x_saccharification_water=self.M301.ins[2].imass['Water'],
            x_tau_saccharification=R303.tau_saccharification,
            x_N_fermenters=D['Number of reactors'],
            x_cell_mass=produced_cells,
            x_residual_sugar=b.imass['Glucose'] + b.imass['Xylose'],
            x_CO2_vent=R303.outs[0].imass['CO2'] + R302.outs[0].imass['CO2'],
            x_F_ethanol_vent=R303.outs[0].imass['Ethanol'],
            x_lignin_to_boiler=self.br.BT.ins[0].F_mass,
        )

    def parameters(self):
        from biosteam.evaluation.evaluation_tools import bounded_triang, triang
        model, params = _default_parameters(self.system, self.feedstock,
                                            exclude_streams=(self.ethanol,))
        R303 = self.R303
        rxn = R303.saccharification[2]

        # As in biorefineries/cornstover/model.py (the ethanol conversion and
        # fermentation time parameters there are TRY and set by the harness)
        @model.parameter(element=R303, kind='coupled', distribution=bounded_triang(rxn.X, addition=0.04),
                         baseline=rxn.X, name='Saccharification conversion', units='glucan to glucose')
        def set_saccharification_conversion(X):
            rxn.X = X

        @model.parameter(element=R303, kind='isolated', distribution=triang(R303.tau_saccharification),
                         baseline=R303.tau_saccharification, name='Saccharification time', units='h')
        def set_saccharification_time(tau):
            R303.tau_saccharification = tau

        # Boiler natural gas price (not in the module model; BioSTEAM default +-10 %)
        BT = self.br.BT

        @model.parameter(element=BT, kind='cost', distribution=triang(BT.natural_gas_price),
                         baseline=BT.natural_gas_price, name='Natural gas price', units='USD/kg')
        def set_natural_gas_price(price):
            BT.natural_gas_price = price
        params = params + list(model.parameters[-3:])
        canonical = {
            ('cornstover', 'Price'): 'p_feedstock_price',
            ('cornstover', 'Flow rate'): 'p_feedstock_capacity',
            ('Electricity', 'Price'): 'p_electricity_price',
            ('TEA', 'Income tax'): 'p_income_tax',
            ('sulfuric_acid', 'Price'): 'p_acid_price',
            ('BT', 'Natural gas price'): 'p_natural_gas_price',
            ('R303', 'Saccharification conversion'): 'p_saccharification_yield',
            ('R303', 'Saccharification time'): 'bp_saccharification_time',
            ('natural_gas', 'Price'): 'bp_natural_gas_heating_agent_price',
        }
        return _wrap(params, canonical)

    @staticmethod
    def _assign(u):
        ID = u.ID
        if ID in ('M301', 'H301', 'R301'):
            return 'conversion'
        if ID in ('T701', 'T702', 'T703', 'P701', 'P702', 'FWT', 'CSL_storage', 'DAP_storage',
                  'Ammonia_storage', 'H2SO4_storage'):
            return 'facilities'
        if ID in ('BT',):
            return 'chp'
        if ID in ('WWTC',) or ID.startswith('WWT'):
            return 'wastewater'
        if ID in ('CWP', 'CT', 'PWC', 'ADP', 'CIP', 'CIP_package', 'FT', 'blowdown_mixer'):
            return 'facilities'
        if ID == 'M701':  # ethanol-denaturant mixer
            return 'separation'
        if ID in ('slurry_mixer', 'gas_mixer'):  # WWT sludge/biogas to boiler
            return 'chp'
        if ID == 'treated_water_mixer':
            return 'wastewater'
        if ID == 'process_water_mixer':
            return 'facilities'
        n = _numeric_area(ID)
        return {1: 'feedstock', 2: 'conversion', 3: 'fermentation', 4: 'separation',
                5: 'wastewater', 6: 'wastewater', 7: 'facilities', 8: 'facilities'}.get(n)

    def extra(self):
        br = self.br
        return dict(x_natural_gas_BT=br.BT.ins[3].F_mass)


class _Sugarcane(_Impl):
    name = 'sugarcane'
    #: Maximum fermentable sugar in the fermentor feed [g glucose-eq./L]
    #: (very-high-gravity limit; module baseline ~300 g/L). The module's own
    #: limit (first-effect evaporation V <= 0.5) corresponds to >1000 g/L.
    S_feed_max = 350.

    def load(self):
        import biosteam as bst
        import biorefineries.cane as cane
        cane.YRCP2025()
        br = cane.Biorefinery('S1')
        self.br = br
        self.system = sys = br.system
        self.tea = br.tea
        s, u = br.flowsheet.stream, br.flowsheet.unit
        self.ethanol = s.advanced_ethanol
        self.feedstock = s.sugarcane
        R = [i for i in sys.units if isinstance(i, bst.Fermentation)]
        assert len(R) == 1
        self.R302 = R = R[0]
        self.E301 = u.E301
        self.M301 = u.M301
        self.BT = u.BT401
        chems = R.chemicals
        self.theoretical_yield = 2 * chems.Ethanol.MW / chems.Glucose.MW
        self.sucrose_glc_eq = 2 * chems.Glucose.MW / chems.Sucrose.MW
        self.fermentor_ins = tuple(R.ins)
        self.fermentor_outs = tuple(R.outs)
        self.broth = R.outs[1]
        self.simulate_system()
        self.snapshot_recycles()
        self.module_baseline = self._snapshot()
        self.baseline = dict(y=self.module_baseline['y'], t=self.module_baseline['t'],
                             r=R.productivity)
        self.inflection_yield = self.theoretical_yield  # sequential cell growth on remainder
        # maximum sugar concentration [g/L broth] at the feed limit S_feed_max
        mb = self.module_baseline
        self.S_max = self.S_feed_max * (mb['S'] / mb['V_broth']) / mb['S_feed']
        # Volume additivity: 1/S_feed = y/t - y*v_EtOH - (1 - X)*k + c0, with
        # v_EtOH the ethanol specific volume [L/g], X = y/theoretical, and
        # k, c0 calibrated at the baseline and at X = 0.1, t = 12 g/L.
        chems = self.R302.chemicals
        self.v_ethanol = 1. / chems.Ethanol.rho('l', 305.15, 101325.)  # m3/kg = L/g
        self.k_volume = 0.
        self.c_volume = self._volume_lhs(mb['y'], mb['t'], mb['S_feed'])
        self._needs_reset = False
        try:
            y1, t1 = 0.1 * self.theoretical_yield, 12.
            self.set_try(y1, t1, self.baseline['r'])
            self.simulate()
            lhs1 = self._volume_lhs(y1, t1, self._feed_sugar_concentration())
            Xb = mb['y'] / self.theoretical_yield
            self.k_volume = (self.c_volume - lhs1) / (Xb - 0.1)
            self.c_volume += (1 - Xb) * self.k_volume
        except Exception:
            pass
        finally:
            self.reset()
        self.unit_groups, self.unassigned = _make_unit_groups(self.tea, self.system, self._assign)

    def _sugar_fed(self):
        return sum(i.imass['Glucose'] + self.sucrose_glc_eq * i.imass['Sucrose'] for i in self.R302.ins)

    def _feed_sugar_concentration(self):
        """Fermentable sugar [g glucose-eq./L] in the combined fermentor feed."""
        V = sum(i.F_vol for i in self.R302.ins if not i.isempty())
        return self._sugar_fed() / V

    def _snapshot(self):
        S = self._sugar_fed()
        E = self.ethanol_produced()
        mpsp, price, purity = self.mpsp()
        return dict(S=S, y=E / S, t=self.titer(), MPSP=mpsp, price=price, purity=purity,
                    MESP_gal=price * KG_PER_GAL_ETHANOL, titer_module=self.R302.titer,
                    S_feed=self._feed_sugar_concentration(), V_broth=self.broth.F_vol)

    def _volume_lhs(self, y, t, S_feed):
        return 1. / S_feed - y / t + y * self.v_ethanol

    def max_titer(self, y):
        X = y / self.theoretical_yield
        den = 1. / self.S_feed_max + y * self.v_ethanol + (1 - X) * self.k_volume - self.c_volume
        return y / den if den > 0 else float('inf')

    def set_try(self, y, t, r):
        if y > self.theoretical_yield * (1 + 1e-9):
            raise InfeasibleSpecification('yield above theoretical')
        if t > 1.05 * self.max_titer(y):
            raise InfeasibleSpecification(f'titer {t:.1f} g/L above max attainable ~{self.max_titer(y):.1f} g/L')
        R = self.R302
        R.efficiency = min(y / self.theoretical_yield, 1.)
        R.titer = t
        R.productivity = r
        self.t_target = t

    def simulate(self):
        if self._needs_reset:
            self.recover()
            self._needs_reset = False
        R = self.R302
        t = self.t_target
        try:
            # The module's titer specification is solved before the yeast
            # recycle converges (error up to ~3 %); correct its target.
            for it in range(6):
                self.simulate_system()
                t_cur = self.titer()
                if abs(t_cur - t) <= 2e-3 * t:
                    break
                R.titer *= t / t_cur
        except InfeasibleSpecification:
            self._needs_reset = True
            raise
        except Exception:
            self._needs_reset = True
            raise
        if abs(self.titer() - t) > 5e-3 * t:
            self._needs_reset = True
            raise RuntimeError(f'titer not converged: {self.titer():.2f} vs {t:.2f} g/L')
        S_feed = self._feed_sugar_concentration()
        if S_feed > self.S_feed_max * 1.005:
            raise InfeasibleSpecification(f'titer {t:.1f} g/L requires {S_feed:.0f} g/L sugar in '
                                          f'fermentor feed > {self.S_feed_max:.0f} g/L')

    def reset(self):
        bl = self.baseline
        self.set_try(bl['y'], bl['t'], bl['r'])
        self.recover()
        self._needs_reset = False
        self.simulate()

    def measure(self):
        R = self.R302
        b = self.broth
        D = R.design_results
        E = self.E301
        return dict(
            F_sugar_fed=self._sugar_fed(),
            V_broth=b.F_vol,
            F_water_broth=b.imass['Water'],
            V_fermenters=D['Reactor volume'] * D['Number of reactors'],
            tau_fermentation=R.tau,
            x_ethanol_conversion=R.efficiency,
            x_evaporator_V=E.V,
            x_evaporator_water_removed=E.outs[1].F_mass,
            x_dilution_water=self.M301.ins[-1].F_mass,
            x_N_fermenters=D['Number of reactors'],
            x_cell_mass=b.imass['Yeast'] - sum(i.imass['Yeast'] for i in R.ins),
            x_residual_sugar=b.imass['Glucose'] + b.imass['Sucrose'],
            x_CO2_vent=R.outs[0].imass['CO2'],
            x_F_ethanol_vent=R.outs[0].imass['Ethanol'],
            x_titer_module_def=R.get_titer(),
            x_feed_sugar_concentration=self._feed_sugar_concentration(),
        )

    def parameters(self):
        model = self.br.model
        want = {
            ('-', 'Cane operating days'): 'p_operating_days',
            ('feedstock', 'Available land'): 'bp_available_land',
            ('crude oil', 'Price'): 'bp_crude_oil_price',
            ('feedstock', 'Price'): 'p_feedstock_price',
            ('Natural gas', 'Price'): 'bp_natural_gas_price_residual',
            ('Electricity', 'Price'): 'bp_electricity_price_residual',
            ('-', 'IRR'): 'p_IRR',
            ('-', 'Income tax'): 'p_income_tax',
        }
        out = []
        for p in model.parameters:  # keep model order (crude oil before residuals)
            key = (p.element_name, p.name)
            if key not in want:
                continue
            name = want[key]
            if name in ('p_IRR', 'p_income_tax'):
                out.append(_pct_parameter(p, name))
            else:
                prm = Parameter.from_biosteam(p, name=name)
                if name == 'p_feedstock_price':
                    prm.description += ' (baseline cane price before transportation-cost adjustment)'
                elif name.endswith('residual'):
                    prm.description += ' (residual added to crude-oil-correlated price)'
                out.append(prm)
        assert len(out) == len(want), [o.name for o in out]
        return out

    @staticmethod
    def _assign(u):
        ID = u.ID
        if ID in ('E301', 'P301', 'H301', 'M301'):
            return 'conversion'
        if ID in ('R302', 'T301', 'C301', 'S301'):
            return 'fermentation'
        if ID in ('T302', 'P305', 'T303', 'P306', 'T305'):
            return 'facilities'  # product/denaturant storage
        n = _numeric_area(ID)
        return {1: 'feedstock', 2: 'conversion', 3: 'separation', 4: 'chp',
                5: 'facilities', 6: 'facilities', 7: 'wastewater'}.get(n)

    def extra(self):
        return dict(x_natural_gas_BT=self.BT.ins[3].F_mass)


class _Corn(_Impl):
    name = 'corn'
    allow_empty_recycles = False
    max_slurry_solids = 0.35  # wt, very-high-gravity limit

    def load(self):
        from biorefineries.corn import Biorefinery
        from biorefineries.corn.tea import default_tea_parameters
        br = Biorefinery()
        self.br = br
        self.system = br.system
        self.tea = tea = br.TEA
        self.module_IRR_breakeven = tea.IRR
        tea.IRR = default_tea_parameters['IRR']  # module's design IRR (0.15)
        self.ethanol = br.ethanol
        self.feedstock = br.corn
        self.V405, self.V307 = br.V405, br.V307
        self.rpw = br.recycled_process_water
        chems = self.V405.chemicals
        self.theoretical_yield = 2 * chems.Ethanol.MW / chems.Glucose.MW
        self.fermentor_ins = tuple(self.V405.ins)
        self.fermentor_outs = tuple(self.V405.outs)
        self.broth = self.V405.outs[1]
        self.simulate_system()
        self.snapshot_recycles()
        self.module_baseline = mb = self._snapshot()
        self.baseline = dict(y=mb['y'], t=mb['t'], r=mb['t'] / self.V405.tau)
        self.inflection_yield = self.theoretical_yield  # growth only on remaining glucose
        # install titer specification on the slurry tank
        from thermosteam.network import ProcessSpecification
        specs = self.V307._specifications
        assert len(specs) == 1 and specs[0].f.__name__ == 'correct_recycle_dilution_water'
        # Knob: total water entering the slurry tank per dry corn; the process
        # water is its complement to the other water inputs (backwater recycle,
        # corn moisture, enzymes), so the backwater feedback does not distort
        # the Newton steps.
        self.W_over_dry = self.W_over_dry_baseline = self.V307.outs[0].imass['Water'] / self._dry_corn()
        specs[0] = ProcessSpecification(self._set_recycled_water, (), ())
        self.t_target = None
        # maximum sugar concentration at the solids limit
        out = self.V307.outs[0]
        solids = out.F_mass - out.imass['Water']
        W_min_total = solids * (1 - self.max_slurry_solids) / self.max_slurry_solids
        dW = out.imass['Water'] - W_min_total
        self.S_max = mb['S'] / (mb['V_broth'] - dW / 1000.)
        self.unit_groups, self.unassigned = _make_unit_groups(
            self.tea, self.system, self._assign, _lang_factor_group_class(self.tea))

    def _dry_corn(self):
        f = self.feedstock
        return f.F_mass - f.imass['Water']

    def _other_slurry_water(self):
        return sum(i.imass['Water'] for i in self.V307.ins if i is not self.rpw)

    def _set_recycled_water(self):
        self.rpw.imass['Water'] = max(self.W_over_dry * self._dry_corn() - self._other_slurry_water(), 0.)

    def _sugar_fed(self):
        return sum(i.imass['Glucose'] for i in self.V405.ins)

    def _slurry_solids(self):
        out = self.V307.outs[0]
        return 1. - out.imass['Water'] / out.F_mass

    def _snapshot(self):
        S = self._sugar_fed()
        mpsp, price, purity = self.mpsp()
        return dict(S=S, y=self.ethanol_produced() / S, t=self.titer(), V_broth=self.broth.F_vol,
                    MPSP=mpsp, price=price, purity=purity, MESP_gal=price * KG_PER_GAL_ETHANOL,
                    slurry_solids=self._slurry_solids())

    def set_try(self, y, t, r):
        if y > self.theoretical_yield * (1 + 1e-9):
            raise InfeasibleSpecification('yield above theoretical')
        if t > 1.1 * y * self.S_max:
            raise InfeasibleSpecification(f'titer {t:.1f} g/L above max attainable ~{y*self.S_max:.1f} g/L')
        self.V405.reaction.X = min(y / self.theoretical_yield, 1.)
        self.V405.tau = t / r
        self.t_target = t
        self.W_over_dry = _predict_water(self, y, t, self.W_over_dry * self._dry_corn()) / self._dry_corn()

    def simulate(self):
        t = self.t_target
        n_no_water = 0
        for it in range(10):
            self.simulate_system()
            if t is None:
                return
            t_cur = self.titer()
            if abs(t_cur - t) <= 1e-3 * t:
                break
            if self.rpw.imass['Water'] <= 1e-6 and t_cur < t:
                n_no_water += 1
                if n_no_water >= 2:
                    raise InfeasibleSpecification(f'titer {t:.1f} g/L unattainable without process water')
            b = self.broth
            rho = b.chemicals.Water.rho('l', b.T, b.P)
            dry = self._dry_corn()
            W = self.W_over_dry * dry
            W_new = W + (1. / t - 1. / t_cur) * b.imass['Ethanol'] * rho
            self.W_over_dry = max(W_new, 0.) / dry
        else:
            if abs(self.titer() - t) > 5e-3 * t:
                raise RuntimeError(f'titer not converged: {self.titer():.2f} vs {t:.2f} g/L')
        if self._slurry_solids() > self.max_slurry_solids + 1e-4:
            raise InfeasibleSpecification(
                f'titer {t:.1f} g/L requires slurry solids {self._slurry_solids():.3f} > '
                f'{self.max_slurry_solids}')

    def reset(self):
        bl = self.baseline
        self.set_try(bl['y'], bl['t'], bl['r'])
        self.W_over_dry = self.W_over_dry_baseline
        self.recover()
        self.simulate()

    def measure(self):
        V = self.V405
        b = self.broth
        D = V.design_results
        br = self.br
        return dict(
            F_sugar_fed=self._sugar_fed(),
            V_broth=b.F_vol,
            F_water_broth=b.imass['Water'],
            V_fermenters=D['Reactor volume'] * D['Number of reactors'],
            tau_fermentation=V.tau,
            x_ethanol_conversion=V.reaction.X,
            x_recycled_process_water=self.rpw.F_mass,
            x_slurry_solids=self._slurry_solids(),
            x_N_fermenters=D['Number of reactors'],
            x_cell_mass=b.imass['Yeast'] - sum(i.imass['Yeast'] for i in V.ins),
            x_residual_sugar=b.imass['Glucose'],
            x_CO2_vent=V.outs[0].imass['CO2'],
            x_F_ethanol_vent=V.outs[0].imass['Ethanol'],
            x_DDGS=br.DDGS.F_mass,
            x_crude_oil=br.crude_oil.F_mass,
            x_module_IRR_breakeven=self.module_IRR_breakeven,
        )

    def parameters(self):
        model, params = _default_parameters(self.system, self.feedstock,
                                            exclude_streams=(self.ethanol,))
        canonical = {
            ('corn', 'Price'): 'p_feedstock_price',
            ('corn', 'Flow rate'): 'p_feedstock_capacity',
            ('Electricity', 'Price'): 'p_electricity_price',
            ('TEA', 'Income tax'): 'p_income_tax',
            ('sulfuric_acid', 'Price'): 'p_acid_price',
        }
        return _wrap(params, canonical)

    _feedstock_units = ('MH101', 'V102', 'MH103', 'M104', 'V105', 'W106', 'V107')
    _fermentation_units = ('E401', 'E402', 'V403', 'P404', 'V405', 'P406', 'MX', 'P407', 'V412', 'PX')
    _storage_units = ('V511', 'P512', 'V509', 'P510', 'V513')

    @classmethod
    def _assign(cls, u):
        ID = u.ID
        if ID in cls._feedstock_units:
            return 'feedstock'
        if ID in cls._fermentation_units:
            return 'fermentation'
        if ID in cls._storage_units:
            return 'facilities'
        if ID in ('V409', 'P410', 'P301', 'E413'):
            return 'separation'  # vent scrubber, beer pump, beer column heat exchange
        if ID == 'HX101':
            return 'conversion'  # mash heater before liquefaction
        if ID == 'MX4':
            return 'separation'  # ethanol-denaturant mixer
        if ID == 'S1':
            return 'wastewater'  # thin stillage splitter (backwater)
        if ID in ('other_facilities', 'PWC', 'CT', 'CWP'):
            return 'facilities'
        n = _numeric_area(ID)
        return {3: 'conversion', 5: 'separation', 6: 'wastewater'}.get(n)

    def extra(self):
        return {}


def _wrap(params, canonical):
    out = []
    used = set()
    from ..harness import _slug
    for p in params:
        key = _param_key(p)
        name = canonical.get(key)
        if name is None:
            name = 'bp_' + _slug(f'{key[0]} {key[1]}')
        base, k = name, 2
        while name in used:
            name = f'{base}_{k}'
            k += 1
        used.add(name)
        out.append(Parameter.from_biosteam(p, name=name))
    return out


_IMPLS = {'ethanol_cornstover': _Cornstover, 'ethanol_sugarcane': _Sugarcane,
          'ethanol_corn': _Corn}


# %% Adapter

class Adapter(BiorefineryAdapter):
    name = 'ethanol'
    area_map = {a: a for a in CANONICAL_AREAS}

    def configurations(self):
        return [dict(c) for c in CONFIGURATIONS]

    def load(self, config):
        warnings.filterwarnings('ignore')
        impl = _IMPLS[config['config_id']](self)
        impl.load()
        self.impl = impl
        self.config = config
        self.system = impl.system
        self.tea = impl.tea
        self.unit_groups = impl.unit_groups
        self.feedstock_stream = impl.feedstock
        self.product_stream = impl.ethanol

    def domain(self):
        impl = self.impl
        th = impl.theoretical_yield
        infl = impl.inflection_yield
        # cornstover: 0.5 % headroom below the inflection (it depends slightly on
        # the parameter-dependent glucose/xylose ratio); sugarcane/corn: the
        # inflection is the theoretical yield (zero cell growth), a degenerate
        # state where the WWT sludge solver fails; use 99 %.
        y_max = 0.995 * infl if isinstance(impl, _Cornstover) else 0.99 * infl
        # lowest yield at which titers >= 1.15 t_min are attainable
        lo, hi = 0.05 * th, y_max
        if 0.98 * impl.max_titer(lo) < 1.15 * impl.t_min:
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if 0.98 * impl.max_titer(mid) < 1.15 * impl.t_min:
                    lo = mid
                else:
                    hi = mid
            lo = hi
        y_min = lo
        t_max = min(impl.t_max_cap, math.floor(0.98 * impl.max_titer(y_max)))
        bl = impl.baseline
        return dict(theoretical_yield=th, inflection_yield=infl,
                    y_min=y_min, y_max=y_max, t_min=impl.t_min, t_max=float(t_max),
                    baseline_productivity=bl['r'], baseline_yield=bl['y'],
                    baseline_titer=bl['t'],
                    x_max_sugar_concentration=impl.S_max,
                    x_module_baseline_MPSP=impl.module_baseline['MPSP'],
                    x_module_baseline_MESP_USD_per_gal_denatured=impl.module_baseline['MESP_gal'])

    def max_titer(self, y):
        return 0.98 * self.impl.max_titer(y)

    def parameters(self):
        return self.impl.parameters()

    def set_try(self, y, t, r):
        self.impl.set_try(y, t, r)

    def simulate(self):
        self.impl.simulate()

    def reset(self):
        self.impl.reset()

    def collect(self):
        import biosteam as bst
        impl = self.impl
        tea = self.tea
        eth = impl.ethanol
        MPSP, price, purity = impl.mpsp()
        m = impl.measure()
        produced = impl.ethanol_produced()
        F_final = eth.imass['Ethanol']
        units = tea.units
        heating, cooling = _heat_duties(self.system.units)
        group = bst.UnitGroup('all', units)
        ng_streams = [s for s in self.system.feeds if 'natural_gas' in s.ID]
        BT = getattr(impl, 'BT', None)
        if BT is None:
            BT = getattr(impl.br, 'BT', None)
        if BT is not None and hasattr(BT, 'natural_gas_price'):
            ng_price = BT.natural_gas_price
        else:
            ng_price = ng_streams[0].price if ng_streams else nan  # corn: dryer NG unpriced (0)
        b = impl.broth
        out = dict(
            MPSP=MPSP,
            recovery=F_final / produced,
            F_final_product=F_final,
            F_ferm_product=produced,
            F_feedstock=impl.feedstock.F_mass,
            feedstock_price=impl.feedstock.price,
            base_cost=0., acid_cost=0.,
            heating_duty=heating, cooling_duty=cooling,
            electricity_consumption=group.get_electricity_consumption(),
            electricity_production=group.get_electricity_production(),
            natural_gas_flow=sum(s.F_mass for s in ng_streams),
            electricity_price=bst.PowerUtility.price,
            natural_gas_price=ng_price,
            GWP=nan, FEC=nan,
            x_yield_measured=produced / m['F_sugar_fed'],
            x_titer_broth=b.imass['Ethanol'] / b.F_vol,
            x_F_ethanol_broth=b.imass['Ethanol'],
            x_product_purity=purity,
            x_F_product_stream=eth.F_mass,
            x_MPSP_denatured=price,
            x_MESP_USD_per_gal_denatured=price * KG_PER_GAL_ETHANOL,
        )
        out.update(m)
        out.update(impl.extra())
        return out
