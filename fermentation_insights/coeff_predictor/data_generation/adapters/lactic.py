# -*- coding: utf-8 -*-
"""
Adapter for the lactic acid biorefinery (``biorefineries.lactic``, Li et al.,
ACS Sustainable Chem. Eng. 2021, 9, 1341-1351): corn stover, dilute-acid
pretreatment, enzymatic hydrolysis, co-fermentation of glucose and xylose to
lactic acid, purification by esterification with ethanol / distillation /
hydrolysis to 88 wt% lactic acid.

Configurations (all corn stover, batch fermentation, 50 C, anaerobic):

=============================  =====  ===============================================
config_id                      conv.  fermentation regime
=============================  =====  ===============================================
lactic_neutral_SSCF_cornstover SSCF   regular strain, lime neutralization (Ca-lactate),
                                      H2SO4 acidulation + gypsum filtration (default)
lactic_neutral_SHF_cornstover  SHF    as above; saccharified sugars can be concentrated
                                      by a multi-effect evaporator (E301)
lactic_lowpH_SSCF_cornstover   SSCF   acid-resistant strain, no neutralization;
                                      acidulation (R401) and gypsum filter (S402) bypassed,
                                      316 SS fermenters
lactic_lowpH_SHF_cornstover    SHF    acid-resistant strain with sugar concentration
=============================  =====  ===============================================

The regimes and TRY handling follow ``analyses/2_fermentation_TRY.py`` (and the
paper-era scripts ``2-1_ferm_diluted.py``/``2-2_ferm_concentrated.py`` at
commit f56692d, which additionally set ``R301.neutralization = False`` for the
acid-resistant strain): yield is set with ``utils.set_yield`` (lactic acid
conversion of glucose and xylose, = g/g since both reactions are mass
conserving), titer by ``R301.target_titer`` with dilution allowed (SSCF, SHF)
and concentration allowed (SHF), productivity by ``R301.productivity`` (seed
train at ``ferm_ratio`` x). A titer above the attainable maximum (no dilution
water for SSCF, 220 g/L initial sugar in batch for SHF) is infeasible (the TRY
script flags ``effluent_titer + 2 < target``).

MPSP is ``tea.solve_price(lactic_acid)`` (iterated 3x from price 0 as in
``models.get_MPSP``) divided by the lactic acid mass fraction of the product
stream (0.88), i.e. USD/kg pure lactic acid. Recovery = lactic acid in product
/ (lactic acid + Ca-lactate as acid) leaving R301 (``models.get_recovery``).

Domain: yields 0.30 g/g (TRY scripts) to the inflection yield
1 - X_acetic - X_cell = 0.91 g/g at the baseline byproduct conversions (above it
``set_yield`` lowers the acetic acid conversion); titers 20-200 g/L (TRY script
max 200 g/L; min lowered from 40-50 g/L so that SSCF still has a usable titer
span at low yields, where the undiluted titer is only ~0.3 x 125 = 37 g/L);
productivity baseline 0.89 g/L/h. ``max_titer(y)`` (98% of the attainable
titer) is calibrated per parameter sample by simulating only the upstream
(pretreatment + conversion) section: SSCF at each yield without dilution water;
SHF at y_min and y_max by bisection up to the 220 g/L initial-sugar limit (the
module's bracketed concentration solver must not get unattainable targets),
interpolating titer/yield in between. The sampled acetic acid yield
(``bp_acetic_acid_yield``) moves the sample's own inflection yield
(``x_inflection_yield``; points above it have ``x_above_inflection`` = True).
Measured yield (fresh glucose + xylose basis) exceeds y by up to ~4% at low y
because residual sugars of the seed train are recycled with the inoculum.

Robustness/consistency measures (see also Bioindustrial-Park branch fi-lactic,
which ports lactic to BioSTEAM 2.54 and fixes the titer specifications): the
separation recycle loop uses Aitken (``LACTIC_SEP_METHOD``), every point is
re-simulated until key outputs change < 0.2% (``x_n_simulations``,
``x_stability_rel_change``), and ethyl acetate is purged from the R402-D402
ethanol recycle (no outlet in the module; ~1-2 kg/h formed).

Parameters are the module's own (``models.create_model``) except yield and
productivity, plus feedstock carbohydrate content and feedstock price with the
ranges of ``analyses/3_feedstock_carb-price.py`` (25-70 dw% carbohydrates,
0-200 $/dry-ton).
"""
from __future__ import annotations

import os
import time

import numpy as np

from ..harness import (BiorefineryAdapter, InfeasibleSpecification, Parameter,
                       describe_feedstock)

THEORETICAL_YIELD = 1.0         # g lactic acid / g glucose (or xylose)
Y_MIN = 0.30                    # g/g, TRY scripts (np.arange(0.3, 1.01, ...))
T_MIN, T_MAX = 40., 200.        # g/L; module design range starts at ~40-50 g/L.
# Below ~40 g/L the evaporator losses make recovery collapse (0.23 at 20 g/L),
# MPSP leaves the unifying-equation form (fit MAPE 16% vs 1.5% for t >= 40).
BASELINE_YIELD = 0.76           # g/g (R301.target_yield default)
BASELINE_PRODUCTIVITY = 0.89    # g/L/h
TITER_TOL_ABS = 1.0             # g/L (TRY script: effluent_titer + 2 < target)
TITER_TOL_REL = 0.01
MAX_TITER_MARGIN = 0.98         # top titer of the design = 98% of attainable
SHF_MAX_SUGAR = 220.            # g/L, initial sugar limit in batch (adjust_ferm_loop)
CARBS = ('Glucan', 'Xylan', 'Arabinan', 'Galactan', 'Mannan')
FERMENTABLE = ('Glucose', 'Xylose')  # only these are fermented in R301/R302
CARB_RANGE = (0.25, 0.70)       # dw fraction, 3_feedstock_carb-price.py
FEED_PRICE_RANGE = (0., 200.)   # $/dry-ton, 3_feedstock_carb-price.py
ETAC_MAX_FRACTION = 0.05        # EtAc/EtOH in the D402 ethanol recycle (sanity)
# Separation recycle loop (S403/D402/D404/F402 recycles) solver. The module sets
# fixed-point globally ('aitken isn't stable' in old BioSTEAM); in BioSTEAM 2.54
# fixed-point stops ~2% short of the slowly converging 97% S403 recycle after
# large TRY jumps (history-dependent MPSP up to +20%) and is ~5x slower, whereas
# Aitken converges it tightly (verified against fixed-point restarts).
SEP_LOOP_METHOD = os.environ.get('LACTIC_SEP_METHOD', 'aitken')
MAX_SIMULATIONS = 4             # per TRY point (re-simulate until stable)
STABILITY_RTOL = 2e-3           # max rel. change of key outputs between simulations

_UTILITY_PRODUCERS = ('BoilerTurbogenerator', 'CoolingTower', 'ChilledWaterPackage')

#: module parameter name (as formatted by biosteam) -> canonical name
_CANONICAL = {
    'Sulfuric acid unit price': 'p_acid_price',
    'Lime unit price': 'p_base_price',
    'Natural gas unit price': 'p_natural_gas_price',
    'Electricity unit price': 'p_electricity_price',
    'Enzyme loading': 'p_enzyme_loading',
    'Enzymatic hydrolysis glucan-to-glucose': 'p_saccharification_yield',
    'Boiler efficiency': 'p_boiler_efficiency',
}
#: excluded: TRY (set by the harness) and the "blank" sensitivity placeholder
_EXCLUDED = ('Blank parameter', 'Lactic acid yield', 'Productivity')
#: parameters without effect when acidulation/gypsum filtration are bypassed
_NEUTRAL_ONLY = ('Gypsum unit price', 'Gypsum split', 'Acidulation time')


def _norm(name):
    return ' '.join(str(name).lower().split())


class Adapter(BiorefineryAdapter):
    name = 'lactic'
    area_map = {a: a for a in ('feedstock', 'conversion', 'fermentation', 'separation',
                               'upgrading', 'wastewater', 'chp', 'facilities')}

    def __init__(self):
        self.config = None
        self._pvals = {}
        self._kcache = {}
        self._tci_ratio = 1.
        self._target = (BASELINE_YIELD, None, BASELINE_PRODUCTIVITY)

    # %% Configurations
    def configurations(self):
        configs = []
        for neutral in (True, False):
            for kind in ('SSCF', 'SHF'):
                regime = 'neutral' if neutral else 'lowpH'
                configs.append(dict(
                    config_id=f'lactic_{regime}_{kind}_cornstover',
                    feedstock='cornstover',
                    feedstock_type='lignocellulose',
                    fermentation_product='LacticAcid',
                    final_product='LacticAcid',  # 88 wt% grade; MPSP per kg pure
                    ferm_product_chemical='LacticAcid',
                    final_product_chemical='LacticAcid',
                    neutralization=neutral,
                    neutralizing_agent='Ca(OH)2' if neutral else '',
                    n_acid_groups=1,
                    pKa1=3.86,
                    final_product_stoich_mass_ratio=1.0,
                    upgrading=False,
                    separation_variant='' if neutral else 'no_acidulation',
                    separation_principle=(
                        ('acidulation(H2SO4)+gypsum_filtration+' if neutral else '')
                        + 'flash+esterification(EtOH)+distillation+hydrolysis'),
                    aerobic=False,
                    fermentation_mode='batch',
                    T_fermentation=323.15,
                    product_water_solubility=1000.,
                    # lactic-specific design descriptors
                    saccharification=kind,  # SSCF: simultaneous; SHF: separate
                    sugar_concentration=(kind == 'SHF'),
                    product_purity_spec=0.88,
                ))
        return configs

    # %% Loading
    def load(self, config):
        import biosteam as bst
        from biorefineries import lactic as la
        self.config = config
        kind = config['saccharification']
        if not la._chemicals_loaded:
            la._load_chemicals()
        flowsheet = bst.Flowsheet(f"lactic_{config['config_id']}")
        bst.main_flowsheet.set_flowsheet(flowsheet)
        bst.settings.set_thermo(la.chemicals)
        la.load_process_settings()
        self._restore_tci_ratio()
        system = la.create_system(kind=kind, flowsheet=flowsheet)
        tea = la.create_tea(flowsheet=flowsheet)
        self.flowsheet, self.system, self.tea = flowsheet, system, tea
        s, u = flowsheet.stream, flowsheet.unit
        self.s, self.u = s, u
        self.R301, self.R302 = u.R301, u.R302
        self.feedstock_stream = s.feedstock
        self.product_stream = s.lactic_acid
        self.base_stream = s.lime              # -> T605 -> lime_R301
        self.acid_stream = s.sulfuric_acid     # -> T602 -> sulfuric_acid_R401
        self.natural_gas = s.natural_gas
        self.BT = u.BT
        self.sep_loop = self._find_sep_loop(system)
        self._add_ethyl_acetate_purge(u.R402)
        if SEP_LOOP_METHOD and self.sep_loop is not None:
            self.sep_loop.method = SEP_LOOP_METHOD
        # Fermentation regime
        neutral = config['neutralization']
        self.R301.neutralization = neutral
        u.R401.bypass = u.S402.bypass = not neutral
        if kind == 'SHF':  # CoFermentation (SSCF unit is batch only)
            self.R301.mode = 'batch'
            self.R301.feed_freq = 1
        # Module-default operation (no dilution/concentration) -> baseline
        t0 = time.time()
        self._simulate_robust()
        self.baseline_titer = float(self.R301.effluent_titer)
        base_stream_price = self._solve_stream_price()
        purity = self.product_stream.imass['LacticAcid'] / self.product_stream.F_mass
        print(f'[lactic] {config["config_id"]}: module-default baseline MPSP '
              f'{base_stream_price:.4f} $/kg product ({purity:.3f} LA) = '
              f'{base_stream_price / purity:.4f} $/kg pure LA; titer '
              f'{self.baseline_titer:.2f} g/L; yield {self.R301.lactic_yield:.3f}; '
              f'TCI {tea.TCI / 1e6:.1f} MM$ ({time.time() - t0:.0f} s)', flush=True)
        self.baseline_MPSP_stream = base_stream_price
        # TRY-analysis operation
        self.R301.allow_dilution = True
        if hasattr(self.R301, 'allow_concentration'):
            self.R301.allow_concentration = (kind == 'SHF')
        self.funcs = la.create_funcs(lactic_tea=tea, flowsheet=flowsheet)
        self.model = la.create_model(flowsheet=flowsheet, kind=kind)  # simulates
        self.unit_groups = self._make_unit_groups(bst)
        self._utility_producers = [i for i in system.units
                                   if type(i).__name__ in _UTILITY_PRODUCERS]
        self._process_units = [i for i in system.units if i not in self._utility_producers]
        self._pvals = {}
        self._kcache = {}
        self._params = None
        self._target = (BASELINE_YIELD, self.baseline_titer, BASELINE_PRODUCTIVITY)
        self._k_default = self._calibrate_max_titer()

    @staticmethod
    def _find_sep_loop(system):
        def search(sy):
            for ss in sy.subsystems:
                if ss.recycle is not None and 'D403_l_recycled' in str(ss.recycle):
                    return ss
                r = search(ss)
                if r is not None:
                    return r
        return search(system)

    def _make_unit_groups(self, bst):
        u, sy = self.u, self.flowsheet.system
        kind = self.config['saccharification']
        conv_units = list(sy.conversion_sys.units)
        upstream = {'M301', 'H301', 'R300', 'S301'}  # enzymatic hydrolysis (SHF: + solids removal)
        groups = {
            'feedstock': list(sy.preprocessing_sys.units),
            'conversion': list(sy.pretreatment_sys.units) + [i for i in conv_units if i.ID in upstream],
            # SSCF R301 also contains the 24-h saccharification tank (see x_R301_saccharification_installed_cost)
            'fermentation': [i for i in conv_units if i.ID not in upstream] + [u.T604, u.T605],
            'separation': list(sy.separation_sys.units) + [u.T602, u.T606, u.T606_P],
            'wastewater': list(sy.wastewater_sys.units) + [u.T603],
            'chp': [u.BT],
            'facilities': [u.HXN, u.CT, u.FT, u.PWC, u.ADP, u.CIP, u.T601, u.T601_P, u.M601],
        }
        # chemical storage tanks are attributed to the consuming area so that the
        # purchased chemical (lime, CSL, H2SO4, ethanol, ammonia) -- which enters
        # the storage tank first -- is counted in that area's material cost
        seen = set()
        for k, v in groups.items():
            groups[k] = [i for i in v if not (i in seen or seen.add(i))]
        missing = [i for i in self.system.units if i not in seen]
        if missing:
            print(f'[lactic] units not grouped -> facilities: {[i.ID for i in missing]}', flush=True)
            groups['facilities'] += missing
        self._kind = kind
        return [bst.UnitGroup(k, v) for k, v in groups.items()]

    # %% Domain
    def domain(self):
        inflection = 1. - float(self.R301._X[1]) - float(self.R301._X[2])
        return dict(
            theoretical_yield=THEORETICAL_YIELD,
            inflection_yield=inflection,
            y_min=Y_MIN,
            y_max=inflection,
            t_min=T_MIN, t_max=T_MAX,
            baseline_productivity=BASELINE_PRODUCTIVITY,
            baseline_yield=BASELINE_YIELD,
            baseline_titer=self.baseline_titer,
        )

    def _param_key(self):
        return tuple(sorted((k, round(float(v), 10)) for k, v in self._pvals.items()))

    def _calibrate_max_titer(self, y=0.6):
        """Maximum attainable titer [g/L] at yield y from an upstream-only
        (preprocessing + pretreatment + conversion) simulation with an
        unattainable titer target (SSCF: no dilution water; SHF: sugars
        concentrated to the 220 g/L batch limit). Returns the ratio
        titer / y [g/L per g/g] (weakly decreasing with y)."""
        from biorefineries.lactic.utils import set_yield
        R301, R302 = self.R301, self.R302
        sy = self.flowsheet.system
        old = (R301.target_yield, R301.target_titer)

        def run_upstream(target):
            R301.target_titer = target
            for sub in (sy.preprocessing_sys, sy.pretreatment_sys, sy.conversion_sys):
                sub.simulate()
            return float(R301.effluent_titer)

        try:
            R301.target_yield = y
            set_yield(y, R301, R302)
            if not self.config['sugar_concentration']:
                # unattainable target -> no dilution water -> maximum titer
                titer = run_upstream(1e4)
            else:
                # The concentration specification needs a bracketed target, so
                # start from a surely attainable titer and scale it by the
                # 220 g/L initial-sugar limit (titer ~ initial sugar), then verify
                # (approached from below: an unattainable target makes the
                # module's bracketed solver extrapolate)
                lo = titer = run_upstream(100. * y)  # attainable (dilution or mild concentration)
                hi = None                            # smallest failed target
                for i in range(10):
                    sugar = float(R301.max_sugar)
                    if sugar >= 0.995 * SHF_MAX_SUGAR or (hi is not None and hi - lo < 3e-3 * lo):
                        break
                    if hi is None:
                        est = lo * SHF_MAX_SUGAR / sugar  # overestimates (sublinear)
                        target = lo + 0.9 * (est - lo)
                    else:
                        target = 0.5 * (lo + hi)
                    new = run_upstream(target)
                    if abs(new - target) > max(0.2, 2e-3 * target) or R301.max_sugar > 1.002 * SHF_MAX_SUGAR:
                        hi = target
                        run_upstream(lo)  # restore the last attainable state
                    else:
                        lo = new
                titer = lo
            k = titer / float(R301.lactic_yield)
        finally:
            R301.target_yield, R301.target_titer = old
            set_yield(old[0], R301, R302)
        return k

    def _k(self, y):
        key = (self._param_key(), round(float(y), 6))
        k = self._kcache.get(key)
        if k is None:
            try:
                k = self._calibrate_max_titer(y)
            except Exception as e:  # never raise from max_titer (harness calls it unguarded)
                print(f'[lactic] max-titer calibration failed ({type(e).__name__}: {e}); '
                      f'using load-time value', flush=True)
                k = getattr(self, '_k_default', 125.)
            self._kcache[key] = k
        return k

    def _attainable_titer(self, y):
        y = float(y)
        if not self.config['sugar_concentration']:
            return self._k(y) * y  # cheap (~0.2 s): exact at each yield
        # SHF (each calibration ~5-10 s): calibrate at the domain ends and
        # interpolate titer/yield linearly (t_max(y) is concave -> conservative)
        y0, y1 = Y_MIN, self.domain()['y_max']
        k0, k1 = self._k(y0), self._k(y1)
        w = min(max((y - y0) / (y1 - y0), 0.), 1.)
        return (k0 + w * (k1 - k0)) * y

    def max_titer(self, y):
        return MAX_TITER_MARGIN * self._attainable_titer(y)

    # %% Parameters
    def parameters(self):
        if self._params is not None:
            return self._params
        params = []
        neutral = self.config['neutralization']
        names = {_norm(p.name): p for p in self.model.parameters}
        for p in self.model.parameters:
            n = p.name
            if any(_norm(n) == _norm(x) for x in _EXCLUDED):
                continue
            if not neutral and any(_norm(n) == _norm(x) for x in _NEUTRAL_ONLY):
                continue
            nn = _norm(n)
            if nn == 'plant uptime':
                params.append(self._wrap(Parameter(
                    name='p_operating_days', baseline=365. * float(p.baseline),
                    inverse_cdf=lambda u, d=p.distribution: 365. * float(d.inv(u)),
                    setter=self._set_operating_days,
                    description='TEA: Plant uptime [fraction] x 365 -> operating days [d/y]')))
            elif nn == 'tci ratio':
                params.append(self._wrap(Parameter(
                    name='bp_tci_ratio', baseline=float(p.baseline),
                    inverse_cdf=lambda u, d=p.distribution: float(d.inv(u)),
                    setter=self._set_tci_ratio,
                    description='TEA: TCI ratio [fraction of baseline] (scales decorated '
                                'cost items; each shared class-level cost item scaled once)')))
            elif nn == 'feedstock unit price':
                continue  # replaced below by the 3_feedstock_carb-price range
            else:
                prm = Parameter.from_biosteam(p, name=_CANONICAL.get(p.name))
                params.append(self._wrap(prm))
        from biorefineries.lactic.utils import feedstock_factor
        fs = self.feedstock_stream
        lo, hi = FEED_PRICE_RANGE
        params.append(self._wrap(Parameter(
            name='p_feedstock_price', baseline=float(fs.price),
            inverse_cdf=lambda u: (lo + (hi - lo) * float(u)) / feedstock_factor,
            setter=self._set_feedstock_price,
            description=f'TEA: Feedstock unit price [USD/kg wet] = U({lo:g}, {hi:g}) $/dry-ton '
                        f'/ {feedstock_factor:.2f} (3_feedstock_carb-price.py; models.py '
                        f'baseline 71.3 $/dry-ton)')))
        self._carb_baseline = self._get_carbs()
        lo, hi = CARB_RANGE
        params.append(self._wrap(Parameter(
            name='bp_feedstock_carbohydrate_content', baseline=self._carb_baseline,
            inverse_cdf=lambda u: lo + (hi - lo) * float(u),
            setter=self._set_carbs,
            description=f'Feedstock: carbohydrate (glucan+xylan+arabinan+galactan+mannan) '
                        f'content [dry-weight fraction] = U({lo:g}, {hi:g}), balance with '
                        f'extractives (3_feedstock_carb-price.py set_carbs)')))
        self._params = params
        for prm in params:
            self._pvals[prm.name] = prm.baseline
        return params

    def _wrap(self, prm):
        setter = prm.setter
        name = prm.name
        def set_and_record(value, setter=setter, name=name):
            setter(value)
            self._pvals[name] = float(value)
        prm.setter = set_and_record
        return prm

    def _set_operating_days(self, days):
        self.tea.operating_days = float(days)

    def _set_feedstock_price(self, price):
        self.feedstock_stream.price = float(price)

    def _set_tci_ratio(self, new_ratio):
        old = self._tci_ratio
        seen = set()
        for unit in self.system.units:
            items = getattr(unit, 'cost_items', None)
            if not items:
                continue
            for item in items.values():
                if id(item) in seen:
                    continue
                seen.add(id(item))
                item.cost = item.cost / old * new_ratio
        self._tci_ratio = float(new_ratio)

    def _restore_tci_ratio(self):
        if getattr(self, 'system', None) is not None and self._tci_ratio != 1.:
            self._set_tci_ratio(1.)
        self._tci_ratio = 1.

    def _get_carbs(self):
        fs = self.feedstock_stream
        return float(fs.imass[CARBS].sum() / (fs.F_mass - fs.imass['H2O']))

    def _set_carbs(self, content):
        # analyses/3_feedstock_carb-price.py: set_carbs
        fs = self.feedstock_stream
        dry_mass = fs.F_mass - fs.imass['H2O']
        old_total = fs.imass[CARBS].sum()
        ratio = fs.get_normalized_mass(CARBS)
        new = dry_mass * float(content) * ratio
        fs.set_flow(new, 'kg/hr', CARBS)
        fs.imass['Extractives'] -= new.sum() - old_total
        if any(fs.mass < 0):
            raise ValueError(f'Carbohydrate content of {content*100:.1f} dw% is infeasible')

    # %% Simulation
    def set_try(self, y, t, r):
        from biorefineries.lactic.utils import set_yield
        y, t, r = float(y), float(t), float(r)
        if not (0. < y <= THEORETICAL_YIELD):
            raise InfeasibleSpecification(f'yield {y:.4g} g/g outside (0, {THEORETICAL_YIELD}]')
        t_att = self._attainable_titer(y)
        # SHF: targets above the attainable titer make the module's bracketed
        # concentration solver extrapolate, so screen them strictly
        t_lim = (t_att * 1.002 if self.config['sugar_concentration']
                 else t_att * (1 + 2 * TITER_TOL_REL) + TITER_TOL_ABS)
        if t > t_lim:
            raise InfeasibleSpecification(
                f'titer {t:.4g} g/L above attainable {t_att:.4g} g/L at yield {y:.3g} '
                f'({"220 g/L sugar limit" if self.config["sugar_concentration"] else "no dilution"})')
        R301, R302 = self.R301, self.R302
        R301.target_yield = y
        R301.target_titer = t
        set_yield(y, R301, R302)
        R301.productivity = r
        R302.productivity = r * R302.ferm_ratio
        self._target = (y, t, r)

    @staticmethod
    def _add_ethyl_acetate_purge(R402):
        # Ethyl acetate (formed from ~1 kg/h acetic acid in R402) has no outlet
        # from the R402-D402 ethanol recycle when all D402 distillate is recycled
        # (ethanol short), so its steady state does not exist and its level
        # depends on the number of recycle iterations / Aitken extrapolation
        # (up to ~1e6 kg/h). Purge it from the recycled ethanol each time R402
        # runs (equivalent to a lights purge of the ~1-2 kg/h formed).
        ethanol1 = R402.ins[1]
        @R402.add_specification(run=True)
        def purge_ethyl_acetate_from_ethanol_recycle():
            if ethanol1.imass['EthylAcetate']:
                ethanol1.imass['EthylAcetate'] = 0.

    def _zero_ethyl_acetate(self):
        # Ethyl acetate (from acetic acid in R402) has no outlet from the
        # R402-D402 ethanol recycle (all D402 distillate is recycled when ethanol
        # is short), so it accumulates with every recycle iteration; reset it before
        # each simulation so results do not depend on the simulation history
        # (it stays at ~1-100 kg/h, negligible vs. ~90,000 kg/h ethanol).
        try:
            self.u.R402.ins[1].imass['EthylAcetate'] = 0.
        except Exception:
            pass

    def _sane(self):
        e1 = self.u.R402.ins[1]
        etac, etoh = e1.imass['EthylAcetate'], e1.imass['Ethanol']
        if etac > ETAC_MAX_FRACTION * max(etoh, 1.):
            raise RuntimeError(f'ethyl acetate accumulated in ethanol recycle ({etac:.0f} kg/h)')
        if not np.isfinite(self.product_stream.F_mass) or self.product_stream.imass['LacticAcid'] <= 0:
            raise RuntimeError('no lactic acid product')

    def _simulate_once(self, zero_etac=True):
        if zero_etac:
            self._zero_ethyl_acetate()
        self.system.simulate()
        self._sane()

    def _simulate_robust(self, zero_etac=True):
        """Simulate with the module's 'bugfix barrage': reset recycles/cache
        (as in the TRY script's exception handler) and switch the separation
        recycle solver."""
        try:
            return self._simulate_once(zero_etac)
        except InfeasibleSpecification:
            raise
        except Exception as e:
            first = e
        sep = self.sep_loop
        method0 = sep.method if sep is not None else None
        for method in (method0, 'wegstein', 'aitken'):
            try:
                self.system.empty_recycles()
                self.system.reset_cache()
                if sep is not None and method:
                    sep.method = method
                return self._simulate_once()
            except InfeasibleSpecification:
                raise
            except Exception:
                continue
            finally:
                if sep is not None and method0:
                    sep.method = method0
        raise first

    def _state_metric(self):
        from biosteam.process_tools import utils as ptu
        hus = [hu for i in self._process_units for hu in i.heat_utilities]
        return np.array([self.product_stream.imass['LacticAcid'], self.natural_gas.F_mass,
                         ptu.get_heating_duty(hus), self.u.S403.outs[0].F_mass,
                         float(self.R301.effluent_titer)])

    def simulate(self):
        """Simulate until the titer specification is met and key outputs are
        stable between consecutive simulations (the dilution/concentration
        specifications act on a single pass of the seed loop, and the 97%
        S403 lactic acid recycle converges slowly), at most MAX_SIMULATIONS."""
        y, t, r = self._target
        prev = None
        change = float('nan')
        for k in range(MAX_SIMULATIONS):
            self._simulate_robust(zero_etac=(k == 0))
            m = self._state_metric()
            titer_ok = t is None or abs(m[-1] - t) <= max(TITER_TOL_ABS, TITER_TOL_REL * t)
            if prev is not None:
                change = float(np.max(np.abs(m - prev) / np.maximum(np.abs(prev), 1.)))
                if change <= STABILITY_RTOL and (titer_ok or k >= 2):
                    break
            prev = m
        self._n_sims, self._stability = k + 1, change
        if t is None:
            return
        titer = float(self.R301.effluent_titer)
        tol = max(TITER_TOL_ABS, TITER_TOL_REL * t)
        if titer < t - tol:
            raise InfeasibleSpecification(
                f'titer target {t:.4g} g/L not attainable at yield {y:.3g} (got {titer:.4g})')
        if titer > t + tol:
            raise RuntimeError(f'titer {titer:.4g} g/L above target {t:.4g} (dilution failed)')
        if abs(float(self.R301.lactic_yield) - y) > 1e-4:
            raise RuntimeError(f'yield {self.R301.lactic_yield:.4g} != target {y:.4g}')

    def reset(self):
        from biorefineries.lactic.utils import set_yield
        R301, R302 = self.R301, self.R302
        try:
            self.system.empty_recycles()
            self.system.reset_cache()
        except Exception:
            pass
        R301.target_yield = BASELINE_YIELD
        R301.target_titer = min(self.baseline_titer, 0.95 * self._k_default * BASELINE_YIELD)
        set_yield(BASELINE_YIELD, R301, R302)
        R301.productivity = BASELINE_PRODUCTIVITY
        R302.productivity = BASELINE_PRODUCTIVITY * R302.ferm_ratio
        self._target = (BASELINE_YIELD, None, BASELINE_PRODUCTIVITY)
        self._simulate_robust()

    # %% Results
    def _solve_stream_price(self):
        product, tea = self.product_stream, self.tea
        product.price = 0.
        for _ in range(3):  # models.get_MPSP / TRY script solve_TEA
            product.price = tea.solve_price(product)
        return float(product.price)

    def _sugar_fed(self):
        """Fresh fermentable sugars (glucose + xylose, kg/h) fed to fermentation:
        SSCF: after in-situ saccharification in R301 minus the residual sugars
        recycled with the inoculum; SHF: saccharified (and concentrated) feeds."""
        R301 = self.R301
        if self.config['saccharification'] == 'SSCF':
            return float(R301.saccharified_stream.imass[FERMENTABLE].sum()
                         - R301.ins[1].imass[FERMENTABLE].sum())
        return float(R301.ins[0].imass[FERMENTABLE].sum() + R301.ins[5].imass[FERMENTABLE].sum())

    def collect(self):
        import biosteam as bst
        from biosteam.process_tools import utils as ptu
        u, s, R301, R302 = self.u, self.s, self.R301, self.R302
        product, tea = self.product_stream, self.tea
        stream_price = self._solve_stream_price()
        F_LA = float(product.imass['LacticAcid'])
        purity = F_LA / product.F_mass
        broth = R301.outs[0]
        F_ferm = float(2 * broth.imol['CalciumLactate'] * 90.07794 + broth.imass['LacticAcid'])
        F_sugar = self._sugar_fed()
        tau = float(R301.tau_cofermentation)
        turnaround = float(getattr(R301, 'tau_batch_turnaround', getattr(R301, 'tau_turnaround', 0.)))
        heat_utilities = [hu for i in self._process_units for hu in i.heat_utilities]
        pu = self.system.power_utility
        y, t, r = self._target
        X = R301.cofermentation_rxns.X
        out = dict(
            MPSP=stream_price / purity,
            recovery=F_LA / F_ferm if F_ferm else float('nan'),
            F_final_product=F_LA,
            F_ferm_product=F_ferm,
            F_sugar_fed=F_sugar,
            F_feedstock=float(self.feedstock_stream.F_mass),
            feedstock_price=float(self.feedstock_stream.price),
            V_broth=float(broth.F_vol),
            F_water_broth=float(broth.imass['Water']),
            V_fermenters=float(broth.F_vol * (tau + turnaround)),
            tau_fermentation=tau,
            base_cost=float(self.base_stream.cost),
            acid_cost=float(self.acid_stream.cost),
            heating_duty=ptu.get_heating_duty(heat_utilities),
            cooling_duty=ptu.get_cooling_duty(heat_utilities),
            electricity_consumption=pu.consumption / 1e3,
            electricity_production=pu.production / 1e3,
            natural_gas_flow=float(self.natural_gas.F_mass),
            electricity_price=float(bst.PowerUtility.price),
            natural_gas_price=float(self.natural_gas.price),
        )
        try:
            gwp, fec = float(self.funcs['get_GWP']()), float(self.funcs['get_FEC']())
            out['GWP'], out['FEC'] = gwp / purity, fec / purity  # per kg pure LA
            out['x_GWP_per_kg_product_stream'] = gwp
            out['x_FEC_per_kg_product_stream'] = fec
        except Exception:
            pass
        inflection = 1. - float(R301._X[1]) - float(R301._X[2])
        fs = describe_feedstock(self.feedstock_stream)
        sacch_cost = 0.
        if self.config['saccharification'] == 'SSCF':
            ic = R301.installed_costs
            sacch_cost = float(sum(v for k, v in ic.items() if 'accharification' in k))
        x = dict(
            x_titer_achieved=float(R301.effluent_titer),
            x_yield_measured=F_ferm / F_sugar if F_sugar else float('nan'),
            x_lactic_conversion=float(X[0]),
            x_acetic_conversion=float(X[1]),
            x_cell_mass_conversion=float(X[2]),
            x_inflection_yield=inflection,
            x_above_inflection=bool(y > inflection + 1e-9),
            x_max_titer_attainable=float(self._attainable_titer(y)),
            x_product_purity=purity,
            x_product_stream_price=stream_price,
            x_product_stream_flow=float(product.F_mass),
            x_cell_mass=float(broth.imass['FermMicrobe']),
            x_acetic_acid_broth=float(broth.imass['AceticAcid']
                                      + 2 * broth.imol['CalciumAcetate'] * 60.05196),
            x_residual_sugar_broth=float(broth.imass[FERMENTABLE].sum()),
            x_lime_fermentation=float(self.base_stream.F_mass),
            x_sulfuric_acid_acidulation=float(self.acid_stream.F_mass),
            x_gypsum=float(s.gypsum.F_mass),
            x_CSL=float(s.CSL.F_mass),
            x_ethanol_makeup=float(s.ethanol.F_mass),
            x_ethyl_acetate_recycle=float(u.R402.ins[1].imass['EthylAcetate']),
            x_dilution_water=float(R301.ins[4].F_mass),
            x_sugar_conc_fermentation=F_sugar / float(broth.F_vol) if broth.F_vol else float('nan'),
            x_max_sugar_batch=float(getattr(R301, 'max_sugar', float('nan'))),
            x_F401_water_evaporated=float(u.F401.outs[0].imass['Water']),
            x_R301_saccharification_installed_cost=sacch_cost,
            x_natural_gas_cost=float(self.natural_gas.cost),
            x_BT_steam_heating_duty=ptu.get_heating_duty(self.BT.heat_utilities),
            x_operating_days=float(tea.operating_days),
            x_feedstock_carbohydrate_content=self._get_carbs(),
            x_fs_w_glucan=fs.get('fs_w_glucan', float('nan')),
            x_fs_w_xylan=fs.get('fs_w_xylan', float('nan')),
            x_fs_w_fermentable_sugar_eq=fs.get('fs_w_fermentable_sugar_eq', float('nan')),
            x_baseline_MPSP_module_default=self.baseline_MPSP_stream,
            x_n_simulations=float(getattr(self, '_n_sims', float('nan'))),
            x_stability_rel_change=float(getattr(self, '_stability', float('nan'))),
        )
        if self.config['saccharification'] == 'SHF':
            x['x_sugar_evaporator_V'] = float(u.E301.V)
            x['x_sugar_evaporator_water_removed'] = float(u.E301.outs[1].F_mass)
        out.update(x)
        return out
