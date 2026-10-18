# -*- coding: utf-8 -*-
"""Forward fixes that can be checked without BioSTEAM: harness natural-gas
pools (D3, harness v3), adapter configuration descriptors and versions."""
import importlib

import pytest

from fermentation_insights.coeff_predictor.data_generation import harness, schema


class _Chems(tuple):
    pass


class _S:
    """Duck-typed stream: F_mass [kg/h], price [USD/kg], composition."""
    def __init__(self, price=0., **mass):
        self.imass = dict(mass)
        self.chemicals = _Chems(mass)
        self.F_mass = float(sum(mass.values()))
        self.price = price

    @property
    def cost(self):
        return self.F_mass * self.price


class _U:
    def __init__(self, utilities):
        self.ins = [s for _, s in utilities]
        self._inlet_utility_indices = {name: i for i, (name, _) in enumerate(utilities)}


class _Sys:
    def __init__(self, units, feeds):
        self.units = units
        self.feeds = feeds


def test_harness_v3_natural_gas_pools():
    assert schema.HARNESS_VERSION >= 3
    fuel = _S(CH4=100.)                         # BoilerTurbogenerator 'Fuel' utility, unpriced
    dryer = _S(CH4=10., price=0.4)              # 'Natural gas' utility with its own price
    biogas = _S(CH4=6., CO2=4.)                 # >= 50 wt% CH4 utility of another name
    ash = _S(Ash=5.)                            # other stream utility
    free_ch4 = _S(CH4=7., price=0.2)            # unregistered CH4 feed
    units = [_U([('Fuel', fuel), ('Ash disposal', ash)]), _U([('Natural gas', dryer)]),
             _U([('Biogas utility', biogas)])]
    sys = _Sys(units, feeds=[fuel, dryer, biogas, free_ch4])
    prices = {'Fuel': 0.25, 'Natural gas': 0.3, 'Biogas utility': 0.1, 'Ash disposal': -0.01}
    d = harness.collect_natural_gas(sys, prices=prices)
    assert d['ng_flow_utility'] == pytest.approx(110.)            # fuel + biogas
    assert d['ng_cost_utility'] == pytest.approx(100 * 0.25 + 10 * 0.1)
    assert d['ng_flow_material'] == pytest.approx(17.)           # priced dryer gas + free feed
    assert d['ng_cost_material'] == pytest.approx(10 * 0.4 + 7 * 0.2)
    assert harness._is_natural_gas_utility('fuel', _S(Water=1.))
    assert not harness._is_natural_gas_utility('Ash disposal', ash)


def test_harness_v4_applies_and_records_the_biosteam_fixes(tmp_path, monkeypatch):
    """D17 (harness v4): run_campaign applies the adapter's BioSTEAM corrections
    (default harness.BIOSTEAM_FIXES) before loading and names them in every
    config and sample_start record; the synthetic adapters declare none."""
    import json
    import types
    from fermentation_insights.coeff_predictor import synthetic
    from fermentation_insights.coeff_predictor.data_generation import biosteam_fixes
    assert schema.HARNESS_VERSION == 4 and harness.BIOSTEAM_FIXES == ('ic_biogas',)
    assert harness.BiorefineryAdapter.biosteam_fixes == harness.BIOSTEAM_FIXES
    assert synthetic.SyntheticAdapter.biosteam_fixes == ()

    class IC:  # stand-in for biosteam.wastewater.high_rate.InternalCirculationRx
        method = 'lumped'

        def _run(self):
            return 'original'
    original = IC._run
    monkeypatch.setattr(synthetic._BST, 'wastewater',
                        types.SimpleNamespace(high_rate=types.SimpleNamespace(InternalCirculationRx=IC)),
                        raising=False)
    ad = synthetic.SyntheticAdapter('SYNA')
    ad.biosteam_fixes = harness.BIOSTEAM_FIXES          # as a BioSTEAM adapter
    with synthetic.biosteam_stub():
        try:
            path = harness.run_campaign(ad, str(tmp_path), n_random=1, n_productivity=0, seed=1,
                                        n_y=3, n_t=3, point_timeout=0, verbose=False)
            assert IC._run is biosteam_fixes._ic_run_degas_once
            assert harness.active_biosteam_fixes() == 'ic_biogas'
        finally:
            biosteam_fixes.remove()
        assert IC._run is original and harness.active_biosteam_fixes() == ''
    recs = [json.loads(line) for line in open(path)]
    starts = [r for r in recs if r['record'] == 'sample_start']
    assert len(starts) >= 2 and all(r['harness_version'] == 4 and r['biosteam_fixes'] == 'ic_biogas'
                                    for r in starts)
    assert all(r['biosteam_fixes'] == 'ic_biogas' for r in recs if r['record'] == 'config')
    # the synthetic adapters (no BioSTEAM units) record no correction
    with synthetic.biosteam_stub():
        path = harness.run_campaign(synthetic.SyntheticAdapter('SYNB'), str(tmp_path), n_random=1,
                                    n_productivity=0, seed=1, n_y=3, n_t=3, point_timeout=0, verbose=False)
    assert all(json.loads(line).get('biosteam_fixes', '') == '' for line in open(path))


def _configs(name):
    try:
        m = importlib.import_module(f'fermentation_insights.coeff_predictor.data_generation.adapters.{name}')
    except ImportError as e:  # pragma: no cover - environment without an adapter dependency
        pytest.skip(str(e))
    return m, {c['config_id']: c for c in m.Adapter().configurations()}


def test_adapter_versions_and_descriptors():
    _, hp = _configs('HP')
    assert all(c['adapter_version'] == 2 for c in hp.values())
    for c in hp.values():
        assert c['separation_principle'].startswith('acidulation+') == c['neutralization']
    tal_m, tal = _configs('TAL')
    for cid, c in tal.items():
        # v3 (D16): only TAL_SA_sugarcane changed (juice recycle fed to the mill once)
        assert c['adapter_version'] == (3 if cid == 'TAL_SA_sugarcane' else 2)
        assert c['neutralizing_agent'] == 'NaOH' and not c['neutralization']
        assert c['product_water_solubility'] == pytest.approx(9.8)
    # BioSTEAM capitalizes only the first character of parameter names (M3)
    assert tal_m._CANONICAL_LC['fermentation y. lipolytica yield'] == 'p_cell_mass_yield'
    _, su = _configs('succinic')
    for c in su.values():
        neutral = c['config_id'].startswith('succinic_neutral_')
        assert c['neutralization'] == neutral and c['neutralizing_agent'] == 'Ca(OH)2'
        # v3 (D15): only the neutral configurations changed (salt formation, acidulation)
        assert c['adapter_version'] == (3 if neutral else 2)
        assert c['separation_principle'].startswith('acidulation(H2SO4)+gypsum_filtration+') == neutral
        assert c['separation_principle'].endswith('evaporation+crystallization')
    assert not su['succinic_corn']['neutralization'] and su['succinic_neutral_corn']['neutralization']
    _, et = _configs('ethanol')
    assert all(c['adapter_version'] == 2 for c in et.values())


def test_succinic_neutral_lime_forms_the_salt(monkeypatch):
    """D15 (adapter v3): loading a neutral configuration makes the fermenter's
    lime neutralize the product (calcium succinate), which switches on the
    acidulation R401 / R401_P and the gypsum filter S405 (their bypass
    specifications read ``R302.base_neutralizes_product``); the low-pH
    configurations keep the module default (no salt, R401/S405 bypassed),
    also when switched within one process."""
    from types import SimpleNamespace as NS
    from fermentation_insights.coeff_predictor import synthetic
    m, su = _configs('succinic')
    ad = m.Adapter()
    ad._feedstock_loaded = 'corn'           # skip the Bioindustrial-Park import
    ad.spec = NS(neutralization=None)
    ad.R302 = NS(base_neutralizes_product=None, neutralization=None, T=303.15)
    ad.model = NS(metrics_at_baseline=lambda: None)
    ad.system = NS(converge_method=None)
    ad._converge_method = 'wegstein'
    monkeypatch.setattr(ad, '_make_unit_groups', lambda bst: [])
    with synthetic.biosteam_stub():
        ad.load(su['succinic_neutral_corn'])
        assert ad.R302.base_neutralizes_product is True and ad.spec.neutralization is True
        ad.load(su['succinic_corn'])
        assert ad.R302.base_neutralizes_product is False and ad.spec.neutralization is False
        ad.load(su['succinic_neutral_corn'])
        assert ad.R302.base_neutralizes_product is True
        with pytest.raises(RuntimeError, match='one feedstock'):
            ad.load(su['succinic_neutral_glucose'])


def test_tal_multi_effect_evaporator_root_search_is_robust():
    """H1: in the region where an effect evaporates all of the chemical left
    (superheated vapor; the next effect fails evaluating its enthalpy of
    vaporization above Tc), the overall evaporated fraction is exactly 1;
    other errors still propagate."""
    m, _ = _configs('TAL')
    msg = "Failed to extrapolate enthalpy of vaporization method 'DIPPR_PERRY_8E' at T=1207.2 K"

    def V_overall(v):  # monotonic; crashes inside the all-evaporated region (v > 0.15)
        if 0.2 < v < 0.9:
            raise RuntimeError(msg)
        return min(1., 0.3728 + 4.05 * v)

    f = m._robust_V_overall(V_overall)
    assert f(0.05) == pytest.approx(0.5753)
    assert f(0.5) == 1. and f(1.) == 1.
    # bisection on the wrapped objective finds the root left of the crash region
    lo, hi = 0., 1.
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) < 0.65 else (lo, mid)
    assert lo == pytest.approx((0.65 - 0.3728) / 4.05, abs=1e-9)

    def other(v):
        raise RuntimeError('could not converge')
    with pytest.raises(RuntimeError, match='could not converge'):
        m._robust_V_overall(other)(0.5)
