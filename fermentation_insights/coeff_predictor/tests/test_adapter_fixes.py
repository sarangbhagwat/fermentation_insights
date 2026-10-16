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
    assert schema.HARNESS_VERSION == 3
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
    for c in tal.values():
        assert c['adapter_version'] == 2 and c['neutralizing_agent'] == 'NaOH' and not c['neutralization']
        assert c['product_water_solubility'] == pytest.approx(9.8)
    # BioSTEAM capitalizes only the first character of parameter names (M3)
    assert tal_m._CANONICAL_LC['fermentation y. lipolytica yield'] == 'p_cell_mass_yield'
    _, su = _configs('succinic')
    for c in su.values():
        assert c['adapter_version'] == 2 and c['neutralizing_agent'] == 'Ca(OH)2'
    assert not su['succinic_corn']['neutralization'] and su['succinic_neutral_corn']['neutralization']
    _, et = _configs('ethanol')
    assert all(c['adapter_version'] == 2 for c in et.values())
