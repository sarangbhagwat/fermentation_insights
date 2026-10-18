# -*- coding: utf-8 -*-
"""
Runtime corrections of BioSTEAM defects found while checking the training
data. Nothing here is applied by the campaigns that built the shipped data
(harness v3): the corrections are opt-in, for the simulation checks in
:mod:`.sim_checks` and for the next data rebuild. BioSTEAM itself is not
modified; the corrections replace a method in memory, in the process that
calls :func:`apply`, until :func:`remove` is called.

``ic_biogas`` (known data defect D17, DESIGN.md section 1)
    ``biosteam.wastewater.high_rate.InternalCirculationRx._run`` (biosteam
    2.54, ``internal_circulation_rx.py``, the default ``method='lumped'``)
    splits the unreacted influent into ``waste`` and ``eff``, runs the
    growth and biogas reactions on the helper copy ``inf`` (used only for the
    effluent COD) and degasses it INTO ``biogas``, then runs the same
    reactions on ``waste`` and ``eff`` and degasses them into ``biogas``
    again. The biogas is therefore counted twice (the unit creates mass:
    +18.5 t/h at the HP_neutral_cornstover baseline) and the boiler burns
    methane that was never produced. The corrected ``_run`` is identical
    except that ``inf`` is degassed into a scratch stream; with
    ``method='separate'`` (no double count of the produced gas) the original
    method runs unchanged.
"""
from __future__ import annotations

__all__ = ('apply', 'remove', 'applied', 'FIXES')

_ORIGINAL = {}


def _ic_class():
    import biosteam as bst
    return bst.wastewater.high_rate.InternalCirculationRx


def _ic_run_degas_once(self):
    """``InternalCirculationRx._run`` of biosteam 2.54 with the biogas of the
    helper copy ``inf`` discarded instead of added to ``biogas``."""
    if str(self.method).lower() == 'separate':
        return _ORIGINAL['ic_biogas'](self)
    import thermosteam as tmo
    inf = self._inf
    inf.copy_like(self.ins[0])
    biogas, eff, waste = self.outs
    degassing = self._degassing
    biogas.phase = 'g'
    biogas.empty()
    inf.split_to(waste, eff, self.q_Qw)
    growth_rxns = self.growth_rxns
    biogas_rxns = self.biogas_rxns
    growth_rxns(inf.mol)
    biogas_rxns(inf.mol)
    # the only change: the helper copy's gas is not part of the outlet
    degassing(inf, tmo.Stream(None, thermo=inf.thermo))
    Se = self.compute_COD(inf)
    Qi, Si, Xi, Qe, Y = self.Qi, self.Si, self.Xi, self.Qe, self.Y_biomass
    Xe = (Qi*Xi + Qi*(Si - Se)*Y) / (Qe + (Qi - Qe)*self.q_Xw)
    Xw = Xe * self.q_Xw
    for rxns in (growth_rxns, biogas_rxns):
        rxns(waste.mol)
        rxns(eff.mol)
    degassing(eff, biogas)
    degassing(waste, biogas)
    eff.imass['WWTsludge'] = Xe*self.Qe
    waste.imass['WWTsludge'] = Xw*self.Qw
    diff = sum(i.imol['WWTsludge'] for i in self.outs) - inf.imol['WWTsludge']
    if diff > 0:
        decay_rxn = self.decay_rxn
        decay_rxn._X = min(1., diff / inf.imol['WWTsludge'])
        for i in (eff, waste):
            decay_rxn.force_reaction(i.mol)
            i.imol['O2'] = max(0, i.imol['O2'])
            degassing(i, biogas)
    if self.T:
        biogas.T = eff.T = waste.T = self.T


#: name -> (class getter, attribute, replacement)
FIXES = {
    'ic_biogas': (_ic_class, '_run', _ic_run_degas_once),
}


def applied(name='ic_biogas'):
    get_cls, attr, new = FIXES[name]
    return getattr(get_cls(), attr) is new


def apply(*names):
    """Apply the named corrections (all if none are named)."""
    for name in names or tuple(FIXES):
        get_cls, attr, new = FIXES[name]
        cls = get_cls()
        if getattr(cls, attr) is not new:
            _ORIGINAL[name] = getattr(cls, attr)
            setattr(cls, attr, new)


def remove(*names):
    """Restore BioSTEAM's original methods."""
    for name in names or tuple(FIXES):
        get_cls, attr, _ = FIXES[name]
        if name in _ORIGINAL:
            setattr(get_cls(), attr, _ORIGINAL.pop(name))
