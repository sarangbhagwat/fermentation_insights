# Writing a biorefinery adapter

An adapter connects one Bioindustrial-Park biorefinery to the generic
data-generation harness (`harness.py`). The harness sweeps a small titer–yield
(TY) design for many *samples* (configuration × uncertain-parameter values ×
productivity) and records, at every TY point, the MPSP plus system-level
(TEA, per-process-area costs/utilities) and process-level (mass and energy
flows) data. These intermediate data are what the physics-informed model uses,
so **correctness and consistency of definitions across adapters matter more
than anything else**.

Adapters live in `adapters/<NAME>.py` and define `class Adapter(BiorefineryAdapter)`.
Run one with:

```bash
PYTHONPATH=<bip_worktree>:/home/user/fermentation_insights \
  <env>/bin/python -m fermentation_insights.coeff_predictor.data_generation.run \
  --adapter <NAME> --out <dir> --n-random 2 --n-productivity 1 --sample-limit 2 --tag _smoke
```

Do **not** edit `harness.py`, `schema.py` or `run.py` (they are shared). If you
find a bug there, work around it in your adapter and report it.

## Required semantics

### Units and definitions (identical across all adapters)
| quantity | definition |
|---|---|
| `y` (yield) | **absolute** g fermentation product / g fermentable sugar (glucose-eq. mass basis as used by the biorefinery's own spec; for salts, acid-eq. mass). The harness passes absolute yields; convert to the biorefinery's own convention (e.g. % of theoretical) inside `set_try`. |
| `t` (titer) | g fermentation product (acid-eq.) / L broth |
| `r` (productivity) | g / L / h |
| `MPSP` | USD / kg **pure** final product, computed exactly as in the biorefinery's TRY-analysis scripts (e.g. `tea.solve_price(product) / purity`). |
| `recovery` (g) | kg pure final product / kg fermentation product (acid-eq.) leaving fermentation, as in the TRY scripts' recovery metric. For upgraded products this includes the stoichiometric mass change (e.g. acrylic acid / 3-HP). |
| `F_final_product` | kg/h pure final product |
| `F_ferm_product` | kg/h fermentation product in broth (acid-eq. mass; salts converted) |
| `F_sugar_fed` | kg/h fermentable sugars (glucose-eq. mass) entering fermentation |
| `F_feedstock`, `feedstock_price` | kg/h wet feedstock and USD/kg wet |
| `V_broth`, `F_water_broth` | m3/h and kg/h water of the broth leaving fermentation |
| `V_fermenters` | m3 total volume of all main fermentation vessels |
| `tau_fermentation` | h |
| `base_cost`, `acid_cost` | USD/h of base dosed in fermentation (neutralization *or* pH control, at its purchase price: the purchased feed, not an intermediate stream) and of acid used for acidulation/regeneration (0 if none). With `neutralization=False` `base_cost` can still be nonzero: TAL (NaOH pH control), succinic low-pH (partial lime pH control), HP low-pH (medium lime bought although it never enters the broth); `neutralizing_agent` names that base. With `neutralization=True` the base must actually convert the product to its salt and the salt must be acidulated downstream where the biorefinery does so (succinic v3: lime -> calcium succinate -> H2SO4 acidulation + gypsum filtration) |
| `heating_duty`, `cooling_duty` | GJ/h, sums over process units (exclude utility-producing facilities such as the boiler) |
| `electricity_consumption`, `electricity_production` | MW |
| `natural_gas_flow` | kg/h purchased natural gas (0 if none) |
| `electricity_price`, `natural_gas_price` | USD/kWh, USD/kg (current values) |
| `GWP`, `FEC` | optional; per kg final product, if the biorefinery has an LCA object (cheap to read) |

Also return any useful adapter-specific quantities with prefix `x_` (e.g.
`x_evaporator_water_removed` kg/h, `x_CO2_vent` kg/h, `x_cell_mass` kg/h).

### `configurations()`
Return one dict per design configuration. Required keys:

```python
dict(
    config_id='HP_neutral_cornstover',  # unique, stable
    feedstock='cornstover',             # glucose | corn | sugarcane | cornstover | <other>
    feedstock_type='lignocellulose',    # pure_sugar | starch | sugar_juice | lignocellulose
    fermentation_product='HP',          # short label
    final_product='AcrylicAcid',        # short label
    ferm_product_chemical='HP',         # thermosteam chemical ID in the system
    final_product_chemical='AcrylicAcid',
    neutralization=True,                # pH control by base addition during fermentation
    neutralizing_agent='Ca(OH)2',       # base dosed in fermentation (also pH control); '' if none
    n_acid_groups=1,                    # carboxylic acid groups on the fermentation product
    pKa1=4.51,                          # nan if not an acid
    final_product_stoich_mass_ratio=0.80,  # theoretical kg final / kg fermentation product
    upgrading=True,                     # chemical/catalytic conversion after fermentation
    separation_variant='',              # e.g. 'hexanol_extraction'; '' = biorefinery default
    separation_principle='evaporation+crystallization',  # free text but consistent
    aerobic=True,
    fermentation_mode='batch',          # batch | fed-batch | continuous
    T_fermentation=305.15,              # K
    product_water_solubility=1000.,     # g/L at 25 C (literature; use 1000 if miscible)
    adapter_version=2,                  # bump when a fix changes recorded semantics
)
```

**Versioning.** A campaign runs each configuration in its own process, so an
edited adapter takes effect at the next configuration a campaign starts.
Whenever a fix changes what a recorded column means or holds, add or bump
`adapter_version` in every configuration dict of the adapter. It is written
into the `config` record. Then add the matching correction for older rows
to `coeff_predictor/data_patches.py` if one can be derived from recorded
columns (DESIGN.md, "Data corrections"). Harness changes bump
`schema.HARNESS_VERSION`, which is written into each `sample_start`.
Changes so far:

* harness v3: the boiler's natural gas (`BoilerTurbogenerator`'s stream
  utility 'Fuel') and any inlet utility of >= 50 wt% CH4 are now natural gas
  in the utility pool (`ng_*_utility`). v2 booked the boiler gas as unpriced
  material (D3).
* harness v4: every campaign applies the BioSTEAM runtime corrections of
  `harness.BIOSTEAM_FIXES` (`biosteam_fixes.py`; so far `'ic_biogas'`:
  BioSTEAM's `InternalCirculationRx` with the default `method='lumped'` no
  longer counts its biogas twice, D17) in its own process before any system
  is loaded, and writes the active ones into the `config` record and every
  `sample_start` (`biosteam_fixes`, ';'-joined). An adapter inherits them
  (`BiorefineryAdapter.biosteam_fixes`); only an adapter without BioSTEAM
  units may set it to `()` (the synthetic test adapters). Rows of the 27
  configurations with such a reactor that do not name `ic_biogas` (harness
  <= 3) cannot be patched and are rejected (`unpatchable_D17`); the
  configurations were regenerated (2026-10-18). The other 14 (lactic,
  succinic, ethanol_corn, ethanol_cornstover) have no such reactor
  (`sim_checks.py` check `ic_biogas`) and keep their harness-v3 rows.
* HP v2: `base_cost` is the purchased `lime_fresh` at the sampled
  `p_base_price`. There is now a per-point `x_inflection_yield` = (1 − sampled
  cell-mass conversion − succinic − xylitol conversions) × theoretical, plus
  `x_above_inflection`. `separation_principle` starts with 'acidulation+' only
  when neutral.
* TAL v2: `p_cell_mass_yield` is mapped (the canonical names now match
  case-insensitively), `neutralizing_agent='NaOH'`, and
  `product_water_solubility` is 9.8 g/L (experimental data).
* succinic v2: `neutralizing_agent='Ca(OH)2'` also in the low-pH configurations.
* succinic v3 (neutral configurations only; the low-pH ones stay v2): the
  fermenter's lime neutralizes the product (`R302.base_neutralizes_product
  = True`): calcium succinate (and acetate) in the broth, acidulation with
  H2SO4 in R401 (`acid_cost` = the purchased 93 wt% `sulfuric_acid_acidulation`)
  and gypsum filtration in S405 (stream `gypsum`, at the sampled gypsum
  price, 0 at baseline; a negative price is a disposal cost and lands in
  `coproduct_sales`). `F_ferm_product`, the titer and the recovery stay on
  the acid-equivalent basis (calcium succinate counted as succinic acid).
  New per-point columns `x_calcium_succinate_broth`, `x_gypsum_flow`;
  `separation_principle` starts with 'acidulation(H2SO4)+gypsum_filtration+'.
  v2 only set `spec.neutralization`, so the purchased lime was discarded and
  R401 / R401_P / S405 were bypassed (DESIGN.md D15): such rows are rejected
  (`unpatchable_D15`), and the four configurations were regenerated.
* TAL v3 (TAL_SA_sugarcane only; the other TAL configurations stay v2):
  Bioindustrial-Park fi-TAL feeds the juice recycle `M201-0` to the crushing
  mill `U201` once (the upstream system appended it a second time: D16), and
  the configuration's max-titer line is recalibrated on the corrected
  system (679 y − 1.3 g/L instead of 614 y − 0.6; TAL_sugarcane: 678 y −
  1.4). Older TAL_SA_sugarcane rows cannot be patched and are rejected
  (`unpatchable_D16`); the configuration was regenerated (2026-10-18).
* TAL, no version change (recorded semantics unchanged): the root search of
  every MultiEffectEvaporator is robust (`TAL._robust_V_overall`). The KSA
  systems' IPA evaporator F407 used to fail ("Failed to extrapolate enthalpy
  of vaporization ... at T=1207 K") whenever a trial first-effect fraction
  made one effect evaporate all the IPA left (from a guess far from the new
  solution, or from the bounds once `reset()` had cleared the guess, so the
  reset and the reload failed the same way; DESIGN.md H1). The overall
  fraction there is exactly 1, which the wrapper now returns; converged
  points are unchanged.
* ethanol v2: the corn yeast growth now conserves mass (BIP fi-ethanol).
  Corn `natural_gas_flow`/`natural_gas_price` cover all natural-gas
  utilities at the price BioSTEAM charges. There is a new corn
  `p_natural_gas_price` (`stream_prices['Natural gas']`, ±10 %), and the
  no-effect heating-agent price parameter is gone. Sugarcane RIN credits are
  zeroed, with the foregone credit recorded as `x_RIN_credit` [USD/h], so the
  MPSP is no longer net of policy credits.

### `load(config)`
Build/load the system; set `self.system`, `self.tea`, `self.unit_groups`
(list of `biosteam.UnitGroup`), `self.feedstock_stream`, `self.product_stream`
(the main final-product stream whose price is solved for; the harness uses it
to exclude the product's own sales from `coproduct_sales`), and anything else
you need. The harness calls it once per configuration (configurations can also be
run in separate processes via `--configs`).

### `domain()`
Return `theoretical_yield`, `inflection_yield` (yield above which the
cell-mass/byproduct regime changes; the fitted region ends here), `y_min`,
`y_max` (= inflection yield, absolute g/g), `t_min`, `t_max` (g/L),
`baseline_productivity`, `baseline_yield`, `baseline_titer`. If the inflection
depends on a sampled parameter (e.g. a cell-mass yield), also return
`x_inflection_yield` (the sample's own) and `x_above_inflection`
(`y > x_inflection_yield`) from `collect()`. `dataset.build` drops flagged
points and uses the per-sample inflection. Use the same ranges
as the biorefinery's own TRY-analysis scripts / the fermentation_insights
publication (e.g. HP: yields 5–75.5% of theoretical, titers 10–200 g/L).

### `max_titer(y)` (strongly recommended)
Analytical upper bound on the feasible titer at yield `y` (e.g. from the
maximum attainable sugar concentration: `t_max_feasible ≈ y × S_max`). Points
above it waste simulations.

### `parameters()`
Uncertain technological/contextual parameters for the loaded configuration,
taken from the biorefinery's own (latest, default) parameter-distribution
files/models, **excluding yield, titer and productivity** (the harness sets
those). Use `Parameter.from_biosteam(p, name=...)` and give canonical names
(`p_feedstock_price`, `p_electricity_price`, `p_natural_gas_price`,
`p_IRR`, `p_income_tax`, `p_operating_days`, `p_feedstock_capacity`,
`p_boiler_efficiency`, `p_turbogenerator_efficiency`, `p_cell_mass_yield`,
`p_base_price`, `p_acid_price`, `p_enzyme_loading`, `p_saccharification_yield`)
where they apply; otherwise the default `bp_<slug>` name. Keep units as in the
biorefinery but document them in the `Parameter.description`. If a setter
needs a re-simulation to take effect, that is fine — the harness simulates
after setting TRY.

### `set_try(y, t, r)` and `simulate()`
Set the fermentation specifications and simulate robustly (reuse the
biorefinery's own "bugfix barrage" logic: reset cache / empty recycles /
switch solver). Raise `InfeasibleSpecification` when the specification is
physically infeasible (e.g. required sugar concentration not attainable,
product-inhibition limits), and let other exceptions propagate.

### `collect()`
Return the per-point quantities above (the harness adds TEA metrics and
per-area metrics itself from `self.tea`, `self.unit_groups` and
`self.area_map`). `area_map` maps every UnitGroup name to one of
`schema.CANONICAL_AREAS`
(`feedstock, conversion, fermentation, separation, upgrading, wastewater, chp,
facilities`); map a group to `None` only if it would double count.

## Verification checklist (report the results)
1. Baseline MPSP of each configuration vs. the value reported by the
   biorefinery (README/TRY scripts) — explain differences (new BioSTEAM).
2. At 2–3 TY points: measured yield `F_ferm_product/F_sugar_fed` ≈ `y`;
   titer `F_ferm_product/V_broth` ≈ `t`; recovery in a plausible range.
3. Sum over areas of `installed_cost` ≈ TEA `installed_cost`.
4. Typical seconds per TY point and failure rate.
5. A smoke campaign (`--sample-limit 2`) completes and the JSONL looks right.
6. Mass balances: compare every unit's outlets with its inlets and explain
   each unit that creates or destroys mass (a module quirk that is
   cost-neutral, or a defect); no unit lists the same inlet stream twice
   (`sim_checks.py` check `duplicate_inlets`); the sugar fed to
   fermentation does not exceed the sugar the feedstock can supply
   (`F_sugar_fed <= F_feedstock x fs_w_fermentable_sugar_eq x sugar
   availability`, `report_checks.sugar_balance_check`, REPORT.md section 1).
   These checks found defects in the shipped data (D16, D17 below). Not
   every unit closes in the current data: the remaining non-closures are
   listed in DESIGN.md section 1 ("Remaining unit-level mass
   non-closures"), judged cost-neutral.
7. Wastewater treatment with BioSTEAM's `InternalCirculationRx`: its biogas
   is counted twice (D17); check the IC unit's mass balance, and simulate
   with `biosteam_fixes.apply('ic_biogas')` to see the effect.
8. Unit design results against the unit classes' own bounds (`_bounds`;
   BioSTEAM does not warn), e.g. distillation columns at 3-24 ft diameter:
   `sim_checks.py` check `unit_bounds` (D18: the HP purification columns are
   27-40 ft wide).

Status of item 1: the comparison of the baseline coefficients of the 32
configurations of the publication with the publication's (REPORT.md section
8 (i); known causes of the differences: DESIGN.md section 1, "Comparison
with the publication") REPLACES the per-configuration comparison with the
module READMEs, which has not been done: the README values are single
points at module-specific TRY and price settings and older BioSTEAM
versions, and the module defaults under the current BioSTEAM (2.54,
Bioindustrial-Park b3988834) are recorded by two adapters only. What is
recorded (module default under current BioSTEAM, USD/kg of the product
stream, vs. the module README):

| configuration | module default now | README | ratio | known cause |
|---|---|---|---|---|
| lactic_neutral_SSCF_cornstover | 1.5220 (`x_baseline_MPSP_module_default`) | 1.414 (`lactic/README.rst:24`, SSCF) | 1.076 | not reconciled (BioSTEAM / module version) |
| lactic_neutral_SHF_cornstover | 1.5669 | 1.454 (`lactic/README.rst:35`, SHF) | 1.078 | not reconciled |
| lactic_lowpH_SSCF_cornstover | 1.1737 | – | – | – |
| lactic_lowpH_SHF_cornstover | 1.2210 | – | – | – |
| ethanol_corn / ethanol_cornstover / ethanol_sugarcane | 0.5650 / 0.7079 / 0.7786 (`x_module_baseline_MPSP`) | – | – | – |

The recorded lactic baseline points equal these values divided by the
product purity (e.g. 1.7287 USD/kg pure lactic acid for neutral SSCF). The
other module READMEs (`HP/README.md:46` 1.500 USD/kg acrylic acid,
`TAL/README.md:59` 4.869, `succinic/README.rst:38` 1.367,
`oxalic/README.rst:35` 1.77) are at TRY points and versions that are not
reconciled with the adapters' baselines.
For items 2-4, every point records the measured flows (`F_ferm_product`,
`F_sugar_fed`, `V_broth`, `F_final_product`) and the recovery, `dataset.py`
rejects samples whose cost-channel decomposition does not close
(`closure_rel`, `high_closure`), and the failure and rejection counts are in
`report.json` and REPORT.md section 1. The per-point check of y and t
against the measured flows (F_ferm_product / F_sugar_fed / y and
F_ferm_product / V_broth / t, min / median / max per configuration) is
reported in REPORT.md section 1 (`report_checks.yt_definition_check`,
`results/tables/yt_definition_check.csv`; report only, it rejects nothing).

## Module and BioSTEAM defects corrected by the regeneration of 2026-10-18

Found in the final data review (2026-10-18) and documented in DESIGN.md
section 1, "Known data defects"; both were in the training data of the v3
models and are corrected in the current data (harness v4, TAL adapter v3;
the 27 affected configurations were regenerated, old raw files in
`raw_archive/D16_D17_v3`). `sim_checks.py` (`scripts/run_sim_checks.sh`)
re-simulates the configurations with and without the corrections and
writes the effect to `checks/` (`checks/README.md`).

* **D16, TAL_SA_sugarcane.** Bioindustrial-Park
  `TAL/systems/system_SA_solubility_exploit_ethanol_sugarcane.py:116`
  (upstream, not a change of the fi-TAL port) appended the juice recycle
  `M201-0` to the crushing mill `U201`, which BioSTEAM 2.54's juicing system
  already connects, so the recycle entered the mill twice and the
  fermentation received 15 % more sugar per kg cane than the other
  sugarcane configurations (more than the cane holds). Fixed on fi-TAL (the
  line is commented out, as in the other six TAL systems; the unused
  `system_SA_THF_Ethanol_solubility_exploit_ethanol_sugarcane.py` keeps the
  upstream line), TAL `adapter_version` 3 for the configuration, max-titer
  line recalibrated. `sim_checks.py` check `duplicate_inlets` finds no unit
  with a repeated inlet in the 41 configurations any more
  (`checks/duplicate_inlets_bias.csv`; the finding before the fix:
  `checks/duplicate_inlets_v3_*.csv`).
* **D17, every system with an `InternalCirculationRx`** (all HP and TAL,
  ethanol_sugarcane, IBO_EtOH_corn, oxalic_neutral_sugarcane). BioSTEAM's
  lumped IC reactor model (`biosteam/wastewater/high_rate/internal_circulation_rx.py`,
  `_run`) degasses both the reacted helper copy of the influent and the
  reacted outlets into the biogas, so the biogas, and the boiler fuel, was
  counted twice. The runtime correction `biosteam_fixes.py` ('ic_biogas':
  only the helper copy's gas is discarded; BioSTEAM itself is untouched) is
  applied by the harness since v4. `sim_checks.py` check `ic_biogas`
  simulates each affected configuration as harness <= 3 did (`as_is`), as
  the current harness does (`harness`) and with the correction (`fixed`):
  `harness` = `fixed` (`checks/ic_biogas_*.csv`). The defect should be
  reported to BioSTEAM.
