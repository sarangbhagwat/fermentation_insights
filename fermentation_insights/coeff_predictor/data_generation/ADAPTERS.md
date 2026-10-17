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
