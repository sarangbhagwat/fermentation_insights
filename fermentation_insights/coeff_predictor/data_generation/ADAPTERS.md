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
| `base_cost`, `acid_cost` | USD/h of neutralization base and of acid used for acidulation/regeneration (0 if none) |
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
    neutralizing_agent='Ca(OH)2',       # '' if none
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
)
```

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
`baseline_productivity`, `baseline_yield`, `baseline_titer`. Use the same ranges
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
