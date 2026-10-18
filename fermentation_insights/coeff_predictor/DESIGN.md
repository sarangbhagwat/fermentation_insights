# coeff_predictor — design

Goal: predict the four coefficients of the unifying equation (Bhagwat et al.,
*Nat. Commun.* 2026) and the recovery term

    MPSP(y, t) = (a' + b'/y + c'/t + d'/(y·t)) / g

directly from specifications of the **feedstock**, **fermentation product**,
**final product** and **fermentation regime** (plus optional contextual
parameters such as prices), using a physics-informed deep-learning model
trained on BioSTEAM simulations that record intermediate process- and
system-level data.

Conventions (identical to the publication): `y` = g fermentation product per g
fermentable sugar (glucose-eq.), `t` = g/L, `g` = kg final product per kg
fermentation product (median recovery over all simulated, status-ok TY points
of the sample, i.e. before clipping to the fitted domain, as in
`MPSP_YT_fit_only.py`), coefficients are fitted by ordinary least squares to
MPSP over the feasible TY domain (yields up to the inflection yield -- the
sample's own inflection yield, which depends on its sampled byproduct /
cell-mass yields), and `a' = g·A` etc. where `A, B, C, D` are the
MPSP-level coefficients.

## 1. Data (from `data_generation/`)

One *sample* = configuration (feedstock × fermentation product × final product
× regime/separation design) × one draw of the biorefinery's uncertain
parameters × one productivity. Each sample has ~25 simulated TY points with
MPSP, recovery, TEA totals, per-process-area costs/utilities, and mass/energy
flows (see `data_generation/schema.py`).

### Data corrections

A review of the seven adapters found defects in data that had already been
recorded. Forward fixes are versioned, so that every recorded row says
whether it contains them:

* harness fixes by `harness_version` (written into every `sample_start`;
  v3 = D3 below);
* adapter fixes by `adapter_version` (a field of each configuration dict,
  written into the `config` record; 2 = fixed, missing = 1). A sample's
  configuration is the latest `config` record that precedes its
  `sample_start` in its own file. A resumed campaign appends a new record,
  so a file that mixes versions is still handled row by row.

`data_patches.py` corrects rows whose version predates a fix. `dataset.build`
applies it right after reading the raw JSONL (`patches=True`, CLI
`--no-patches`). A patch only changes rows that still show the defect, and it
is idempotent: patched rows and records get the patch id in `data_patches`.
The counts per patch go to `report.json['data_patches']`, and each sample
carries `harness_version`, `adapter_version` and `data_patches`.

| id | affected rows | correction |
|---|---|---|
| D3 | harness < 3: HP, oxalic, isobutanol, TAL, succinic, ethanol cornstover/sugarcane | v2 booked the boiler's natural gas (BioSTEAM `BoilerTurbogenerator` stream utility 'Fuel', which the TEA charges at `stream_prices['Fuel']` = `BT.natural_gas_price`) as an unpriced material feed. The boiler gas `B` moves to the utility pool at `natural_gas_price`: `ng_flow_utility += B`, `ng_cost_utility += B·natural_gas_price`, `ng_flow_material -= B`. `B` = `ng_flow_material` for HP, oxalic and isobutanol (all their material gas is the boiler's), and `x_natural_gas_BT` for the others. Unaffected: lactic (its boiler gas is priced) and ethanol corn (no boiler). |
| D4 | ethanol_sugarcane, adapter < 2 | The MPSP was net of the module's RIN credit, recorded as `coproduct_sales` (0.2593 USD per kg of denatured product, checked to be constant before patching). Patched: `x_RIN_credit = coproduct_sales/h`, `MPSP += coproduct_sales/(h·F_final_product)`, `x_MPSP_denatured`/`x_MESP_…` raised by the credit, `coproduct_sales = 0`. `product_stream_sales` (the harness's sales of the product stream) is not updated and stays net of the credit in the patched rows (58 % off in the ethanol_sugarcane reproduction); no module reads it (`dataset.py`, `physics.py` and the model use `MPSP` and the cost channels), so this is harmless. |
| D6 | HP, adapter < 2 | `base_cost = x_lime_fermentation · p_base_price`. The TEA charges `lime_fresh`, which is linked to `fermentation_lime`, at the sampled price. v1 read `fermentation_lime`, whose price was frozen. |
| D7 | HP, adapter < 2 | `x_inflection_yield = (1 − p_cell_mass_yield)·theoretical_yield`, `x_above_inflection = y > x_inflection_yield`. Basis checked: the config inflection is 0.755 = 1 − 0.245, and `x_cell_mass_conversion` equals `p_cell_mass_yield` at low y. |
| D11 | ethanol_corn, adapter < 2 | The corn natural gas missed the thermal oxidizer X611 and its price was the stream's own (0). Fixed in adapter v2 (sum over the 'Natural gas' and 'Fuel' utilities, `bst.stream_prices`); the affected samples were regenerated with D2. |
| D12 | ethanol_corn, adapter < 2 | The real corn natural-gas price was never sampled (the sampled heating-agent price had no effect). Fixed in adapter v2 (`p_natural_gas_price` sets `stream_prices['Natural gas']`); regenerated with D2. |
| D2 | ethanol_corn, adapter < 2 | The molar yeast growth destroyed sugar mass. This cannot be patched, so the samples are rejected (`unpatchable_D2`) and the configuration was regenerated. |
| D15 | succinic_neutral, adapter < 3 | The salt formation, acidulation and gypsum filtration were never simulated (see "Known data defects"). Not patchable: such samples are rejected (`unpatchable_D15`); the four configurations were regenerated with adapter v3 (2026-10-17). |
| M1 | TAL, succinic low-pH, adapter < 2 | `neutralizing_agent` = 'NaOH' (TAL pH control) or 'Ca(OH)2' (succinic partial pH control). `neutralization` stays False. |
| M3 | TAL, adapter < 2 | `bp_fermentation_fermentation_y_lipolytica_yield` → `p_cell_mass_yield`, in the sample values and in the config's parameter list. |
| M12 | TAL, adapter < 2 | `product_water_solubility` = 9.8 g/L water at 25 °C, interpolated from the TAL module's experimental data. v1 used 5.41 from the fitted model. |
| M23 | HP non-neutral, adapter < 2 | `separation_principle` loses its 'acidulation+' prefix. |
| S1 | every known configuration, `downstream_spec_version` < 1 (all recorded data) | Not a defect fix: adds the structured final-product / downstream-process specification of `downstream.TABLE` (§5) to the config record (the table wins; a differing recorded `coproducts` / `neutralization_degree` is kept as `<field>_recorded`). Versioned by `downstream_spec_version`; `downstream.apply_to_frame` does the same for processed tables built before S1. |

Resumed samples append duplicate rows (M11). `read_raw` already keeps only
the points of each sample's chosen, complete attempt. Points above the
sample's own inflection are dropped (`x_inflection_yield`, `x_above_inflection`;
§7).

### Known data defects

Defects found in the data review that are not (or not fully) corrected by
the patches above. "Effect" is what they do to the training targets and the
model; none of them was corrected for the v2 models. D15 and H1 were fixed
afterwards by regenerating the affected data (2026-10-17); models trained
before that (v2) still carry them. D16 and D17 were found in the final data
review (round 4, 2026-10-18) and are in the data of the shipped (v3)
models. The v3 pipeline trains on the corrected
data and masks the D8-affected feedstock-channel targets in the loss (§6,
round 5).

Disclosure (test split, data verification only): while verifying the D15
regeneration and the D8 identity, the data work of round 5 read recorded
simulation points of a few test-split samples (succinic_glucose|s1|0,
test-OOD configuration; succinic_neutral_glucose|s1|0 and the in-distribution
test samples of TAL_corn, TAL_cornstover and TAL_glucose): to check that an
unchanged configuration reproduces its recorded MPSPs and to read the
deterministic baseline TY grid. No model was fit, tuned, selected or
calibrated with them; the adapter fix came from the Bioindustrial-Park code
and every number quoted here was recomputed on development rows only.

| id | configurations / fields | defect | effect on the coefficients and the model | status |
|---|---|---|---|---|
| D1 | lactic, random samples | the feedstock-price bounds were rebound after sampling (prices outside the intended range) | none now: the random samples were regenerated with the bound fix (`adapters/lactic.py` 408-421; old files in `raw_lactic_D1`, prices 0.0005-0.175 USD/kg = U(0, 200) USD/dry t). The fix was not versioned (lactic configs carry no `adapter_version`) | fixed by regeneration |
| D5 | HP low-pH (8 configurations) | the cation-exchange regeneration H2SO4 (`H2SO4_fresh` → T611) and the anion-exchange NaOH are not in `acid_cost` / `base_cost`; `acid_cost` = 0 in every low-pH HP configuration | their cost is in materials_other instead of the neutralization channel (lactic and oxalic book such reagents as neutralization); totals unaffected, channel labels and the `acid_regeneration` field ('none') understate HP's reagent use | open |
| D8 | TAL, TAL_SA | the spiked acetate (10 g/L AcOH-eq. by default, up to 0.45 g per g sugar) also makes TAL, at the sugars' yield fraction and with a 1.5× higher theoretical yield per g (AcOH → 1/3 TAL, 0.700 g/g, vs glucose 0.4667 g/g), but `F_sugar_fed` counts the sugars only: F_ferm/(y·F_sugar) reaches 1.65, while F_ferm = y·(F_sugar + 1.5·F_acetate) within 0.5 % (median; development rows) | y is the module's (and the publication's) yield on all substrates, glucose-eq. by theoretical TAL yield, so it agrees with ADAPTERS.md; but the sugar demand per kg TAL is 1/y − 1.5·a/t (a = acetate per L broth, 1.8-10.4 g/L between points), so the feedstock channel carries a negative c′ (all 188 development TAL samples; median −3.0, −7.3 to −0.7 USD/kg·g/L; in the intensity, kg feedstock per kg TAL, −98 to −9) and a positive d′ (median 0.33), where the other products have c′ ≈ d′ ≈ 0 (medians 0.000 / 0.002); no sugar-only feedstock description explains that, and, through the feedstock-only sub-network shared by all products of a feedstock, plausibly biases the corn feedstock intensity of other products (REPORT.md section 8 (ii-b) of v2: corn feedstock sensitivities 19-24 % low). Round 5 (§6): those low corn sensitivities came mainly from the coefficient-space channel loss (every corn configuration, ethanol and isobutanol included, 19-27 % low in-sample, 18-34 % out of fold) and vanish with the surface-space channel loss; TAL's large OLS c′ is mostly an ill-conditioned split (a pure b′/y fit is within 3.6-4.6 % RMS of TAL's feedstock intensity, its b′ within ±2-4 % of the product-blind value of the feedstock; TAL_SA_sugarcane −10 %, which is defect D16) | analysed, y kept ("D8 analysis" below); mitigated in the v3 models: the TAL / TAL_SA feedstock-channel targets are masked in the loss (`experiments.DEFECTIVE_TARGETS`, §6 round 5); a carbon-basis regeneration or an acetate co-feed input remain open |
| D9 | lactic | `R403._setup` rebuilds the hydrolysis reaction at X = 0.8 on every simulation, so the sampled `bp_hydrolysis_reactor_r403_hydrolysis_conversion` (0.745-0.850) is never applied | that `bp_*` has no effect (it is not a model input); no bias | open (harmless) |
| D10 | isobutanol | the starch-content setter rescales the current composition, so the result depends on history | sampled and actual dry starch differ (sampled 0.59-0.85, actual 0.67-0.76, corr. 0.87), moisture drifts 13.4-17.1 % (baseline 15 %); the model sees the actual `fs_w_*`, so the inputs are right, but the starch range covered is narrow | open |
| D13 | oxalic | the module's unit groups put lime, H2SO4 and CSL/MgCl2/DAP storage in facilities (99.6 % of the non-feedstock material cost) | channel labels only (facilities instead of fermentation / separation materials); totals unaffected | open |
| D14 | oxalic | the reaction `Glucose -> 2 OxalicAcid + 2 AceticAcid + O2` creates mass and carbon: (OA + AcOH + cells + CO2)/sugar = 1.61-1.70, air feed 5-8 kg/h for 133 t/h sugar, acetic acid 0.667 g/g OA | physically invalid fermentation balance; `aerobic=True` is misleading (no aeration cost), and the oxalic yields are not comparable with the other products' | open (kept for fidelity to the module; patch the reaction and regenerate for physical consistency) |
| D15 | succinic_neutral (4 configurations, 72 samples) | succinic adapter v2 set only `spec.neutralization` (as the publication's TRY script); the reactor kept `base_neutralizes_product=False`, which removes the purchased lime from the broth and bypasses R401 / R401_P / S405 (acidulation, gypsum filtration) | the v2 data equal succinic (low pH) except the base cost (acid, gypsum and acidulation costs 0), but `downstream.TABLE` describes acidulation + gypsum filtration (neutralization degree 1, H2SO4): the acid / gypsum features claimed costs the targets lacked (v2 models). Regenerated data: MPSP +0.22 USD/kg vs. succinic (low pH) at identical TY points and baseline parameters (+0.10 vs. the v2 rows): lime 1 mol/mol SA (+0.5 mol/mol acetic acid) instead of 0.21 (+0.13 USD/kg), 0.89-0.90 kg H2SO4 (93 wt%) per kg SA (+0.09 USD/kg), 1.44-1.45 kg wet gypsum (20 % moisture) per kg SA (0 USD/kg at baseline, sampled −0.029 to +0.008), R401 + S405 +1.0-1.6 M USD installed, while the fermenters are 304 instead of 316 stainless steel (−5 to −9 M USD, in v2 already); recovery unchanged (the module loses no acid in R401 / S405) | fixed by regeneration: adapter v3 (`base_neutralizes_product=True` for the neutral configurations), old files in `raw_archive/succinic_neutral_v2`, v2 rows rejected (`unpatchable_D15`); the v3 models (shipped) were retrained on the corrected data |
| H1 | TAL_SA_glucose, TAL_SA_sugarcane | after a convergence failure whose `reset()` also failed, the old harness wrote the later samples as complete with no simulated point (27 and 19 samples; TAL_SA_glucose's sample 1 was cut after 3 points). Cause: the IPA evaporator F407 of the KSA systems ("Failed to extrapolate enthalpy of vaporization ... at T=1207 K"): BioSTEAM's root search for the first-effect evaporation tried fractions at which an effect evaporates all the IPA left and superheats its vapor; once the search started from the bounds (reset() clears its guess), it failed every time | TAL_SA_glucose kept only its baseline sample, TAL_SA_sugarcane 10 of 29 (v2 models) | fixed: harness reloads after a failed reset and writes `retry` records (resume re-runs them); TAL adapter returns the exact overall fraction 1 in that region (`_robust_V_overall`; re-simulated recorded points change by < 1e-4); the empty and the cut attempts were removed (originals in `raw_archive/TAL_SA_v2`) and re-run: 29 / 29 samples each; F407 failed once more without the mitigation (after 13 re-run samples), never with it; remaining failures 9 of 700 re-run points of TAL_SA_glucose (1.3 %; wastewater recycle M503 not converged at the maximum-titer edge; TAL_glucose 1.0 % in the original campaign), 0 of the 150 points of TAL_SA_sugarcane re-run with the mitigation |
| D16 | TAL_SA_sugarcane (29 samples; 2 of them in the in-distribution test set) | Bioindustrial-Park `TAL/systems/system_SA_solubility_exploit_ethanol_sugarcane.py:116` (upstream master b3988834, not a change of the fi-TAL port) runs `u.U201.ins.append(u.M201-0)`, but in BioSTEAM 2.54's juicing system the juice recycle `M201-0` already feeds the crushing mill: `U201.ins` = [s19, s1, s1], so the recycle enters the mill twice. The six other TAL systems comment this line out (e.g. `system_TAL_solubility_exploit_ethanol_sugarcane.py:124`); a scan of the systems of all 41 configurations (`sim_checks.py` check `duplicate_inlets`) finds no other unit with a repeated inlet | the mill creates mass (+121 t/h at the baseline point, 99 % of the system's +1.38 % imbalance; data review, round 4) and sugar: every simulated point feeds fermentation 1.088 x the sugar the cane can supply (`report_checks.sugar_balance_check`, REPORT.md section 1; every other configuration 0.85-1.00, the other sugarcane ones 0.946), i.e. 15 % more sugar per kg cane than the other sugarcane configurations, so the feedstock demand per kg product is understated. Re-simulated with the duplicate reference removed (in memory; `data_generation/sim_checks.py` check `duplicate_inlets`, `data_generation/checks/duplicate_inlets_bias.csv`), the sugar fed falls by 13 % (81.4 → 70.8 t/h at the baseline point; sugar fed / cane sugar 0.946, as the other sugarcane configurations) and the MPSP rises: +2.1 % at the baseline point, +2.3 to +5.2 % on a 3 x 3 TY grid (+4.8 % at its lowest yield), and +1.3 to +6.3 % (median +3.3 %) on the baseline sample's recorded points (data review). The refitted b′ is 9-12 % low in the data (grid fit 0.299 → 0.326; baseline-sample fit a′ b′ c′ d′ g 2.978 0.286 32.30 0.561 0.704 → 2.894 0.327 33.26 0.404 0.704); the d′ change depends on the fit (+3 % on the grid, −28 % on the baseline sample: an ill-conditioned split). A feedstock-channel proxy puts the MPSP understatement over the 29 samples at about 0.8-8 % (median ~3 %; data review). The corrected b′ gives TAL_SA / TAL = 0.327 / 0.315 = 1.04 on sugarcane, as on corn (1.05), corn stover (1.05) and glucose (1.04); recorded 0.91. This is the TAL_SA_sugarcane b′ deviation (−10 %) of the D8 analysis below, and most of the release model's largest in-sample feedstock-sensitivity error (+19.5 %, REPORT.md section 8 (ii-b)): the product-blind feedstock head predicts about the physical cane demand. The adapter's titer frontier for this configuration (`adapters/TAL.py:152`, 614 y − 0.6 vs. 678 y − 1.4 for TAL_sugarcane) was calibrated on the defective system. The publication's own TAL_SA_sugarcane b′ (0.293 vs. 0.305 for TAL_sugarcane, REPORT.md section 8 (i)) shows the same anomaly, so it very likely carries the same defect | open: the data are not corrected (it would change training data: a pipeline re-run). At the next data rebuild: comment out line 116 on the fi-TAL branch (as in the other TAL systems), bump the TAL `adapter_version`, recalibrate the TAL_SA_sugarcane max-titer line and regenerate its 29 samples (~1 h) |
| D17 | 27 configurations with BioSTEAM's `InternalCirculationRx` (anaerobic IC reactor of the wastewater treatment, default `method='lumped'`): all 16 HP, all 8 TAL, ethanol_sugarcane, IBO_EtOH_corn (small) and oxalic_neutral_sugarcane (negligible); not succinic, lactic, ethanol_corn, ethanol_cornstover | `biosteam/wastewater/high_rate/internal_circulation_rx.py` (`_run`, lines 170-215 of biosteam 2.54; unchanged since the file was added in 2023) splits the unreacted influent into `waste` / `eff`, runs the growth and biogas reactions on the helper copy `inf` (used for the effluent COD) and degasses it INTO `biogas` (line 188), then runs the same reactions on `waste` and `eff` and degasses them into `biogas` again (lines 205-206): the biogas is counted twice. The boiler then burns methane that was never produced | the IC unit creates mass equal to its biogas outlet's second half (+18.5 t/h at the HP_neutral_cornstover baseline; −45 kg/h after the correction), so the steam / power credit is too high and the MPSP too low. Re-simulated with the corrected `_run` (only `inf`'s gas discarded; `data_generation/biosteam_fixes.py`, `sim_checks.py` check `ic_biogas`, every parameter at baseline, the baseline point and a 3 x 3 TY grid; `data_generation/checks/ic_biogas_bias.csv`), the MPSP rises at the baseline point by +1.1-1.8 % (HP glucose / sugarcane), +2.8-3.8 % (HP corn), +4.1-5.4 % (HP corn stover), +0.9 / +1.5 % (TAL_SA glucose / sugarcane), +2.8 / +3.0 % (TAL glucose / sugarcane), +3.4 / +8.3 % (TAL_SA_corn / TAL_corn), +7.2 / +13.3 % (TAL_SA_cornstover / TAL_cornstover), +0.9 % (ethanol_sugarcane), +0.6 % (IBO_EtOH_corn), 0 % (oxalic); on the grid up to +12.9 % (HP corn stover), +15.3 % (TAL_cornstover) and +10.2 % (ethanol_sugarcane, +9.7 % at the lowest yield vs. +0.9 % at baseline). The bias depends on the yield (the wastewater organic load, and with it the phantom gas, grows as the yield falls for corn, corn stover and ethanol_sugarcane), so it does not scale like a constant: the refitted b′ is 10-20 % too low for HP / TAL on corn and corn stover and 16 % for ethanol_sugarcane (≤ 2.4 % for glucose and sugarcane), d′ and a′ change by up to tens of % (ill-conditioned splits), and the electricity / natural-gas intensities that the price heads learn are too favourable (the electricity channel is −23 % of the TAL_cornstover MPSP). The publication's modules very likely ran the same code, so the comparison of section "Comparison with the publication" cannot reveal it | open: the data are not corrected (it would change training data: regeneration of HP, TAL, ethanol_sugarcane and IBO and a pipeline re-run). Predicted MPSPs of these families are optimistic by the amounts above. Prepared for the next rebuild: the opt-in runtime correction `data_generation/biosteam_fixes.apply('ic_biogas')` (BioSTEAM itself untouched, as the TAL `_robust_V_overall` precedent), to be enabled in the harness with a `harness_version` bump; the defect should also be reported to BioSTEAM |

#### D8 analysis (2026-10-17; development rows only)

* **What y is.** The TAL module's yield specification (`spec.load_yield`)
  sets one conversion X (fraction of theoretical) for glucose, xylose and
  acetic acid (`acetate_to_TAL_rxn.X = X`), and the harness passes
  y = X · 0.4667 g TAL/g glucose. The module's theoretical yields are
  0.4667 g/g for the sugars (Glucose → 2/3 TAL + 2 CO2; xylose the same mass
  yield) and 0.700 g/g for acetic acid (AcOH → 1/3 TAL + 3 H2O), so 1 g of
  acetic acid counts as 1.5 g glucose-eq. The recorded rows confirm
  F_ferm_product = y · (F_sugar_fed + 1.5 · x_F_acetate_fed): ratio median
  0.995 (0.975-0.999; the rest is the seed train) over the 3 116 ok points of
  the development samples of TAL_corn, TAL_cornstover, TAL_glucose,
  TAL_SA_corn and TAL_SA_cornstover. y is thus the yield on all substrates,
  glucose-eq. by theoretical TAL yield. That is the module's own definition
  ("weight fraction of substrates converted to product over theoretical
  maximum yield") and the TAL axis of the publication's TY maps, and
  ADAPTERS.md defines y on "the biorefinery's own spec" basis: y is not a
  recording error. The inconsistent pieces are `F_sugar_fed` (sugars only)
  and `x_yield_measured` = F_ferm/F_sugar_fed (up to 1.65 y).
* **A carbon-basis yield is derivable but not better.** Glucose, xylose and
  acetic acid all hold 40.0 wt% carbon, so y_C = F_ferm_product /
  (F_sugar_fed + x_F_acetate_fed) is the yield per g of fermentable carbon
  (glucose-eq.). But y_C/y = 1.00-1.13 on the development rows, set by the
  acetate share: 1/y_C = 1/y − (a/2)/t, with a = x_F_acetate_fed/V_broth
  varying between points (1.8-10.4 g/L; 10 g/L AcOH-eq. by default, changed
  by the module's titer loader for low-titer, high-yield points). In y_C the
  TY design becomes irregular and the MPSP leaves the unifying form (the
  substitution adds 1/t² and acetate-loading terms); sampling a carbon-basis
  yield directly would need a different yield specification in the module.
  Side note: the acetate reaction is carbon-balanced but creates mass
  (C2H4O2 → 1/3 C6H6O3 + 3 H2O: 60 g → 96 g), a module quirk akin to D14.
* **How a user specifies TAL.** A user quotes a yield per sugar consumed (or,
  like the module, as a fraction of theoretical on all substrates) and
  treats an acetate co-feed as a separate process choice, which the model
  has no input for. Under either basis the training rows carry the
  acetate's sugar substitution: sugar demand per kg TAL is 1/y − 1.5·a/t, a
  negative titer term in the feedstock channel whose coefficient changes with
  a between points: the fitted feedstock-channel c′ is negative in every
  development TAL sample (188 samples after the regeneration; median −3.0,
  range −7.3 to −0.7 USD/kg·g/L) and d′ is positive (median 0.33), against c′ ≈ d′ ≈ 0 for the other products
  (medians 0.000 / 0.002).
* **Decision.** y semantics unchanged: the recorded y is the module's own and
  consistent with ADAPTERS.md, and neither the carbon basis nor a sugar-only
  basis is clearly more correct (both change the functional form or need
  module changes). The data carry what a mitigation needs
  (`x_F_acetate_fed`, `F_sugar_fed`, `V_broth` at every point). Options were
  to mask the TAL / TAL_SA feedstock-channel targets (`ch__feedstock__*`,
  `int__fs__*`), to normalise that channel by the sugar-equivalent feed
  F_sugar_fed + 1.5·x_F_acetate_fed, or to add an acetate co-feed input.
  Round 5 (§6) adopted the mask (`experiments.DEFECTIVE_TARGETS`): the
  product-blind feedstock head cannot represent the substitution term
  anyway, and keeping it out of the shared head is what the other products
  of corn and corn stover need; the normalisation and the co-feed input
  stay open.

Minor items, also not corrected: oxalic recoveries slightly above the
stoichiometric ratio in 26 of 211 samples (max 1.003; lime loading > 1 gives
negative free acid in the module; ≤ 0.3 % on the coefficients; the
`max_g_over_stoich = 1.05` filter does not catch it); succinic sugar
concentrations above the module's own 600 g/L limit at 12-19 % of the points
of the non-corn-stover configurations (the module's check is never called;
refitting without them changes MPSP on the feasible points by a median of
0.1-0.9 %); the same for the TAL glucose configurations, whose module skips
its sugar-concentration check for glucose
(`system_TAL_solubility_exploit_ethanol_glucose.py`: `check_sugar_concentration`
only `if not feedstock_ID == 'Glucose'`): sugar per broth volume above
600 g/L at 201 / 656 TAL_glucose and 189 / 650 TAL_SA_glucose points (max
about 1220 g/L; corn and sugarcane ≤ 660 g/L), and at the low-yield /
high-titer corner (y 0.07-0.09, t 83-100 g/L) the recovery collapses 15-30 %
below the sample median and MPSP is 12-24 % above the 4-coefficient fit
(e.g. TAL|TAL_SA_glucose|s1|6 at y 0.088, t 100: g 0.470 vs 0.679 at t 71,
MPSP 14.1 vs fit 11.1; samples TAL_SA_glucose s1|6, 7, 16, 21, 28 and
TAL_glucose s1|14, 25, 26; dropping the corner point changes the fitted MPSP
on the other points of these samples by 1-3 %, max 7 %; module behaviour,
as in the publication, kept); 13 retained samples (counted before the H1
regeneration of TAL_SA_glucose, whose new samples add the corner anomalies
above) with an isolated interior MPSP spike of more than 10 % along titer
(TAL_glucose 7, TAL_corn 2, lactic 3, HP_neutral_corn 1; likely convergence
or sizing artefacts, kept); fixed module reagent prices not recorded as
inputs: succinic_neutral's H2SO4 (0.0948 USD/kg at every point; the
neutralization channel is 12-14 % of its MPSP) and TAL's NaOH (5-11 % of
MPSP) leave `p_acid_price` / `p_base_price` NaN (imputed, missing indicator
set) while HP, lactic and oxalic sample them -- to be recorded by the
adapters at the next data rebuild; the sampled gypsum price puts the gypsum
value (−1.1 to +3.7 % of MPSP) in the coproducts channel although
`downstream.TABLE` lists only electricity as a coproduct of every
gypsum-producing configuration (the `ds_cop__gypsum` term is 0 for every
training row; the convention is the same for all of them); the lactic SSCF random
samples rejected for too few feasible points are the low-carbohydrate ones
(retained median carbohydrate content 0.54, rejected 0.25-0.59), so the
retained data are censored toward carbohydrate-rich feedstock.

### Product-specific cost structure (not defects)

Two configurations break priors that hold for every other one; the data are
the modules' own physics, kept as simulated, but they matter for the model:

* **Isobutanol (`IBO_EtOH_corn`): the MPSP rises with the yield.** The
  adapter (`adapters/isobutanol.py`) sends the sugar not converted to
  isobutanol to ethanol, sold at a fixed price (`x_ethanol_price` 0.79
  USD/kg), so its yield axis is a product split, unlike the other families,
  whose yield deficit goes to unsold cell mass or byproducts. At low yield
  the ethanol credit exceeds the feedstock cost (baseline sample, y = 0.05,
  t = 5 g/L: coproducts −8.58 USD/kg, feedstock +5.57 USD/kg). Consequences:
  b′ < 0 in 102 of the 109 samples (median −0.021) and d′ < 0 in all 109
  (median −0.169), against b′ = "the sugar needed per kg product" for the
  others (README section 1); the simulated MPSP rises from y = 0.05 to 0.36
  on the t = 5 g/L line in 99 % of the samples (median +33 %), and the
  relative impact of yield vs. titer (README section 4.4) is negative at
  the isobutanol baseline (a yield gain RAISES the MPSP).
* **Isobutanol: broth alcohol titers without an inhibition limit.** The
  adapter replaces the module's kinetic model by a stoichiometric one
  without product inhibition: 172 of the 2712 retained isobutanol points
  (106 samples) have more than 120 g/L ethanol in the broth (max 148 g/L),
  283 more than 150 g/L total alcohol (max 168 g/L; glucose loading up to
  578 g/L). Refitting each sample without its points above 120 g/L ethanol
  (OLS of g·MPSP on the unifying basis, which reproduces the stored
  coefficients) changes the mean d′ from −0.168 to −0.118 (−30 %); the
  fitted MPSP on the kept points changes by 1.1 % on average, the
  per-sample maximum by a median 5.1 % (up to 13 %). Not corrected (a
  titer cap in the adapter's `max_titer`, or flagging those points, is an
  option for the next data rebuild).
* **Oxalic acid (`oxalic_neutral_sugarcane`): the MPSP rises with the titer
  at high titer.** Reaching a higher titer means evaporating the cane
  juice, which cuts the exported electricity (baseline sample, y = 0.24:
  electricity credit −0.600 → −0.516 USD/kg and MPSP 1.604 → 1.678 USD/kg
  from t = 52.5 to 100 g/L); every oxalic sample has a yield line on which
  the simulated MPSP rises with the titer (between neighbouring titers by
  up to +9.7 %; from the line's minimum to its highest titer by a median
  +8.9 % of the samples' largest rise, up to +15 %). 26 oxalic
  samples also have a recovery slightly above the stoichiometric ratio 1
  (g up to 1.003; the predictions are bounded by g_max).

Effect on the physics residuals: `train.check_constraints` enables a
residual only if the training targets satisfy it in ≥ 99 % of the samples.
dMPSP/dy ≤ 0 holds in 93.25 % (the 108 violating samples are all
`IBO_EtOH_corn`) and dMPSP/dt ≤ 0 in 98.1 % (the 30 violating samples are all
oxalic), so both residuals are disabled for EVERY product
(`models/release/meta.json`, `constraint_report`). The REPORT.md
constraint-violation tables (sections 5 and 6) therefore show the data's own
rates (`data (fit)`) next to the methods': a `viol_mono_y` near the data's is
the isobutanol share of the held-out samples, not a model error. Enabling the
monotonicity residuals per product (off only for isobutanol / oxalic) would
change training and is left for a future model round (it needs a re-run and
a new evaluation).

The natural-gas intensity has the same kind of issue. The simulated
natural-gas use `q__ng` is never negative (0 of 35 583 points), but it is a
hinge, max(0, ·): zero wherever the biomass fuel covers the steam demand
(44 % of all points; 6-83 % of the HP points per configuration, 80 % of
ethanol_cornstover's, 92 % of ethanol_sugarcane's, 89 % of oxalic's) and
positive only at low titer or high yield. The unifying basis cannot
represent a hinge, so the 4-term fits `int__ng__{a,b,c,d}` (the targets of
the natural-gas intensity head) go below zero at some of the sample's own
points in 923 of the 1600 samples (`nonneg_ng` satisfied in 42 %): the
non-negativity residual of the natural-gas intensity is disabled for that
reason, and a predicted natural-gas intensity (`int__ng__*`) can be negative
in parts of the domain, as the fitted targets are. The channel identities
are unaffected (`price_linear` holds exactly), so the totals and the price
linearity are; clip the natural-gas use at 0 when reporting it as a
quantity.

### Measured vs. nominal yield and titer

ADAPTERS.md checklist item 2 asks that y and t follow the shared
definitions. REPORT.md section 1 reports the per-point ratios
F_ferm_product / F_sugar_fed / y and F_ferm_product / V_broth / t per
configuration (`results/tables/yt_definition_check.csv`, report only). HP,
succinic, oxalic and isobutanol stay within ±0.6 % (ethanol_cornstover within
0.9 %). The departures above 1 %:

| configurations | yield ratio | titer ratio | cause |
|---|---|---|---|
| TAL, TAL_SA (8) | 0.998-1.645 | 1.000 | D8 (the spiked acetate counts as substrate in y, not in F_sugar_fed) |
| lactic SHF (2) | 0.999-1.056 | 0.997-1.002 | residual seed-train sugars recycled with the inoculum (adapter docstring) |
| lactic SSCF (2) | 0.999-1.047 | 0.976-1.001 | the same for the yield; the titer ratio below 1 is not explained |
| ethanol (3) | 1.000 | 1.001-1.026 | ethanol leaving with the CO2 vent counts in F_ferm_product but not in the module's broth titer (adapter docstring) |

The coefficients are fitted on the nominal y and t (the modules' own
definitions), the inputs of the unifying equation; the ratios matter only
when comparing yields or titers across modules.

### Comparison with the publication

REPORT.md section 8 (i) compares the baseline coefficients with the
publication's (`results/data/<config>_coefficients.npy`): b′, c′, d′ and g
agree within a few per cent for most configurations, a′ is higher (median
ratio 1.215). Two causes are known. Succinic acid at neutral pH differs by
design (D15: the regenerated data simulate the salt formation and
acidulation; a′ ratio 1.24-1.34 against 1.02-1.15 at low pH). HP's a′ is
1.08-1.63 x the publication's, mostly fermentation capital (a median 48 %
of the HP baseline a′, `ch__capex_fermentation__a`), which scales with 1/r:
a review comparison with the publication's productivity coefficient sets
(`results/data/<cfg>_0.2bp_*` and `_5.0bp_*`, at the median y and t of each
baseline sample, release model) found the MPSP increase at 0.2 x the
baseline productivity 1.24-1.64 x the publication's for the 16 HP
configurations (HP_corn +1.148 vs. +0.825 USD/kg) and the decrease at 5 x
1.22-1.63 x, while TAL and succinic agree within 0.89-1.11; the data show
the same (HP_corn random sample at r = 0.057: a′ 0.82 vs. the publication's
0.2bp a′ 0.596). The patches change no process value
(`data_generation/bip_patches/README.md`), so the difference comes from the
module or BioSTEAM version; a candidate, not verified, is the fermenter
design / costing of `BatchBioreactor` in BioSTEAM 2.54, whose vessel-size
argument the HP patch adapts to (`V` → `V_max`). The specific change has
not been identified. HP's productivity response therefore carries this
uncertainty.

### Data set (after the regeneration of 2026-10-17)

Processed tables (`experiments.stage_data`, i.e. `dataset.build` with the
default filters; data hash c345efda41c4fdca): 1600 samples, 35 583
feasible points, 41 configurations of 7 biorefineries (HP 464, TAL 232,
ethanol 327, isobutanol 109, lactic 113, oxalic 211, succinic 144);
development set 1353 samples / 37 configurations. Every configuration has
completed the samples planned by `run_all.sh` except
`lactic_lowpH_SHF_cornstover` (27 of 36); all 23 quality-filter rejections
are lactic (11 too few points, 9 too few yields, 2 low R², 1 incomplete).
The v2 models were trained on the tables before the regeneration: 1550
samples (TAL 185: 47 TAL_SA samples rejected for too few points, the H1
empty attempts; succinic 141 with the D15 neutral data).

Regenerated configurations (same seed, sample ids and campaign settings as
`run_all.sh`): `succinic_neutral_{corn, cornstover, glucose, sugarcane}`
18 / 18 samples each (adapter v3, 288 ok points each, no failed point),
`TAL_SA_glucose` and `TAL_SA_sugarcane` 29 / 29 (28 and 19 samples re-run).
Fit MAPE, median per configuration on the development rows: succinic_neutral
0.86-1.16 % (max 1.53 %; low pH 1.01-1.38 %), TAL_SA_glucose 1.12 % (max
3.38 %), TAL_SA_sugarcane 0.90 % (other TAL_SA 1.01-1.15 %). Patches: the
re-run TAL_SA attempts (adapter v2, harness v3) and the succinic_neutral
v3 rows receive S1 only; the kept original TAL_SA attempts (adapter v1,
harness v2: 1 + 10 samples) keep D3 / M1 / M3 / M12 / S1. The two TAL_SA
files now hold an adapter-v1 and an adapter-v2 config record each
(`stats['config_conflicts']` = 2, expected: each sample uses the record that
precedes it in its file).

## 2. Physics layer 1 — exact cost-channel decomposition (system-level data)

At each TY point `i` of a sample, the product revenue required at the MPSP is
`R_i = MPSP_i · F_final,i · h` (USD/yr; `h` operating hours). BioSTEAM's
discounted-cash-flow MPSP is (to within tax non-linearities) linear in annual
costs and capital, so

    R_i = Σ_c  E_c,i  +  FOC_i  +  φ_K·I_i  +  ε_i,     FOC_i = F + (FOC_i − F)

where `E_c,i` are annual operating cost channels measured directly from the
simulation, `FOC_i` the TEA fixed operating cost (measured), `I_i` the fixed
capital investment (TEA `FCI`; installed equipment cost is used only if FCI
is unavailable — FCI/installed varies by 1.69–1.82 across TY points), `F` the
labor-type fixed costs, `FOC_i − F` the capital-related O&M (maintenance,
insurance, property tax, administration), and `φ_K` the capital charge per
USD of FCI (annualized capital recovery at the IRR incl. taxes and working
capital) from the no-intercept least squares of `R − ΣE − FOC` on `I`
(financial parameters are constant within a sample). `φ = φ_K + d(FOC −
F)/dI` is reported. `ε_i` is the closure residual (kept as its own channel
so that channels sum exactly to MPSP); it absorbs DCF weights ≠ 1 (e.g.
startup fractions) instead of letting a free intercept turn them into a
negative or productivity-dependent `F`.

The fixed channel `F` is measured (harness v2 records the TEA's
`labor_cost` and `labor_burden`), following how the Bioindustrial-Park TEA
classes build FOC:

* `CellulosicEthanolTEA` (HP, TAL, lactic, succinic, oxalic, cellulosic
  ethanol): `FOC = FCI·property_insurance + ISBL_DPI·maintenance +
  labor_cost·(1 + labor_burden)` → `F = labor_cost·(1 + labor_burden)`
  (`diag_fixed_source = 'labor'`);
* `ConventionalEthanolTEA` (corn ethanol, isobutanol): `FOC =
  FCI·(property_tax + property_insurance + maintenance + administration) +
  labor_cost·(1 + fringe_benefits + supplies)`; it has no `labor_burden`, so
  `F = m·labor_cost` with `m` from the least squares `FOC = s·I + m·labor`
  (exact: FOC is affine in FCI; real data give m = 1.6) (`'labor_foc'`);
* legacy (v1) data without labor: `F` = intercept of FOC regressed on `I`
  (`'foc_intercept'`). This was biased for the cellulosic TEA, where
  maintenance is charged on ISBL installed cost (not ∝ FCI): on TAL v2
  data the intercept was −1.35 × labor (all samples had a negative `F`), on
  oxalic 1.08 × labor instead of 1.9 ×.

Without FOC the legacy free fit `R − ΣE = φ·I + F` is used. Samples with a
negative `F` or with measured labor above FOC (`diag_fixed_exceeds_foc`) are
flagged (`channels_ok = False`).

Channels (contribution to MPSP = channel USD/yr ÷ (F_final·h)):

| channel | content | price-explicit? |
|---|---|---|
| `feedstock` | wet feedstock purchase | yes: `p_fs × feedstock intensity` |
| `electricity` | net purchased (+) / exported (−) electricity | yes: `p_el × net electricity intensity` |
| `natural_gas` | purchased (priced) natural gas, both cost pools | yes: `p_ng × priced NG intensity` |
| `neutralization` | base for pH control + acid for acidulation | no |
| `materials_other` | all other purchased materials (enzymes, nutrients, catalysts, solvents, boiler chemicals …) | no |
| `utilities_other` | non-electric utility agents (net, after boiler offsets) | no |
| `capex_upstream` | capital charge (φ_K·FCI + FOC − F) × weight share of feedstock + conversion areas | no |
| `capex_fermentation` | capital charge × weight share of fermentation | no |
| `capex_downstream` | capital charge × weight share of separation + upgrading | no |
| `capex_wastewater` | capital charge × weight share of wastewater | no |
| `capex_utilities` | capital charge × weight share of CHP + facilities (+ uncovered installed cost) | no |
| `fixed` | F: labor-type fixed costs (recorded labor × (1 + burden)) | no |
| `coproducts` | − coproduct sales | no |
| `closure` | ε (should be small) | no |

**Capital weights.** The capital charge is split over the capex channels by
each area's contribution to FCI. In `CellulosicEthanolTEA`, DPI =
installed + `f`·ISBL installed (`f` = warehouse + site development +
additional piping, ISBL only) and every indirect cost (proratable, field,
construction, contingency, other) is proportional to DPI, so `FCI_a ∝
I_isbl,a·(1 + f) + I_osbl,a`. Harness v2 records each area's OSBL installed
cost (`<area>__installed_cost_osbl`, from `tea.OSBL_units`); `f` is
estimated per sample from `FCI_i = c·(I_i + f·ISBL_i)` (real data: 0.175).
Without OSBL data (v1, or a TEA without `OSBL_units` such as
`ConventionalEthanolTEA`, whose FCI/installed = 1) installed costs are used
(`diag_capex_weights = 'tdc' | 'installed'`). The ISBL-proportional part of
`FOC − F` (maintenance on ISBL DPI, from the least squares `FOC − F =
s_K·FCI + s_S·ISBL`) is allocated by ISBL installed cost
(`diag_foc_alloc = 'isbl'`, `diag_foc_isbl_rate` = 0.03 on real data); the
rest goes by the FCI weights. Installed cost not covered by the areas goes to
`capex_utilities`.

**Natural gas (harness v2).** BioSTEAM charges an inlet stream registered
as a 'natural gas' utility at `bst.stream_prices[name]` if the stream has
no price (cost pool *utility* → TEA `utility_cost`). Otherwise the stream's
own price applies (pool *material* → TEA `material_cost`), as for any other
mostly-methane feed, which may be unpriced (e.g. TAL's dryer gas). The
harness records `ng_flow_utility, ng_cost_utility, ng_flow_material,
ng_cost_material` (kg/h, USD/h, at the charged prices). These replace the
adapters' `natural_gas_flow/price`, which proved unreliable: TAL priced its
unpriced dryer gas, and corn ethanol reported price 0 for utility gas
charged at 0.218 USD/kg.

* channel `natural_gas` = `ng_cost_utility + ng_cost_material`. Each part is
  removed from its own remainder (`utilities_other`, `materials_other`), so
  nothing is counted twice;
* intensity = flow of the **priced** gas only. Per pool, priced flow =
  `cost / p_pool`, where `p_pool` is the pool's highest unit price
  (cost/flow) over the sample's points. A pool holding only priced (or only
  unpriced) streams gives exactly its flow if its cost > 0, else 0.
  *Approximation:* in a pool that mixes priced and unpriced streams, the
  unit price varies with the mix, and its maximum is the priced streams'
  price only if some point is nearly free of unpriced gas. Otherwise it is a
  lower bound, and the intensity an upper bound (`diag_ng_mixed_pool`). The
  plain pool flow would inflate succinic's intensity 4–40× (a 0.2527 USD/kg
  CH4 feed beside ≤ 1e5 kg/h of unpriced CH4); with `p_pool` the price
  matches the priced stream;
* `natural_gas_price` = total NG cost / total priced flow of the sample (the
  effective price). `int__ng` uses the cost-equivalent priced flow `cost /
  price`, so `price × intensity` reproduces the NG cost exactly at every
  point. This equals the priced flow when a single price applies; otherwise
  see `diag_ng_q_equiv_rel_dev`. A sample with no priced NG has channel = 0
  and intensity = 0, and the price falls back to the adapter's value
  (`diag_ng_price_source = 'adapter'`).

Legacy (v1) data, with one adapter-recorded price, use a pool decided once
per configuration (CHP-area material-cost evidence, else the sign of the
remainders). Samples whose chosen pool leaves a negative remainder
(mispriced or split NG) get `channels_ok = False`. The vote is skipped for v2
samples.

Each channel's contribution `m_c(y, t)` is fitted with the same OLS design
matrix `[1, 1/y, 1/t, 1/(y·t)]` as the total; by linearity of OLS the channel
coefficients sum **exactly** to the total coefficients. For the price-explicit
channels the *intensity* (e.g. kg feedstock per kg final product, kWh per kg)
is fitted instead, and the channel coefficient = price × intensity
coefficient. All channel coefficients are converted to the primed (×g) basis.

## 3. Physics layer 2 — structured model outputs (hard constraints)

A shared encoder maps standardized features → latent `h`. Heads output:

* intensity coefficients `q_fs, q_el, q_ng` (4 each) for the price-explicit
  channels → coefficients `p_j · q_j` (exact linear price dependence: costs
  are linear in prices at fixed process design);
  - feedstock intensity is structured: `b'_fs = exp(z) / w_sugar` where
    `w_sugar` is the feedstock's glucose-eq. fermentable sugar mass fraction
    and `exp(z) = 1/η ≥ 1` (η = fraction of potential sugar reaching
    fermentation), other feedstock-intensity coefficients are free but small;
* coefficients of every other channel (4 each);
* recovery `g = g_max · sigmoid(z_g)`, `g_max` = stoichiometric mass ratio of
  final to fermentation product (hard upper bound).

Totals `a', b', c', d'` = Σ channels. The **equation layer** turns totals into
MPSP at any (y, t). Raw head outputs are scaled by per-target training
standard deviations for conditioning.

**Water-removal gate** (`ModelConfig.water_gate`, off by default; needs
`FeatureSpec(downstream=True)`): the c′ and d′ coefficients of the channels
that scale with the broth water processed downstream (`capex_downstream`,
`capex_utilities`, `capex_wastewater`, `utilities_other` and the natural-gas
intensity) are multiplied by the raw structural input `struct_water_removal`
(0 when the broth itself is the product, route `none`). This is a hard
physical constraint: without evaporation, distillation or treatment of the
broth water these titer-proportional costs vanish (oxalic acid, the only
broth product: NG c′ 0.24, downstream/utility/wastewater capex c′ ≤ 0.10,
versus 5–18 and 1–2 for products whose water is removed). The structural
input is never an encoder input, so a value unseen in training (0 when the
broth product is held out) cannot shift the latent features.

**Input clipping** (`ModelConfig.clip_inputs` ≥ 0 = margin in training
standard deviations; default −1 = off; an HPO search dimension): the standardized continuous encoder
features are clipped to their training range (`x_lo`/`x_hi` in
`target_scales.json`), so a query whose descriptors lie outside the training
products' range (e.g. oxalic acid's degree of reduction or pKa when no other
product is like it) is predicted as at the nearest edge of the training
range instead of by unbounded linear extrapolation of the encoder.
Productivity (`log r`, `1/r`) and the price-explicit prices are never
clipped, so the productivity physics and the exact price linearity are
unchanged; the applicability-domain report still flags such queries.
Checkpoints saved before this option existed load with clipping off. In
the development CV it helped leave-configurations-out (oxalic) but hurt
leave-one-product-out for HP and ethanol (§6), so it is not the default.
`ModelConfig.clip_exempt` (round 5) lists columns that are never clipped:
the pipeline exempts the exact final-product stoichiometric ratio and its
MW-ratio twin, whose clipping caused the HP failure (§6, round 5, 5d-5e).

**Channel-wise input masking** (`ModelConfig.channel_inputs`; library default
`None` = one shared encoder; `experiments.py` default `'physics'` since
2026-10-17, CLI `--no-masking`). With a shared encoder every channel can
depend on every product-identity proxy, including channels whose process
area never handles the product (feedstock purchase, feedstock handling and
pretreatment, labor). With masking, every feature column belongs to a
process-area group (`features.feature_group`, `FeatureSpec.column_groups`):

| group | columns |
|---|---|
| feedstock | `fs_w_*`, `feedstock_type`, `phys_fs_nonsugar_solids` |
| pretreatment | `p_enzyme_loading`, `p_saccharification_yield` |
| scale / financial | `phys_log10_capacity`, `p_operating_days` / `p_IRR`, `p_income_tax` |
| fermentation | `T_fermentation`, `aerobic`, `fermentation_mode`, `p_cell_mass_yield`, `phys_inflection_over_theoretical` |
| neutralization | neutralization flag / agent / degree, `n_acid_groups`, `pKa1`, `phys_{base,acid,gypsum}_per_product`, `ds_sep__acidulation`, `ds_sep__gypsum_filtration` |
| reagent_prices | `p_{base,acid}_price` (split from neutralization: they set what the reagents cost, not how much product is recovered) |
| product | `fp_*`, `phys_fp_*`, theoretical / inflection yield, `phys_ty_over_*`, solubility |
| final_product | `pp_*`, `phys_pp_*`, stoichiometric ratio, phase, broth flag, `phys_log10_w_final` |
| upgrading | `upgrading`, `ds_upg__*`, `phys_n_upgrading_ops` |
| separation | `ds_sep__*` (the operations, which identify the product), `phys_n_separation_ops` |
| water_route | `water_removal_route`, `phys_water_removal_{a,c}`, `phys_evap_water_c`, `phys_distill_c` |
| coproducts / utilities | `ds_cop__*` / boiler and turbogenerator efficiencies |

(missing indicators and one-hot levels go with their field; unknown columns
are 'other', used by no preset output, with a warning; `struct_*` never). Each
distinct set of groups gets its own sub-network (`masked_width` 64,
`masked_depth` 2; the pipeline uses width/2 and depth − 1 of the HPO
candidate) fed only those columns, plus log r and 1/r for the pseudo-group
'productivity'; every output (intensity, channel, g) is a head on the
sub-network of its groups (`CHANNEL_INPUT_PRESETS['physics']`; `'physics_v2'` differs only in the
recovery, which does not see the reagent prices -- §6, round 4):

| output | groups |
|---|---|
| feedstock intensity | feedstock, pretreatment |
| capex_upstream | feedstock, pretreatment, scale, financial |
| fixed, closure | feedstock, scale, financial |
| capex_fermentation | fermentation, productivity, neutralization, reagent_prices, feedstock, scale, financial |
| neutralization | neutralization, reagent_prices |
| capex_downstream | product, final_product, upgrading, separation, water_route, neutralization, reagent_prices, scale, financial |
| utilities_other | product, final_product, upgrading, separation, water_route, fermentation, utilities, scale |
| natural-gas intensity | product, final_product, upgrading, separation, water_route, feedstock, utilities, scale |
| materials_other | feedstock, pretreatment, fermentation, product, final_product, upgrading, separation, water_route, scale |
| capex_wastewater | feedstock, fermentation, water_route, coproducts, scale, financial |
| capex_utilities | feedstock, water_route, upgrading, coproducts, utilities, scale, financial |
| electricity intensity | feedstock, fermentation, productivity, water_route, upgrading, coproducts, utilities, scale |
| coproducts | feedstock, fermentation, coproducts, scale |
| g | product, final_product, upgrading, separation, water_route, neutralization, reagent_prices (`'physics_v2'`, the pipeline default since round 4: without reagent_prices) |

An excluded input has exactly zero influence on an output (∂output/∂input = 0,
tested); the water gate, the structured feedstock intensity, the g bound and
the sum identity apply unchanged. A price channel resolved as an ordinary
channel (`price_explicit=False`) also sees its price. The hybrid
(`residual_width > 0`, `residual_depth`) adds the standardized outputs of a
shared residual network on all inputs (zero-initialized head; its weights
get the AdamW weight decay `residual_weight_decay`).

Data checks (development set, configuration medians of each channel's
numerator at y = 0.2, t = 50 g/L; R² of the across-configuration variance
explained by the feedstock alone / by the product family alone): feedstock
intensity 0.99 / 0.27, capex_upstream 0.90 / 0.16, capex_utilities 0.73 /
0.13, versus capex_downstream 0.03 / 0.91, neutralization 0.03 / 1.00,
natural gas 0.09 / 0.90, capex_fermentation 0.09 / 0.79, materials_other 0.25 /
0.82; electricity 0.53 / 0.43 and capex_wastewater 0.51 / 0.56 depend on both.
Over the productivity samples of a configuration (only r varies) only
capex_fermentation changes materially (range up to 47 % of the total
numerator) and electricity ≤ 1.8 %, every other channel ≤ 1.1 %, hence
'productivity' goes to those two only. Product-specific exceptions that the
masked channels cannot represent (modelling conventions or known defects):
TAL's feedstock intensity has c′ = −17 … −90 (substrate basis, D8) and HP's a
d′ of 0.8–5 (all others ≈ pure b′ = 1/(η w_sugar)); the cellulosic and
sugarcane ethanol modules book less / more capital in the upstream areas
(0.15 / 0.17 vs 0.27 / 0.06 USD/kg for the other products); labor is set per
module (TAL ≈ 1/3 of HP at the same capacity, succinic corn / cornstover
1.8·10⁷ USD/yr vs 1.1·10⁶ for succinic sugarcane); the neutralization channel
of HP (also non-neutral) and lactic has a broth-volume term c′ ≈ 2.5–3 that
the stoichiometric features do not explain.

Evidence (§6, round 3; 5 seeds): leave-one-product-out improves for HP
(24 → 15 %) and TAL (61 → 42 %), paired over seeds beyond the seed noise
(t = −2.6, −10), the other families are unchanged within noise (ethanol 74 →
83 %, IBO 86 → 96 %, oxalic 74 → 79 %, succinic 33 → 35 %, lactic 56 → 54 %);
leave-configurations-out is not worse (configuration mean 29.5 → 27.9 %,
geo-mean 21.2 → 18.4 %, better on 4 of 5 seeds). The channels whose
physics the masks encode transfer better to held-out products (out-of-fold
error relative to the total numerator, LOPO): feedstock 11.3 → 7.1 %,
capex_upstream 3.3 → 1.7 %, coproducts 29.6 → 13.0 %, materials_other 16.5 →
13.4 %, capex_fermentation 5.0 → 4.3 %, fixed 4.6 → 3.9 %; neutralization
gets worse (6.7 → 9.2 %: its broth-volume term above), electricity, natural
gas and the other capital channels are unchanged. The hybrid (residual width
32, weight decay 1 or 10) is not better than pure masking (ethanol LOPO 104 /
110 %), so it stays opt-in. Products without a close relative (ethanol, IBO,
oxalic, lactic) are not helped by masking; their errors come from the
recovery g and from electricity / natural gas / downstream channels that
legitimately depend on the product. Checkpoints saved before the option load
with one shared encoder.

**Hard productivity structure** (`ModelConfig.r_structure`; library default
off; on in the v2 pipeline, `experiments.R_STRUCTURE`, CLI `--no-r-structure`). Each listed output (preset `'physics'`:
`capex_fermentation` and the electricity *intensity*; a price-explicit
channel is structured before the price, so the price linearity is exact)
has coefficients

    θ_k(x, r) = θ0_k(x) + θ1_k(x) · (r_ref / r)^β,   k = a′, b′, c′, d′,

with `θ1_k = s_k · softplus(·) ≥ 0`, `β` fixed (`r_beta` = 1) or learned in
`r_beta_bounds` (one per output), `r_ref` = exp(training mean of log r)
(`FeatureSpec.r_stats`; only a reparametrization: the head of θ0 gives the
coefficient at `r_ref`), and θ0, θ1 computed without the productivity inputs
('productivity' removed from those outputs' groups in the masked model; the
shared encoder and the hybrid's residual network get no productivity input at
all). With `g` independent of r (the preset never gives 'productivity' to
the recovery) every structured term of MPSP = Σ/g falls with r for all y, t
> 0, so ∂MPSP/∂r ≤ 0 holds by construction, also far outside the training
range (the soft `mono_r` residual is then identically 0 for the preset; it
still applies if a custom `channel_inputs` gives another output the
productivity). θ1 is initialized at the pooled within-configuration slope of
the training coefficients on `r_ref/r` (target scales `*_r_slope`,
`train.r_slopes`), at least `r_slope_floor` = 1 % of the target's standard
deviation, so a coefficient without productivity dependence starts (and,
without data pressure, stays) with a negligible term: a θ1 of the order of
the cross-configuration standard deviation is amplified up to (r_ref/r_min)
≈ 20× at the low-productivity edge (a first prototype initialized at 0.3 σ
had a much larger out-of-fold electricity error in a single-member run:
38 vs 23 % of the numerator). The gate multiplies the whole
structured coefficient (≥ 0, monotonicity kept). The outputs expose
`r_theta0`, `r_theta1` (n, S, 4; before the gate). Checkpoints saved before
the option load with it off.

Data (development set only; `sample_kind` productivity/baseline = only r
differs within a configuration; channel numerators on a 5 × 5 TY grid of the
configuration, relative to the total numerator):

* `capex_fermentation` is the only channel with material r dependence:
  range over r/r_base ≈ 0.25–4 of median 13 % (max 41 %, oxalic) of the
  total numerator, monotone decreasing in 35 of 36 configurations (one
  step up of 0.04 % of the numerator, lactic). Fitting
  θ0 + θ1 (r_ref/r)^β per configuration and coefficient gives β = 1.0
  (IQR 0.95–1.0; a′ 1.0, b′ 0.8, c′ 0.9) with R² ≈ 1.000; a joint fit with
  one β for the four coefficients: R² 0.9999 at β = 1.0 vs 0.997 at 0.8 /
  1.2 (RMS error 0.06 % of the numerator vs 0.2–0.3 %). The dependence is
  almost purely a′ (median a′-term 12 % of the numerator; b′, c′ ≤ 4 %), as
  V ∝ product rate / r makes the capital per kg ∝ 1/r independently of y
  and t; β = 1 rather than a capital exponent < 1 because the fermenters
  are trains of vessels of bounded size. θ1 ≥ 0 for every coefficient (NNLS)
  costs R² 0.0004 (d′ has a slightly negative slope in 28 / 36
  configurations, |θ1_d| ≤ 0.4 % of the numerator).
* Electricity: range ≤ 1.4 % (median 0.4 %), non-monotone between
  consecutive r in 27 / 36 configurations (per-sample fit noise); its a′
  follows ~1/r (β ≈ 1.2, R² 0.90, θ1 > 0 in 30 / 36), consistent with
  agitation / aeration power ∝ V. All other channels ≤ 1.2 % without a
  consistent sign. The recovery g varies by median 0.2 % (max 4.5 %) and is
  uncorrelated with r.
* LHS random samples (r sampled jointly with 17–48 other parameters; ridge on
  the varying parameters + θ1 (r_ref/r)^β, β by leave-one-out): fermentation
  capital r effect median 23 % (max 50 %) of the numerator, partial R² 0.99,
  β median 0.9 (IQR 0.8–1.03), θ1 > 0 in all 28 configurations, a free cubic
  in log r monotone decreasing in 25 / 28 (in the other 3 the rising part is
  ≤ 2 % of its range); with θ1 depending on the other parameters
  (large configurations) β = 1.0 (oxalic), 0.8–0.9 (ethanol, IBO).
  Electricity: partial R² ≤ 0 in most configurations (r explains nothing
  beyond the other parameters), median effect 0.5 %. The closure channel
  depends weakly on r (0.15 % of the numerator; DCF weights on the
  r-dependent capital) and is left unstructured.

Evidence (development CV, round 4, 2026-10-17; PI-NN with the selected
hyperparameters, water gate and 'physics' masking, ± r-structure with β = 1;
paired by seed; MPSP MAPE of `compute_metrics`, pooled / configuration
balanced; Δ = paired difference, t over seeds):

| scheme (members, seeds) | v1 | + r-structure | Δ (t) |
|---|---|---|---|
| productivity, sample-level inner hold-outs, no calibration hold-out (2, 5) | 5.16 / 5.27 | 4.29 / 4.33 | −0.87 (−4.5) / −0.94 (−7.4) |
| … low-r side, configuration balanced | 6.34 | 4.39 | −1.96 (−12.5) |
| … high-r side, configuration balanced | 4.22 | 4.20 | −0.01 (−0.1) |
| productivity, v1 protocol (2, 5) | 15.0 / 11.1 | 15.2 / 9.7 | +0.3 (+0.3) / −1.4 (−2.4) |
| leave-configurations-out, 5 folds (2, 5) | 32.9 / 17.7 | 33.4 / 17.3 | +0.5 (+0.5) / −0.3 (−0.6) |
| random 5-fold (1, 5) | 4.34 / 4.27 | 4.27 / 4.18 | −0.07 (−1.1) / −0.09 (−2.1) |

Ridge / HGB on the productivity scheme: 8.9 / 7.5 % (sample-level hold-outs),
20.9 / 24.3 % (v1 protocol). The finite-difference ∂MPSP/∂r > 0 rate of
`compute_metrics` is 0 with the structure (v1: 6.4 % of the held-out
samples in the productivity scheme, 3.8 % LOCO, 0.8 % random). Out-of-fold
error of the fermentation capital (RMS over the points / RMS total
numerator, configuration balanced): 1.57 → 1.14 % (productivity, t −13),
2.27 → 1.96 % (LOCO, t −2.2), 0.90 → 0.95 % (random, t +1.9: the 1/r form
with θ1 ≥ 0 costs a little in interpolation, cf. β ≈ 0.9 and d′ slopes < 0
above). Electricity is unchanged where the configurations are in training
(productivity 7.2 → 7.4 %, t +0.3; random 7.6 → 7.5 %) and not
significantly worse in LOCO (19.3 → 20.8 %, t +0.9), but worse when
whole configurations are held out from the members (productivity scheme,
v1 protocol: 14.1 → 16.6 %, t +5.1, mostly corn ethanol and IBO). The
same loss appears when the electricity head merely loses its productivity
inputs without any structure (14.1 → 16.1 %, t +3.7), and adding the
product groups to the electricity inputs does not recover it (16.3 %):
r (whose baseline differs per configuration, 0.12–13 g/L/h) acted as a
configuration proxy, not as physics. A learned β (bounds 0.5–1.5) stays at
1.01 ± 0.01 and gives the same errors as β = 1. The v1-protocol
productivity scheme holds out ~30 % of the configurations from every member
(calibration + early-stopping groups by configuration), so its error is
mostly that of unseen configurations (5 % with sample-level hold-outs vs 15
%); r-extrapolation is measured by the sample-level variant (the v2
protocol of §6 uses sample-level hold-outs in every scheme).

**Extrapolation safeguards** (2026-10-17, round 4; library defaults off,
`experiments.py` default on, CLI `--no-extrapolation-safeguards`,
`experiments.EXTRAPOLATION_SAFEGUARDS`; all need `FeatureSpec(downstream=True)`,
which appends the raw structural inputs below after `struct_water_removal`;
physics model only, the plain MLP is unchanged; the sklearn baselines see the
two new structural columns in their design matrix: their development fold
means change by ≤ 0.6 (ridge), 2.4 (HGB) and 0.0 (kNN) points under LOFO, LOPO
and 5-fold LOCO):

* *Sugar-availability prior* (`ModelConfig.fs_prior`): the structured
  feedstock intensity becomes `b'_fs = exp(z)/w_sugar` with
  `z = log(1/s_avail) + Δz`, where `s_avail = struct_sugar_availability =
  (free sugars + starch + s · structural polysaccharides)/w_sugar`
  (glucose-eq., `features.sugar_availability`; `s` = `p_saccharification_yield`,
  default 0.9, for a lignocellulosic feedstock, 0 otherwise: sugarcane bagasse
  is burned) and the network predicts only the residual `Δz` (centered on
  its training median `fs_prior_center`; `fs_eta_bound` does not apply). On
  the development set the prior is 1 (corn, glucose), 1/0.90 (corn stover)
  and 1/0.59 (sugarcane) against configuration medians of `b'_fs·w_sugar` of
  0.93–1.01, 1.04–1.19 and 1.71–1.80; the residual's spread (sd 0.12, range
  −0.8 … +0.07, the low end lactic's b′/d′ split) is 40 % of that of `z` (0.29).
  Without it, `z` of a new feedstock is whatever the feedstock sub-network
  extrapolates (corn: 1/η = 1.50 for every corn sample, truth 0.93–1.01).
* *Feedstock-input clipping* (`ModelConfig.clip_groups`, here `('feedstock',
  'pretreatment')`): those groups' continuous columns are clipped to their
  training range (margin `max(clip_inputs, 0)`), independently of the HPO
  dimension `clip_inputs`. Minor mass fractions have a tiny training spread
  (`fs_w_other` ∈ {0, 0.0145, 0.015}, sd 0.006), so corn (0.10–0.13) enters
  at +18 sd (+9.5 sd beyond the training maximum after the 10-sd clip of the
  FeatureSpec) and saturates every sub-network that sees the feedstock group
  (11 of the 15 outputs). Clipping alone does not fix the feedstock intensity
  (corn 1/η̂ 1.30 vs truth 0.98, sugarcane 0.97 vs 1.76–1.80; 3 seeds), the prior alone neither (one seed: the saturated residual gives
  1/η̂ 0.77); together corn gets 1.03 and sugarcane 1.63 (§6, round 4).
* *Neutralization gate* (`ModelConfig.neutralization_gate`): the whole
  neutralization channel is multiplied by `struct_neutralization` (0 when the
  fermentation product has no acid group and no neutralizing / pH-control
  agent or neutralization degree is given: nothing to neutralize or
  acidulate; `features.neutralization_factor`; unknown → 1). Data check
  (`check_constraints` row `neutralization_gate`, the gate is disabled when
  the data violate it): the channel is exactly 0 for all 399 gate-0
  development samples (ethanol, isobutanol) and nonzero for TAL (NaOH pH
  control, no acid group), HP low-pH and neutralized lactic (low-pH lactic, an
  acid without agent, is 0 but ungated). Without it, a held-out
  alcohol inherits TAL's or HP's neutralization cost (ethanol LOPO: +13–18 %
  of the true numerator).

**Power gate** (`ModelConfig.power_gate`, round 5; library default off;
needs `FeatureSpec(downstream=True, power_export=True)`, which
`train.train_ensemble` requests itself when the gate is on): the structural
input `struct_power_export` (`features.power_export_factor`) is 1 if
'electricity' is among the coproducts (an on-site boiler-turbogenerator can
export surplus power, so the net purchased electricity may be negative), 0 if
the coproducts are known without it (no generation, nothing to export), and
unknown → 1 (ungated). Where it is 0, every coefficient of the net
electricity intensity is mapped to `s·softplus(q/s) ≥ 0` (`s` =
`POWER_GATE_SCALE` = 0.1 × the coefficient's training sd), so the purchased
electricity and ∂MPSP/∂p_el are ≥ 0 at every y, t > 0. The map acts after the
productivity structure (a monotone increasing map of a function decreasing in
r keeps ∂MPSP/∂r ≤ 0) and before the price (exact price linearity); the
reported `r_theta0` / `r_theta1` are before the gate. Rows with
`struct_power_export` = 1 are unchanged (identical outputs to an ungated
twin, tested). Data check (`check_constraints` row `power_gate`; the gate is
disabled when the data violate it): on the development set the net
electricity intensity is ≥ 0 at all 3667 points of the 177 samples without an
electricity coproduct (ethanol corn, three lactic configurations; minimum
0.086 kWh/kg), and a per-coefficient non-negative fit (NNLS) of their
per-point intensities costs at most 0.09 % of MPSP (OLS 0.03 %).

## 4. Loss

    L = w_mpsp · mean_i [(MPSP̂_i − MPSP_i)/MPSP_i]²       (equation layer; simulated points)
      + w_tot  · Σ_k [(θ̂_k − θ_k)/s_k]²                    (total coefficients a',b',c',d')
      + w_ch   · Σ_c Σ_k [(θ̂_ck − θ_ck)/s_ck]²             (channel coefficients — intermediate data)
      + w_g    · [(ĝ − g)/s_g]²
      + w_phys · (physics residuals at labelled AND collocation points)
      + weight decay

Physics residuals (soft constraints; each enabled only if the training data
satisfy it, which is checked and reported):

1. ∂MPSP/∂y ≤ 0 and ∂MPSP/∂t ≤ 0 over the TY domain
   (`b' + d'/t ≥ 0`, `c' + d'/y ≥ 0` at domain grid points);
2. MPSP > 0 over the domain;
3. ∂MPSP/∂r ≤ 0 (higher productivity → smaller fermenters) via autograd
   w.r.t. the productivity input (identically satisfied with the hard
   productivity structure of §3, `ModelConfig.r_structure`);
4. non-negative feedstock and natural-gas intensities over the domain.

**Collocation points**: unlabeled inputs created by perturbing training inputs
(productivity log-uniform over a range wider than the data; the prices can be
widened by `LossConfig.colloc_price_range`, but the pipeline and the release
keep them inside the training box, `[1, 1]`); only physics residuals 1–4 are
evaluated there (PINN-style). Residuals the training data violate in more
than 1 % of the samples are disabled (`train.check_constraints`): in the
release only positive MPSP and the non-negative feedstock intensity are
active. Monotonicity in y is violated by 7 % of the samples, all of them
isobutanol, whose yield deficit is sold as ethanol, and monotonicity in t by
2 %, all of them oxalic acid (§1, "Product-specific cost structure"), so
both are off for every product; monotonicity in r is replaced by the hard
structure of §3; the natural-gas non-negativity is disabled because the
4-term fits of the hinge-shaped natural-gas demand go negative inside the
domain in 58 % of the samples although the simulated use never does (§1),
so a predicted natural-gas intensity can be negative where the fitted
targets are (`models/release/meta.json` → `constraint_report`,
`loss_config`).

**Settings of the v3 pipeline and the shipped checkpoints**
(`models/release/meta.json` → `loss_config`): channel term in surface space
(`ch_space='surface'`) with the total in coefficient space
(`coef_space='coef'`), `w_ch = 3` (fixed, chosen on out-of-fold channel
fidelity, §6 round 5), `w_mpsp = 0.3`, `w_tot = w_g = 1`, `w_phys = 0.03`,
collocation ratio 2, no configuration balancing; the TAL / TAL_SA
feedstock-channel targets are masked (D8). The formula above shows the
coefficient-space channel term, the default of `LossConfig`.

Options (`LossConfig`; defaults = the loss above):

* `coef_space='surface'` compares the total and every channel/intensity in
  *surface* space instead of coefficient space: at each simulated point of
  the sample, `(N̂ − N)/N` for the total and `Σ_c ((m̂_c − m_c)/N)²` for the
  channels (price channels as price × intensity), with `N` = target total
  numerator `g·MPSP_fit`. This is well conditioned (the per-sample OLS
  estimates of a', b', c', d' are strongly correlated, §6), but in the
  development CV it made in-family configurations worse (§6), so it is off;
* `ch_space='surface'` (round 5; default None = `coef_space`) moves only the
  channel / intensity term to surface space and keeps the total in
  `coef_space`: per sample Σ_c mean_i ((m̂_c − m_c)/N)² over its simulated
  points, price channels as price × intensity. This is the squared
  out-of-fold channel-fidelity metric itself (`evaluate.channel_fidelity`,
  §6 round 5), so every channel weighs by its economic size at the sample's
  own TY points, instead of by its cross-configuration coefficient spread
  (in coefficient space an electricity error of 25 % of MPSP for corn
  ethanol, carried by a′, is 1.6 training sd of one of the 56 standardized
  entries: 2.6 / 56 × w_ch = 0.005 of loss at the v2 weight 0.1, while
  other channels can absorb it without changing the total). Its `w_ch` is on
  a different scale from the coefficient-space one;
* `train.mask_channel_targets(samples, {channel: [biorefinery | config_id |
  product, ...]})` sets the channel's targets (a price channel: its
  intensity targets) of those rows to NaN, i.e. masks them in the loss and
  the target scales; totals, g and the other channels keep theirs (for
  known defective decompositions, §1);
* `config_balance=α` weights every sample by `n_config^(−α)` (batch-normalized),
  α = 1 → every configuration counts once (off by default; not better in CV).

## 5. Features (specifications → vector)

* feedstock: type (pure sugar/starch/sugar juice/lignocellulose), composition
  mass fractions (fermentable sugar glucose-eq., lignin, ash, protein, lipid,
  water, …), price, capacity;
* fermentation product: MW, C/H/O/N counts, degree of reduction, theoretical
  yield, Tb, Tm, ΔHvap, Psat at fermentation T, infinite-dilution relative
  volatility vs. water (UNIFAC via thermosteam), aqueous solubility, number of
  acid groups, pKa1, LHV;
* final product: MW, stoichiometric mass ratio to fermentation product,
  upgrading flag, Tb, Tm, LHV;
* final product and downstream process (structured, `FeatureSpec(downstream=True)`;
  replaces the free-text `separation_principle`, see below);
* regime: neutralization flag and agent, aerobic flag, fermentation mode,
  temperature, productivity (log r and 1/r), inflection/theoretical yield ratio,
  separation principle;
* contextual (optional, baseline-imputed with missing-indicators): electricity
  price, natural gas price, IRR, income tax, operating days, boiler/
  turbogenerator efficiency, etc.

### Downstream specification (`downstream.py`)

`separation_principle` is free text, i.e. the product identity in disguise,
and the final product's water content and the separation route -- which set
the titer-dependent c′ (and part of a′) -- were not specified at all. Per kg
product the broth water to remove is ≈ 1000/t − (1 − w_f)/w_f (t in g/L, w_f
the final-product mass fraction; recorded broth water per kg product ∝ 1/t
with slope 1000 for ethanol, succinic, oxalic, IBO), a c′ term 1000 × [water
removed] and an a′ term −(1 − w_f)/w_f × [water removed], times a unit cost
that depends on the route. Fields (units/vocabularies in `downstream.FIELDS`,
merged into `features.FIELD_UNITS`; multi-valued fields '+'-joined; empty =
'none', which survives CSV round trips unlike ''):

| field | values |
|---|---|
| `final_product_mass_fraction` | w_f at baseline (purity of the sold stream; MPSP is per kg pure product); None for a broth product |
| `final_product_is_broth` | the broth is sold (w_f ≈ t/1000: 0.005 / 0.031 / 0.117 at 5 / 30 / 100 g/L for oxalic) |
| `final_product_phase` | liquid, solid, aqueous_solution |
| `water_removal_route` | evaporate_water, crystallize, distill_product, none |
| `separation_ops` | centrifugation, filtration, acidulation, gypsum_filtration, adsorption, ion_exchange, membrane, flash, evaporation, distillation, LLE, decantation, molecular_sieves, crystallization, drying, esterification, hydrolysis |
| `upgrading_ops` | dehydration, hydrogenation, ring_opening_hydrolysis, esterification, hydrolysis, decarboxylation, oxidation, amination, polymerization |
| `coproducts` | electricity (on-site turbogenerator generating electricity), DDGS, corn_oil, ethanol, lignin, gypsum, CO2, fertilizer |
| `neutralization_degree` | mol base-eq. per mol product acid groups converted to the salt (0–1) |
| `acid_regeneration` | H2SO4, none |

Values for all 41 planned configurations (`downstream.TABLE`, sources in its
`downstream_sources` field): from the adapters (separation principle /
variant, neutralization, lactic `product_purity_spec` 0.88, oxalic
`product_form` and `neutralization_degree` 0.91, succinic partial pH control
0.42 mol OH-eq./mol → degree 0.21, isobutanol `coproducts`), the
Bioindustrial-Park flowsheets (read-only: unit operations of the separation
and upgrading areas, boiler-turbogenerator setup) and the recorded baseline
data (w_f = median purity of the final-product stream; electricity listed
where the CHP generates at baseline):

| family | w_f | phase | route | separation ops | upgrading | coproducts | neutr. degree / acid |
|---|---|---|---|---|---|---|---|
| HP (improved) | 0.9967 | liquid | evaporate_water | centrifugation, adsorption, ion_exchange, evaporation, distillation (+ acidulation, gypsum_filtration if neutral) | dehydration | electricity | 1 / H2SO4 (neutral), 0 / none |
| HP hexanol | 0.9967 | liquid | evaporate_water | centrifugation, evaporation, LLE, distillation (+ acid./gypsum if neutral) | dehydration | electricity | as above |
| TAL | 0.82–0.90 (per feedstock; varies with titer) | solid | crystallize | centrifugation, evaporation, crystallization, drying, distillation | – | electricity | 0 / none |
| TAL_SA (KSA) | 0.949 | solid | crystallize | as TAL | hydrogenation, dehydration, ring_opening_hydrolysis | electricity | 0 / none |
| succinic | 0.8487 | solid | crystallize | centrifugation, evaporation, crystallization, filtration, drying (+ acid./gypsum if neutral) | – | electricity | 0.21 / none; neutral 1 / H2SO4 |
| lactic | 0.88 | aqueous_solution | evaporate_water | filtration, flash, distillation, esterification, hydrolysis (+ acid./gypsum if neutral) | – | – | 1 / H2SO4 (neutral), 0 / none |
| ethanol corn / cornstover / sugarcane | 0.9505 / 0.9693 / 0.9704 | liquid | distill_product | distillation, molecular_sieves (+ corn: centrifugation, evaporation, drying [DDGS]; cornstover: filtration) | – | corn: DDGS, corn_oil; others: electricity | 0 / none |
| oxalic | broth | aqueous_solution | none | centrifugation, acidulation, gypsum_filtration | – | electricity | 0.91 / H2SO4 |
| IBO_EtOH_corn | 1.0 | liquid | distill_product | centrifugation, decantation, distillation, drying, evaporation, molecular_sieves | – | DDGS, corn_oil, electricity, ethanol | 0 / none |

succinic_neutral: the data of the v2 models did not contain the acidulation
the table describes (`acid_cost` 0, R401 / S405 bypassed; D15). Since
succinic adapter v3 the lime forms calcium succinate, R401 acidulates it with
H2SO4 and S405 filters the gypsum; the four configurations were regenerated
(baseline product purity unchanged, 0.8486), so table and data agree for
models trained on the regenerated data.

Derived features (`features.DS_PHYS`): `phys_water_removal_c` = 1000 ×
[water removed], `phys_water_removal_a` = −(1 − w_f)/w_f × [water removed]
(0 for a broth product), `phys_evap_water_c` = 1000 × [evaporate_water or
crystallize], `phys_distill_c` = 1000 × [distill_product],
`phys_log10_w_final`, neutralization stoichiometry per kg fermentation
product `phys_base_per_product` = n_acid × degree × (g base per eq.) / MW,
`phys_acid_per_product` (H2SO4) and `phys_gypsum_per_product` (Ca base +
H2SO4), `phys_n_separation_ops`, `phys_n_upgrading_ops`; plus multi-hot
`ds_sep__*`, `ds_upg__*`, `ds_cop__*`, `neutralization_degree`,
`final_product_is_broth`, one-hot `final_product_phase`,
`water_removal_route`, and the raw structural input `struct_water_removal`
(model gate, §3). With `downstream=False` (default) none of these are used and
the feature matrix is exactly the previous one.

Descriptors of known feedstocks/products are stored in a registry so users can
specify them by name (composition and baseline capacity per feedstock; yields,
acid properties per product; per product pair the stoichiometric ratio, the
fermentation regime where it is unique in training — otherwise the user is
asked to specify it — and the baseline value of every contextual `p_*`
parameter, or "not applicable" which keeps its missing indicator at 1; the
downstream fields are regime fields, and `config_id` of a known configuration
fills them from `downstream.TABLE`); new
chemicals can be described with thermosteam or by supplying descriptor values.
Omitted categoricals of a query are encoded as the unknown level (trained by
categorical dropout), and every imputed / defaulted field is reported per row.
Input units and bases are listed in `features.FIELD_UNITS` (feedstock prices
per kg WET feedstock). The applicability domain (kNN + range) uses encoder
features and productivity only: price-explicit channels are exactly linear in
price, so prices outside the training range are reported separately.

## 6. Training & validation protocol

* Deep ensemble (M members, different seeds/bootstraps) → mean ± epistemic
  spread; split-conformal intervals scale the spread, one nonconformity
  score per sample (|target − mean| / (std + eps); MPSP: max over its points,
  so coverage is simultaneous over the sample's surface); if there are too
  few scores for the level the quantile is +inf (no interval is reported).
* **Training protocol of the pipeline (v2, `experiments.PROTOCOL`; round 4
  below):** every network fit -- HPO, every CV scheme (PI-NN, plain MLP,
  ablations), `models/final_dev`, `models/release` -- uses the same
  procedure, so the CV estimates the training procedure of the final models
  with the selected hyperparameters (the HPO chose them on the same
  development set and is not nested in the CV, so the CV is slightly
  optimistic for HPO + fit): no calibration hold-out (`calib_frac` 0); every member early-stops on its
  own random 15 % *sample-level* hold-out (`group_col` None; the hold-outs
  differ between members, so the ensemble sees every sample; patience 40,
  20–300 epochs); no refit (`TrainConfig.refit` 'none'). The baselines fit
  exactly the rows the members train on. CV, final_dev and release ensembles
  have the same number of members (`--members`, 5), because the scores divide
  by the member spread.
* **Intervals of the final models** are calibrated from the out-of-fold
  scores of the development leave-one-configuration-out CV of that same
  procedure (`train.conformal_scores` per held-out sample,
  `train.calibrate_conformal_oof`): one score per development sample, from
  the fold ensemble that did not see its configuration; `experiments.
  oof_calibration_scores` refuses a score table that misses or repeats a
  development sample or contains a test sample (for the PI-NN this stops a
  full run; an unusable plain-MLP table only leaves the plain `final_dev`
  without intervals). The release ensemble (dev + test samples) uses the
  same development scores. The scores of a′–d′ and g divide by (member
  spread + eps), eps = 0.1 × the training target sd; the score tables keep
  the raw |error| and spread (`<key>_abserr`, `<key>_sd`), and
  `calibrate_conformal_oof` re-scores them with the final ensemble's own eps
  (`train.rescore_conformal`; `conformal['eps_rescored']`), so that the
  quantile and the interval use the same floor (the LOCO-fold / full-dev
  target-sd ratio is 0.74–1.08; MPSP has a fixed 1 % relative floor).
  `ens.conformal` records the source, n, the number of configurations, the
  weighting ('sample': the plain quantile over all scores, configurations
  count by their numbers of samples, so the coverage is marginal over the
  samples of new configurations; 'group' = every configuration weighted
  equally, pooled-CDF quantile of Dunn et al. 2022, available),
  `coverage_guaranteed = False` (the fold ensembles are not the final
  ensemble; `<key>_level` is then only nominal, `level_kind`) and a
  cross-fitted coverage estimate (`crossfit`: quantile of all other LOCO
  folds applied to each fold; mean = fold-balanced, pooled = over samples,
  min, q10; compare pooled with the realized pooled test coverage and mean
  with the realized mean over test configurations). The random-CV calibration is stored for reporting only
  (`meta['conformal_random_cv']`). The test split is never used to calibrate.
* v1 protocol (`--v1-protocol`, `experiments.PROTOCOL_V1`, commit fd3f386):
  calibration on a held-out 15 % of whole groups and early stopping on
  held-out groups of the same kind (`experiments.INNER_GROUP_V1`: configurations
  for LOCO, LOFO and productivity, fermentation product for LOPO, samples for
  random CV; the library default `train.SCHEME_GROUP` groups like the outer
  scheme); `final_dev` with a configuration-grouped 15 % calibration hold-out.
  `train.calibrate_conformal` (one hold-out) and `train_ensemble(calib_frac,
  group_col)` keep supporting it; `coverage_guaranteed` there needs ≥ 10
  calibration groups.
* Splits are **grouped** (never split TY points of one sample; MC samples of a
  configuration are correlated):
  - random sample-level K-fold (interpolation);
  - leave-one-configuration-out;
  - leave-one-feedstock-out;
  - leave-one-product-out (hardest; generalization to new products);
  - productivity extrapolation (train on the central productivity range).
* Hyperparameters chosen on inner validation only; a final test set that is
  never used to fit, tune, select or calibrate, reported once per pipeline
  run. Its first use was the v1 run (fd3f386); the v1 test failures then
  directed which failure modes round 4 investigated (development data
  only), so the v2 run re-uses it adaptively: the v2 test numbers and the
  v1 → v2 test change are likely optimistic, not an independent
  confirmation (the development CV is the primary evidence; v2 was fixed
  before its test evaluation and is released whatever that comparison
  shows). The v3 run (2060257) is its third use: the review that found the
  v2 → v3 data defects and set the round-5 agenda had read the v2 report,
  test results included, so the v3 test numbers are likewise optimistic;
  the data corrections also replaced the test rows of the regenerated
  configurations (REPORT.md notice). It is **frozen** in
  `results/split.json` (committed 2026-10-16, before any modelling) from the
  41 *planned* configurations (`experiments.PLANNED_CONFIGS`, = the adapters'
  `configurations()`) and the seed 2026 only, never from the data that
  happen to exist: one whole configuration per biorefinery with ≥ 4 planned
  configurations (HP_corn, TAL_sugarcane, succinic_glucose,
  lactic_neutral_SSCF_cornstover; out-of-distribution) + every sample of the
  other configurations whose seeded sample-id hash is < 0.10
  (in-distribution; membership is a function of the sample id, so stable as
  campaigns add samples). `experiments.py` loads the file and the final
  stage fails if a test configuration has no data.
* HPO (`experiments.py`): `hpo_members` = 2 networks per candidate and fold
  (full mode), trained with the pipeline protocol above (v2: sample-level
  early stopping, so the selection is for the procedure actually used),
  selection by the fold mean of the **configuration-balanced**
  held-out MPSP MAPE (`mpsp_mape_cfg`: mean over held-out configurations of
  their point MAPE; `compute_metrics` also reports its geometric mean
  `mpsp_mape_cfg_gmean`). Reasons (development diagnostics below): the
  point-weighted fold MAPE is dominated by the fold's largest configuration,
  and single-member results differ between seeds by up to ~50 % (relative).
  Round 5 added the option (1 − w)·(that LOCO fold mean) + w·(mean over the
  seven leave-one-product-out folds of the same statistic)
  (`--hpo-lopo-weight w`; the LOPO CV estimates of such a run are then
  post-selection, like LOCO's), but the pipeline keeps w =
  `HPO_LOPO_WEIGHT` = 0, the v2 criterion: the clip exemption removes the
  failure that motivated it, after which w = 0.5 selects on seed noise
  (5e-5f). Under the round-5 channel loss `w_ch` is fixed (`CHANNEL_LOSS`)
  and no longer searched.
* Baselines: nearest-configuration lookup, ridge, gradient boosting, plain MLP
  (coefficients only); ablations removing each physics component. Baselines
  are fit on exactly the samples the networks train on (same calibration and
  early-stopping hold-outs); results report `n_fit`, `n_fit_configs`,
  `n_calib`, and summaries average only folds in which every method
  succeeded (`n_folds_ok`).
* Collocation draws prices inside the training price box (where the
  price-dependent residuals were validated on the data).
* Metrics: per-coefficient MAE/R² (natural units), MPSP-surface MAPE on held
  out simulated TY points (vs. the irreducible error of the per-sample OLS
  fit), physics-constraint violation rates, interval coverage. Additions
  (2026-10-17):
  - MAPE split by sample kind (`*_spec`, `*_random`, also configuration-
    balanced): *spec-determined* samples (`sample_kind` baseline /
    productivity: every biorefinery-specific `bp_*` parameter at baseline,
    so the inputs determine the target) vs *random* samples (unobserved
    `bp_*` draws, which the inputs do not describe);
  - `evaluate.noise_floor`: within each configuration (≥ 10 random
    samples), a 5-fold ridge on the observable inputs (`p_*`, prices, log r,
    1/r) predicts a held-out random sample's coefficient vector (also
    reported: + `bp_*` inputs, and the configuration mean). It was meant as
    an estimate of the error the unobserved `bp_*` cause, but it is NOT a
    noise floor: it is the model-form error of a linear model on 20-60
    samples per configuration (v3: 6.3 % configuration geo-mean, while the
    PI-NN predicts random samples at 3.3 % in random CV, and adding the
    `bp_*` draws makes the ridge worse, 7.2 %; REPORT.md section 1);
  - `coef_gram_rel`: conditioning-aware coefficient error, the numerator
    error vector e in the norm of the sample's own design, √(eᵀGe) with
    G = ΦᵀΦ/n over its points (= RMS error of the numerator surface),
    relative to the RMS target numerator; per coefficient √G_kk·|e_k|.
    Errors along the ill-determined b′–d′ direction, which leave the surface
    unchanged, are not penalized (unlike MAE/R²).

### Development-set diagnostics (2026-10-16; never on test data)

Data: 1102 samples / 29 configurations; development set 911 / 25 (test split
removed). 5-fold leave-configurations-out (folds stratified by biorefinery),
one network per fit, default hyperparameters, seed 0:

| method | fold-mean MAPE ± sd | config-mean | config median | config geo-mean |
|---|---|---|---|---|
| PI-NN | 61 ± 58 % | 28.5 % | 13.7 % | 16.3 % |
| plain MLP | 62 ± 40 % | 30.0 % | 14.0 % | 16.4 % |
| ridge | 124 ± 183 % | 47.5 % | 13.0 % | 16.3 % |
| HGB | 103 ± 145 % | 41.9 % | 16.8 % | 17.6 % |
| kNN | 65 ± 60 % | 38.5 % | 17.3 % | 25.0 % |

The per-sample fit floor is 1.5 % (`fit_mape`). Findings:

* Errors are dominated by held-out configurations without a close relative:
  oxalic (the only fermentation product of its kind; 189 samples; PI-NN 216 %,
  predicted c′ ≈ 20–25 vs 4.6), IBO (89 %, a′ ≈ 0 vs 0.93) and lactic (66 %,
  c′ ≈ 6 vs 28, g 1.0 vs 0.63). Configurations within HP/TAL/succinic are
  typically at 4–20 %. Leave-one-product-out on the families with ≥ 3
  configurations (PI-NN / HGB / kNN / ridge): succinic 20 / 51 / 67 / 115 %,
  TAL 64 / 27 / 52 / 68 %, HP 28 / 16 / 33 / 69 %, ethanol 114 / 92 / 137 / 281 %.
* Conditioning: on the TY domain the basis [1, 1/y, 1/t, 1/(yt)] is
  collinear (column-normalized condition number median 10, up to 53;
  lactic, whose domain starts at y = 0.3, t = 40 g/L: 35). The per-sample OLS
  estimates are strongly correlated (median correlation b′–d′ −0.90, a′–b′
  −0.79, a′–c′ −0.71), and lactic's d′ reaches 45 with c′ down to −5.4.
  Out-of-fold R² (PI-NN): a′ 0.55, b′ 0.31, c′ −0.09, d′ −1.11, g 0.59. Models
  that predict coefficients independently (ridge, HGB) combine them into
  surfaces that are badly off (oxalic fold: 449 % / 356 %).
* Channels: RMS error of each predicted channel relative to the total
  numerator: coproducts 0.47, electricity 0.26, materials_other 0.20,
  capex_utilities 0.18, natural gas 0.16, feedstock 0.14, the rest < 0.09.
* Seed variance: a second seed (different initialization and early-stopping
  hold-out) gave fold MAPEs of 64 / 20 / 208 / 98 % instead of 41 / 21 / 160 /
  63 % (folds 1–4). Averaging the two seeds (a 2-member ensemble) landed in
  between (fold 3: 183 %). The large errors are systematic extrapolation bias,
  not variance.
* Variants (PI-NN, same folds and seeds; fold-mean / config geo-mean):
  baseline 61 / 16.3 %; `coef_space='surface'` 52 / 20.9 % (better on oxalic,
  103 %, worse in-family, stops early: about 83 epochs); w_ch = 0 70 / 16.5 %
  (channel supervision in coefficient space helps slightly or is neutral);
  `config_balance=0.5` + early stop on MPSP 76 / 17.9 %; `clip_inputs=0`
  55 / 16.2 % (oxalic 173 %), but under leave-one-product-out it was worse for
  HP (63 % vs 28 %) and ethanol (155 % vs 114 %). No variant is clearly better
  given the seed variance, so the defaults are unchanged. Clipping is an HPO
  dimension.

### Development-set diagnostics, round 2 (2026-10-17; never on test data)

Same development tables (911 samples / 25 configurations; the patched tables
add only the S1 fields, and `downstream=False` reproduces the previous
round's fold results exactly), same folds and protocol (5-fold
leave-configurations-out stratified by biorefinery, one network per fit,
default hyperparameters), now 3 seeds (0, 1, 2) per fit for every method.
"fold" = fold-mean point MAPE; cfg = configuration-balanced; spec / random =
configuration geo-means over spec-determined / random samples; mean ± sd over
seeds:

| method | features | fold | cfg mean | cfg geo-mean | spec | random | oxalic | IBO |
|---|---|---|---|---|---|---|---|---|
| PI-NN | without | 54 ± 8 | 30.8 ± 4.3 | 20.7 ± 4.6 | 19.2 | 22.2 | 131 ± 74 | 110 ± 29 |
| PI-NN | downstream | 51 ± 10 | 30.8 ± 7.6 | 21.3 ± 6.1 | 19.7 | 23.6 | 89 ± 56 | 101 ± 35 |
| PI-NN | downstream + gate | 43 ± 5 | 26.9 ± 2.3 | 18.9 ± 2.5 | 17.6 | 21.1 | 53 ± 37 | 74 ± 17 |
| plain MLP | without | 79 ± 21 | 37.5 ± 8.9 | 19.7 ± 3.9 | 17.7 | 21.8 | 138 ± 25 | 195 ± 97 |
| plain MLP | downstream | 47 ± 11 | 26.4 ± 5.0 | 18.0 ± 3.0 | 16.2 | 19.4 | 95 ± 18 | 60 ± 23 |
| ridge | without | 109 ± 35 | 44.1 ± 8.5 | 16.8 ± 1.0 | 13.1 | 19.2 | 460 ± 361 | 113 ± 32 |
| ridge | downstream | 79 ± 9 | 35.9 ± 3.2 | 17.6 ± 2.2 | 14.8 | 21.7 | 134 ± 79 | 135 ± 55 |
| HGB | without | 87 ± 39 | 37.2 ± 10.4 | 17.6 ± 0.6 | 13.5 | 19.3 | 367 ± 280 | 133 ± 19 |
| HGB | downstream | 86 ± 42 | 37.1 ± 11.1 | 18.1 ± 1.0 | 14.3 | 20.3 | 347 ± 264 | 126 ± 15 |
| kNN | without | 58 ± 24 | 35.4 ± 8.7 | 23.4 ± 2.6 | 23.2 | 24.2 | 49 | 34 |
| kNN | downstream | 40 ± 17 | 24.7 ± 6.3 | 18.7 ± 1.7 | 18.5 | 19.4 | 49 | 34 |

Spread of the seed-averaged fold MAPEs across folds (sd): PI-NN 34 → 24 → 19
(gate), plain 52 → 24, ridge 131 → 23, HGB 107 → 100, kNN 49 → 15. Gram-norm
coefficient error (configuration geo-mean): PI-NN 17.0 / 19.0 / 15.4 % (gate),
plain 14.6 → 13.6, ridge 13.2 → 14.7, HGB 16.7 → 16.0, kNN 18.3 → 15.4 --
while the per-coefficient out-of-fold R² of d′ is −1 to −15, i.e. most of the
coefficient error lies along the ill-determined b′–d′ direction.

Leave-one-product-out (families with ≥ 3 configurations; mean over seeds,
mean of the 4 families ± seed sd): PI-NN 62 ± 5 → 50 ± 14 → 49 ± 12 % (gate;
ethanol 113 → 88 → 80, HP 34 → 26 → 25, succinic 30 → 23 → 28, TAL 70 → 63 →
61), plain 70 ± 17 → 62 ± 12, ridge 222 ± 168 → 109 ± 52, HGB 60 ± 15 → 66 ±
15 (TAL 37 → 67), kNN 77 ± 12 → 71 ± 12.

Noise floor (16 configurations): random samples 6.1 % (geo-mean; mean 7.2 %),
configuration mean 10.7 %, 4-coefficient fit 1.4 %. Spec-determined samples
are only ~3 points easier than random ones for every method, so the
remaining cross-configuration error (≈ 3× the floor) is configuration-level
bias, not noise.

Findings: the structured fields help the models that generalize through
feature similarity (plain MLP, kNN, ridge; oxalic 460 → 134 % for ridge, IBO
195 → 60 % for the MLP); trees cannot use a route never seen in training
(HGB unchanged). For the PI-NN they help only together with the
water-removal gate, which enforces the vanishing broth-water costs of the
broth product (oxalic 131 → 53 %) and halves the seed spread. A compact set
without the per-operation multi-hot columns was worse for every method (PI-NN
fold 65 %, plain 77 %, gate 55 %). The experiment pipeline therefore uses the
downstream features and the gate by default (`--no-downstream`,
`--no-water-gate`); the library defaults (`FeatureSpec(downstream=False)`,
`ModelConfig(water_gate=False)`) are unchanged. Seed 0 of the previous round
was favourable: its single-seed PI-NN geo-mean 16.3 % compares with 20.7 ±
4.6 % over three seeds.

### Development-set diagnostics, round 3: channel-wise input masking (2026-10-17; never on test data)

Same development tables, folds and protocol as round 2 (downstream features;
PI-NNs with the water gate; one network per fit, default hyperparameters,
masked sub-networks of width 64 and depth 2), seeds 0–4 for the shared and
the masked PI-NN, 0–2 otherwise; leave-one-product-out now over all seven
families (inner validation by configuration). Leave-configurations-out
(5-fold, stratified by biorefinery), mean ± sd over seeds:

| method | seeds | fold | cfg mean | cfg geo-mean | spec | random | oxalic | IBO | lactic | ethanol |
|---|---|---|---|---|---|---|---|---|---|---|
| PI-NN shared encoder | 5 | 48.6 ± 8.5 | 29.5 ± 3.9 | 21.2 ± 3.8 | 19.8 | 22.6 | 61 ± 31 | 92 ± 31 | 56 ± 12 | 59 ± 16 |
| PI-NN masked ('physics') | 5 | 45.7 ± 5.4 | 27.9 ± 5.8 | 18.4 ± 6.2 | 16.7 | 19.5 | 65 ± 14 | 96 ± 26 | 56 ± 7 | 51 ± 9 |
| shared, seeds 0–2 | 3 | 42.9 ± 4.6 | 26.9 ± 2.3 | 18.9 ± 2.5 | 17.6 | 21.1 | 53 ± 37 | 74 ± 17 | 63 ± 2 | 51 ± 17 |
| masked, seeds 0–2 | 3 | 43.1 ± 5.6 | 24.2 ± 3.7 | 14.5 ± 1.8 | 12.6 | 15.8 | 68 ± 16 | 93 ± 29 | 57 ± 9 | 48 ± 9 |
| hybrid, residual weight decay 1 | 3 | 46.8 ± 2.3 | 25.7 ± 1.8 | 15.1 ± 0.6 | 13.4 | 16.7 | 65 ± 20 | 104 ± 14 | 56 ± 9 | 56 ± 12 |
| hybrid, residual weight decay 10 | 3 | 43.5 ± 3.5 | 24.3 ± 3.2 | 14.4 ± 1.4 | 12.1 | 15.5 | 66 ± 17 | 94 ± 3 | 55 ± 10 | 51 ± 8 |
| plain MLP | 3 | 47.1 ± 11.4 | 26.4 ± 5.0 | 18.0 ± 3.0 | 16.2 | 19.4 | 95 ± 18 | 60 ± 23 | 70 ± 3 | 55 ± 25 |
| kNN | 3 | 39.8 ± 17.4 | 24.7 ± 6.3 | 18.7 ± 1.7 | 18.5 | 19.4 | 49 | 34 | 59 ± 15 | 61 ± 51 |

Masked − shared, paired by seed (same folds and inner splits): fold −2.9
(t = −0.8), configuration mean −1.7 (t = −1.5), geo-mean −2.8 points (t =
−1.5; better on 4 of 5 seeds, seed 3 is an outlier with TAL_SA at 60–70 %).
Leave-one-product-out (point MAPE of the held-out family):

| method | seeds | ethanol | HP | IBO | lactic | oxalic | succinic | TAL | mean of 4 | mean of 7 |
|---|---|---|---|---|---|---|---|---|---|---|
| PI-NN shared encoder | 5 | 74 ± 37 | 24 ± 5 | 86 ± 46 | 56 ± 6 | 74 ± 40 | 33 ± 15 | 61 ± 3 | 48.2 ± 10.0 | 58.3 ± 10.6 |
| PI-NN masked | 5 | 83 ± 10 | 15 ± 7 | 96 ± 30 | 54 ± 6 | 79 ± 27 | 35 ± 15 | 42 ± 2 | 43.6 ± 7.1 | 57.5 ± 9.0 |
| hybrid, weight decay 1 | 3 | 104 ± 19 | 15 ± 2 | 81 ± 4 | 53 ± 11 | 76 ± 20 | 32 ± 15 | 48 ± 2 | 49.9 ± 3.0 | 58.6 ± 3.9 |
| hybrid, weight decay 10 | 3 | 110 ± 25 | 15 ± 2 | 87 ± 19 | 53 ± 10 | 84 ± 10 | 31 ± 14 | 44 ± 2 | 50.2 ± 3.3 | 60.6 ± 1.3 |
| plain MLP | 3 | 108 ± 37 | 35 ± 8 | 82 ± 71 | 66 ± 2 | 138 ± 115 | 53 ± 15 | 53 ± 11 | 62.4 ± 12.2 | 76.4 ± 20.2 |
| kNN | 3 | 164 ± 53 | 32 ± 1 | 87 ± 93 | 62 ± 12 | 48 ± 0 | 54 ± 1 | 34 ± 2 | 71.0 ± 12.1 | 68.9 ± 19.2 |

Paired masked − shared: TAL −19.4 (t = −10.2), HP −9.4 (t = −2.6), ethanol
+8.9, IBO +10.3, oxalic +5.2, succinic +1.3, lactic −2.1 (all |t| ≤ 1.1),
mean of 4 −4.7 (t = −1.0), mean of 7 −0.7. The seed spread of ethanol drops
from 37 to 10 points.

Per-channel out-of-fold error (RMS over a sample's points of the channel
error, relative to the RMS target numerator; configuration-balanced; shared →
masked, 5 seeds), leave-configurations-out / leave-one-product-out: feedstock
11.1 → 7.0 / 11.3 → 7.1, capex_upstream 3.7 → 2.0 / 3.3 → 1.7, coproducts
15.4 → 13.7 / 29.6 → 13.0, materials_other 11.7 → 9.5 / 16.5 → 13.4,
capex_fermentation 4.5 → 3.5 / 5.0 → 4.3, capex_wastewater 2.1 → 1.7 / 2.3 →
2.0, fixed 3.7 → 3.8 / 4.6 → 3.9, electricity 15.3 → 15.7 / 16.4 → 14.6,
natural gas 5.6 → 6.2 / 8.2 → 7.2, capex_utilities 8.2 → 7.6 / 7.1 → 7.7,
capex_downstream 2.8 → 2.8 / 5.1 → 4.6, neutralization 3.9 → 3.3 / 6.7 → 9.2.

Findings: the masks transfer what is physically product-independent
(feedstock, upstream, fermentation and coproduct channels), which helps the
families that share feedstocks and process areas with others (HP, TAL). The
errors of the products without close relatives come from g (e.g. LOPO:
lactic ĝ 0.99 vs 0.63, TAL 1.02 vs 0.71, oxalic 0.86 vs 0.997) and from
electricity, natural gas and the downstream channels, which legitimately
depend on the product, so masking cannot help them; the remaining bias of the
masked feedstock intensity (≈ 7 %) is the product-specific structure listed
in §3. The hybrid's residual absorbs these product-specific effects in-family
but extrapolates them to new products (ethanol LOPO 104–110 %). The masked
model (no residual) is therefore the pipeline default (`experiments.py`; CLI
`--no-masking`); the library default is unchanged.

### Development-set diagnostics, round 4 (2026-10-17; never on test data)

The v1 run (commit fd3f386, `results/v1/`) was the first use of the frozen
test split: OOD test PI-NN 14.5 % vs ridge 7.3 %, HGB 10.2 %; in-distribution
test 9.9 % vs random-CV 4.0 %; recovery g R² 0.04 on the OOD test; LOFO 66 % and
LOPO 75 % against kNN 36 / 60 %; productivity extrapolation 12.7 % vs HGB 6.8 %.
Round 4 looked for the causes on the development set only (every diagnostic
loads the samples and keeps `experiments.development_set(samples,
path=results/split.json)`; no decision used a test sample) and changed the
pipeline in four places (4a–4d). The agenda itself came from those v1 test
numbers, so the v2 run is an adaptive second use of the test split: its test
numbers are likely optimistic and do not independently confirm round 4 (§6).
Equally, the v2 development CV falls on the folds (LOCO, LOFO, LOPO,
productivity are deterministic; random CV seed 0) and the proxy draws 0–5 on
which round 4 was diagnosed and selected, so those figures are post-selection
estimates (4e gives fold-level statistics and a confirmation on fresh draws). 4e validates all changes together; 4f lists
what was adopted, what was not, and the open issues. Each change has a CLI
flag that restores v1 (`--v1-protocol`, `--v1-g-inputs`, `--no-r-structure`,
`--no-extrapolation-safeguards`). The v2 run uses the frozen test split a
second time.

#### 4a. Training protocol of the final ensembles

Diagnosis. v1 `final_dev` held out 15 % of the configurations (6 of 37) for
split-conformal calibration and every member early-stopped on another ~15 %
of the configurations (members stopped at 48–172 epochs). The in-distribution
test samples of the 6 calibration configurations came from configurations no
member had trained on, whereas the random CV used sample-level hold-outs:
hence ID 9.9 % vs random-CV 4.0 %. The baselines were fit on the same rows.

Dev proxy of the frozen test: on each of 6 draws one whole configuration of
each of HP, TAL, succinic and lactic is held out (OOD proxy, 78–102 samples)
plus 10 % hashed samples of the other configurations (ID proxy, 108–132); the
rest (~1100 samples) trains. v1 pipeline settings (selected hyperparameters,
water gate, 'physics' masking 32/3, downstream features), 3 members, the same
draws and seeds for every variant; MPSP MAPE %, mean ± sd over draws; `cfg` =
configuration-balanced; train s = CPU seconds for 3 members:

| variant | protocol | OOD | OOD cfg | ID | ID cfg | g R² OOD | best epoch | epochs run | train s |
|---|---|---|---|---|---|---|---|---|---|
| V0 | v1 final: calib 0.15 cfg-grouped + cfg-grouped ES | 11.95 ± 6.87 | 11.16 ± 5.89 | 13.94 ± 5.53 | 10.08 ± 4.15 | 0.58 ± 0.46 | 74 (min 6) | 114 | 50.5 |
| V1 | calib 0 + cfg-grouped ES | 9.19 ± 1.93 | 9.17 ± 2.87 | 8.15 ± 3.26 | 7.62 ± 3.09 | 0.67 ± 0.24 | 74 (min 6) | 114 | 59.4 |
| V2t | V1 + refit e* (truncated cosine) | 8.31 ± 2.07 | 7.94 ± 2.64 | 5.93 ± 1.81 | 6.34 ± 2.64 | 0.76 ± 0.16 | 74 | 114+74 | 102 |
| V2a | V1 + refit e* (annealed over e*) | 8.88 ± 2.40 | 8.48 ± 3.11 | 5.99 ± 1.37 | 6.23 ± 1.64 | 0.71 ± 0.20 | 74 | 114+74 | 102 |
| V2m | V1 + refit at the members' median e* | 7.86 ± 0.96 | 7.57 ± 1.43 | 4.78 ± 1.37 | 4.58 ± 1.33 | 0.88 ± 0.02 | 74 | 114+63 | 99 |
| **V3** | calib 0 + sample-level ES | **7.20 ± 1.20** | **6.69 ± 0.75** | **3.94 ± 0.22** | **3.87 ± 0.51** | 0.885 ± 0.02 | 172 (min 80) | 212 | 113 |
| V3r | V3 + refit e* (truncated) | 7.37 ± 1.33 | 6.80 ± 0.91 | 4.11 ± 0.33 | 3.96 ± 0.50 | 0.883 ± 0.02 | 172 | 212+172 | 210 |

Paired differences (pp, mean ± sd, t over 6 draws): V3 − V0 OOD −4.75 ± 6.37
(t −1.8), ID −9.99 ± 5.34 (t −4.6); V3 − V1 OOD −1.99 ± 2.79 (t −1.75), ID
−4.21 ± 3.28 (t −3.1); V3 − V2t ID −1.99 ± 1.84 (t −2.7); V3 − V2m OOD −0.65
(t −1.0), ID −0.84 (t −1.6); V2t − V1 OOD −0.88 ± 0.49 (t −4.4); V2a − V2t
+0.57 (t 1.3); V3r − V3 +0.16 ± 0.59 (t 0.7). Configuration-balanced OOD per
draw, V0 / V3: 21.2 / 6.7, 7.5 / 7.4, 15.7 / 7.2, 8.0 / 7.1, 6.9 / 6.4, 7.7 / 5.3.
Configuration-grouped early stopping chooses e* from about five
configurations: best epochs 6–226, and in 4 of 6 draws at least one member
stopped at ≤ 20 epochs, which poisons the ensemble mean. Part of the recovery
problem (lead's cause 2) is protocol-driven: g R² on the OOD proxy 0.58 →
0.885. V3 costs 2.2× V0 per fit (212 vs 114 epochs, 932 vs 778 samples per
member). Pass 1 of a refit run reproduces the plain run exactly (identical
state_dicts, tested). **Adopted: V3** (`experiments.PROTOCOL`) for every fit
of the pipeline; the refit code (`TrainConfig.refit`, 'truncate' / 'anneal',
`refit_epochs` 'member' / 'median') stays as a tested, default-off option.

Intervals without a hold-out (V3 procedure; 3 draws, α = 0.1). Out-of-fold
scores from an 8-fold configuration-grouped CV (stratified by biorefinery) of
the same procedure on each draw's training rows calibrate the proxy final
model; coverage is sample-simultaneous (a sample counts as covered only if all
its TY points are):

| key | q (sample wt) | cross-fit est. | per-cfg q10 | cov OOD proxy | cov OOD cfg-bal | cov ID proxy | q (group wt) | cov OOD (group wt) | cov ID (group wt) |
|---|---|---|---|---|---|---|---|---|---|
| a | 4.89 ± 0.78 | 0.929 | 0.85 | 0.979 ± 0.011 | 0.982 | 0.994 | 4.72 | 0.979 | 0.994 |
| b | 5.40 ± 1.80 | 0.932 | 0.92 | 0.846 ± 0.116 | 0.850 | 0.992 | 5.09 | 0.832 | 0.992 |
| c | 3.79 ± 0.34 | 0.926 | 0.76 | 0.919 ± 0.024 | 0.930 | 0.980 | 3.68 | 0.916 | 0.980 |
| d | 2.80 ± 0.64 | 0.940 | 1.00 | 0.873 ± 0.082 | 0.885 | 0.988 | 1.66 | 0.808 | 0.969 |
| g | 3.75 ± 0.90 | 0.918 | 0.68 | 0.879 ± 0.043 | 0.894 | 0.963 | 3.63 | 0.879 | 0.963 |
| MPSP | 8.57 ± 0.41 | 0.901 ± 0.028 | 0.67 | 0.783 ± 0.071 | 0.812 | 0.933 ± 0.045 | 9.16 | 0.808 ± 0.065 | 0.945 |

V0 calibrated on its own 6 calibration configurations covers the OOD proxy at
0.94 ± 0.10 for MPSP, but with q = 10.4 ± 6.8 (very wide); d 0.72, g 0.75.
Calibration from random sample-level 5-fold CV scores: MPSP q 7.2, ID 0.93 /
OOD 0.66 (g ID 0.92 / OOD 0.51). Coverage per configuration is uneven: the TAL
proxy configurations are under-covered in every draw (TAL_glucose 0.44,
TAL_corn 0.33) although TAL configurations in the out-of-fold scores are
covered at 0.85. The 8-fold grouped CV MAPE is 10–12 % (it includes the
isobutanol, oxalic and ethanol folds) against 6.7–7.4 % on the proxy (four
large families only). **Adopted:** the final intervals come from the
leave-one-configuration-out out-of-fold scores of the development CV (same
procedure, same ensemble size; `experiments.stage_final`), sample weighting;
`transfer_conformal`, which used the test configurations, is removed.

Cost (per-fit factors from the table: V3 / V0 2.23, V3 / V1 1.90, refit
+0.85×): v1 CV 4.06 CPU-h → about 8.7 CPU-h with V3 everywhere and every
ablation on every LOCO fold, 5.3 CPU-h with the ablations on 12 stratified
LOCO folds (`--ablation-loco-folds`, default); HPO 0.58 → ~1.1 CPU-h; final +
release ~0.2 CPU-h.

#### 4b. Inputs of the recovery g


Data (development set: 1306 samples, 37 configurations, 29 284 points;
scripts in the round-4 scratch directory `g_inputs/`). What drives g:

* Between versus within configurations: 96.9 % of the variance of ln g lies
  between configurations (RMS 0.20), 3.1 % within (RMS 0.036). The
  fermentation product alone explains R² = 0.94 (configuration-balanced
  0.91); product + separation variant + neutralization + final product
  0.999. The feedstock alone explains 0.18 / 0.08, and that only through
  confounding (lactic only on cornstover, oxalic only on sugarcane).
* Feedstock within a product specification (configuration medians): HP with
  hexanol extraction 0.675–0.695 (2.8 %, cornstover highest), HP neutral +
  hexanol 2.3 %, HP without hexanol and HP neutral ≤ 0.4 %, TAL / TAL_SA ≤
  1.0 %, succinic ≤ 0.1 %; ethanol corn 0.9937 vs 0.9995 (cornstover,
  sugarcane), which the separation inputs already distinguish (the dry-grind
  DDGS train). Neutralization changes g within HP (+4 %: 0.600 → 0.625,
  hexanol 0.680 → 0.707) and lactic (−4 %); it is a g input already.
  Simultaneous saccharification and co-fermentation (lactic lowpH SSCF vs
  SHF: 0.627 vs 0.659, −5 %; residual solids go with the broth to
  separation) is represented by no feature: `fermentation_mode` is 'batch'
  for both and the `saccharification` field exists only in the
  configurations table of the lactic module.
* Within configurations no observed specification moves g: pooled
  within-configuration correlation of ln g (random samples, configuration-
  demeaned) with log capacity −0.01, operating days −0.03, cell-mass yield
  −0.03, base price −0.05, acid price −0.02, log r −0.06 (slopes ≤ 0.25 %
  per sd). The within-configuration spread comes from unobserved `bp_*`
  parameters: TAL sd(ln g) = 0.11 (`bp_separation_tal_decarboxylation_conversion`,
  r = −0.76 … −0.96; a 10 % MPSP noise floor for TAL from g alone),
  succinic 0.005 (crystallization pressure-filter recovery, r ≈ +0.75),
  lactic 0.01–0.03, oxalic 0.003 (lime loading, r = 1.00), HP ≤ 0.003,
  ethanol ≤ 0.0005. Scale and productivity therefore have no place among the
  g inputs, and the reagent prices (in the 'neutralization' group until now)
  neither.
* Titer: the per-point recovery depends on the titer within a sample
  (d ln recovery / d ln t: lactic 0.27–0.32, i.e. 26–33 % over the domain;
  TAL 0.06–0.08; HP with hexanol 0.04; HP 0.02; IBO 0.03; ethanol, succinic,
  oxalic ≈ 0), but g is the median over a TY domain that is fixed per
  product family (HP 10–200 g/L, TAL 10–100, lactic 40–200, succinic 10–120,
  oxalic 5–100, IBO 5–80, ethanol 10–58 / 140 / 150), so the median has no
  within-family dependence on the domain, and across products none either
  (lactic: highest domain titer, lowest g; IBO: lowest, g 0.985). A domain
  feature is a product-identity proxy on these data. The TY-domain bounds
  are not features (only the broth product's `phys_log10_w_final` uses a
  reference titer).
* Consistency of the target: a′ = g·A to machine precision (relative
  1e−13), the median recovery over all ok points equals the median over the
  fitted-domain points (lactic −1.0 / −1.3 %, all others 0), g ≤ g_max
  except oxalic (24 samples, g/g_max ≤ 1.0029, absorbed by `g_margin`). The
  equation is invariant to the g convention (the primed coefficients carry
  g); fitting MPSP × per-point recovery instead and dividing by the median
  g would double the fit error (median 1.35 → 2.45 %, lactic 7–9 %), so the
  titer dependence of the recovery is correctly left to the coefficients.
* Where the v1 g errors are (v1 development CV, out-of-fold, PI-NN):
  configuration-mean |ln(ĝ/g)| random 2.1 %, leave-configurations-out
  4.9 % (MPSP 14.4 %; 11.9 % with the true g), leave-one-product-out 22.3 %
  (MPSP 47.8 %; 34.1 % with the true g), leave-one-feedstock-out 6.8 %, of
  which nearly all is the product-unique feedstocks (cornstover held out:
  every lactic configuration is a new product, +45 %; sugarcane: oxalic,
  −28 %); otherwise TAL 7–12 % (its within-configuration noise), HP ≤ 4.6 %,
  succinic and ethanol ≤ 1.4 %. v1's
  HP_sugarcane +8 % (LOCO) is not an input problem: its g inputs are those
  of HP_glucose / HP_cornstover (0.600 / 0.603); with the protocol below it
  is +1.3 %.

Candidate input groups of the recovery sub-network, paired development CV
(seeds 0–2; one network per fit, no conformal hold-out, early-stopping
hold-out of 15 % grouped as `experiments.INNER_GROUP`, selected
hyperparameters, downstream features, water gate, masked width 32 / depth
3; leave-configurations-out = 5 folds of whole configurations stratified by
biorefinery (`stratified_group_kfold`, seed 0, every configuration held out
once), leave-one-product-out 7 folds, leave-one-feedstock-out 4 folds).
Reference 'physics' (identical, bit for bit, before and after the
reagent-price split); MPSP MAPE = mean over held-out configurations (LOCO) or
held-out groups (LOPO, LOFO); g error = configuration mean of |ln(ĝ/g)|;
mean ± sd over seeds, paired Δ (t over 3 seeds):

| g inputs | LOCO MPSP | LOCO g MAE | LOCO g R² | LOCO g err | LOPO MPSP | LOPO g err | LOFO MPSP | LOFO g err | LOPO degenerate stops |
|---|---|---|---|---|---|---|---|---|---|
| 'physics' | 12.8 ± 3.6 | 0.057 | 0.62 | 3.9 % | 75.5 ± 15.6 | 24.9 % | 75.5 ± 17.4 | 7.9 % | 19 % |
| + feedstock | 11.8 (−1.0, t −2.7) | 0.041 (t −7.1) | 0.78 (t 2.9) | 3.9 % (t 0.8) | 82.8 (+7.2, t 0.7) | 23.7 % (t −1.2) | 76.7 (+1.2, t 0.3) | 8.8 % (t 2.1) | 38 % |
| + feedstock + fermentation (seed 0 only; reference seed 0 in brackets) | 12.4 (11.7) | 0.050 (0.046) | 0.73 (0.76) | 3.7 % (3.2) | 106 (61) | 28.2 % (27.4) | 64.7 (70.8) | 12.5 % (9.9) | 29 % (0 %) |
| + titer_domain (`domain_features`) | 12.3 (−0.5, t −1.3) | 0.054 (t −1.3) | 0.67 (t 1.4) | 3.9 % (t 0.4) | 87.3 (+11.8, t 1.2) | 25.0 % (t 0.4) | 77.5 (+2.0, t 0.5) | 7.9 % (t −0.1) | 33 % |
| − reagent_prices ('physics_v2') | 12.5 (−0.3, t −1.1) | 0.065 (t 1.5) | 0.47 (t −1.4) | 4.1 % (t 1.5) | 65.1 (−10.5, t −1.2) | 24.5 % (t −2.9) | 71.1 (−4.4, t −1.5) | 8.0 % (t 0.5) | 29 % |
| 'physics', w_g = 5 | 11.5 (−1.4, t −1.8) | 0.044 (t −1.8) | 0.76 (t 1.7) | 3.4 % (t −1.8) | 81.2 (+5.7, t 2.0) | 25.4 % (t 0.6) | 75.2 (−0.3, t −0.4) | 7.1 % (t −1.4) | 33 % |

* The feedstock's apparent LOCO gain is the broth product alone (oxalic,
  only on sugarcane: MPSP 80 → 53 %, g error 21 → 9 %); it is an
  extrapolation shift for an unseen water route, not recovery physics. The
  families with feedstock variants get worse: HP g error 0.9 → 1.7 % (MPSP
  4.0 → 4.2 %, t 3.4), TAL 10.5 → 11.0 %, ethanol 30.5 → 32.3 %; under LOFO
  the held-out sugarcane g error grows from −4.6 to −11.5 % and under LOPO
  twice as many fits stop with their best epoch at 0 (38 vs 19 %). The network does not
  resolve the ≤ 3 % within-family feedstock effect; fermentation, scale and
  the titer domain have no within-family variation that the median g
  follows (above). None of them is adopted.
* Without the reagent prices g is unchanged in-family (HP, TAL, succinic,
  ethanol within 0.1 point) and LOPO / LOFO are better within the seed
  noise (configuration-balanced LOFO −2.5 points, t −2.4); the oxalic LOCO
  fold is worse (MPSP 80 → 100 %, t 1.6; the same unseen-route
  extrapolation). Adopted as a hard constraint (∂g/∂price = 0, as observed:
  within-configuration r = −0.05 / −0.02): preset `'physics_v2'`, pipeline
  default; `'physics'` is kept so that saved checkpoints load unchanged.
* w_g = 5 helps LOCO a little (IBO g error 23 → 11 %) and hurts LOPO (t
  2.0; configuration-balanced t 3.8): w_g stays 1. The g target and its
  bound need no change (consistency checks above).
* Early stopping under LOPO: the inner hold-out is one whole training
  product, and on it the criterion `data` is dominated by the extrapolated
  coefficient and recovery terms (TAL as hold-out: tot 2.2, g 2.3, channels
  1.5, MPSP 0.2 at epoch 0). When the recovery of the hold-out product
  drifts away during training (e.g. with feedstock inputs: g term 2.4 →
  3.6–4.0) the best epoch is 0 and the fit is essentially untrained (19 %
  of the reference LOPO fits, 29–38 % for the variants). LOPO
  comparisons of g variants are therefore dominated by this selection noise
  (ethanol and IBO MPSP swing by 30–100 points between variants and seeds).
* Missing information, not features: the TAL recovery noise floor (10 %,
  `bp_separation_tal_decarboxylation_conversion`) and the saccharification
  mode (SHF / SSCF, −5 % for lactic) — the latter would belong in the regime
  specification (`saccharification`, recorded only by the lactic module);
  with one SSCF development configuration it cannot be validated here. The
  opt-in titer-domain feature (`FeatureSpec(domain_features=True)`, group
  'titer_domain', `phys_log10_t_ref` = log10 sqrt(t_min·t_max)) is kept for
  reproducibility only.
#### 4c. Hard productivity structure

Data analysis and paired CV evidence: §3 ("Hard productivity structure").
Only the fermentation capital depends materially on r, as 1/r (β = 1.0, R²
0.9999); with the structure ∂MPSP/∂r ≤ 0 holds by construction (v1: 0.8–8 %
of held-out samples violated it) and the productivity-extrapolation MAPE falls
from 5.16 to 4.29 % (t −4.5; low-r side 6.34 → 4.39 %, t −12.5) with
sample-level hold-outs, LOCO and random unchanged. **Adopted** in the
pipeline (`R_STRUCTURE = 'physics'`, β fixed at 1).

#### 4d. Catastrophic extrapolation (new feedstocks and products)

Development set 1306 samples / 37 configurations (`experiments.development_set`);
pipeline settings (downstream features, water gate, masking 'physics' with
masked width/depth = width/2, depth − 1), HPO-selected hyperparameters
(width 64, depth 4, w_ch 1, w_phys 0.03, ...), `calib_frac` 0.15 and the v1
`INNER_GROUP` unless stated; one member per fit (two for the reproduction);
the same folds, inner splits and seeds for every variant (paired);
leave-configurations-out as the 5-fold biorefinery-stratified split of rounds
1–3 (the pipeline's 37-fold LOCO was too expensive to repeat per variant).

*Reproduction* (2 members, seeds 0–2): leave-one-feedstock-out corn 145 ± 41 %
(v1 150), leave-one-product-out ethanol 154 ± 65 % (v1 177), isobutanol
210 ± 81 % (v1 148). Per-channel error of the held-out numerator (% of the
true numerator, mean over the configuration's samples and seeds):

| held out | configuration | MAPE | ĝ / g | 1/η̂ / 1/η | feedstock | electricity | nat. gas | neutral. | materials | capex ferm. | capex util. | coproducts |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| corn | HP / TAL / succinic (7 cfgs) | 23–47 | ≈ g (TAL +4–7 %) | 1.50 / 0.93–1.01 | +26…+40 | −6…−12 | ≈ 0 | ≈ 0 | −2…+11 | +1…+7 | +3…+7 | ≈ 0 |
| corn | ethanol_corn | 147 | 0.98 / 0.99 | 1.50 / 1.00 | +77 | −50 | −10 | +4 | +13 | +10 | +41 | +61 |
| corn | IBO_EtOH_corn | 287 | 0.92 / 0.98 | 1.50 / 0.98 | +136 | −93 | −43 | +3 | +34 | +12 | +17 | +181 |
| product | ethanol (corn / stover / cane) | 202 / 183 / 87 | 0.81 / 1.00 | 1.06 / 1.18 / 1.16 vs 1.00 / 1.18 / 1.80 | +12 / +6 / −20 | −29 / +18 / +14 | −3 / +17 / +9 | +18 / +18 / +13 | +43 / +25 / +32 | +21 / +21 / +17 | +31 / +11 / −14 | +56 / 0 / −2 |
| product | IBO_EtOH_corn | 210 | 0.92 / 0.98 | 1.09 / 0.98 | +35 | −49 | −34 | +13 | +17 | +18 | +3 | +164 |

Causes. (i) *Leave-one-feedstock-out corn*: corn's minor mass fraction
`fs_w_other` (0.10–0.13) lies +9.5 training sd beyond the other feedstocks'
(0–0.015, sd 0.006; standardized +18, cut at the FeatureSpec's 10-sd clip);
`fs_w_starch`/`fs_w_fiber` are constant in training and dropped, the
`feedstock_type` level is unseen. The feedstock sub-network saturates and
returns one constant 1/η = 1.50 for every corn sample (truth 0.93–1.01): +26
to +40 % of the numerator for HP/TAL/succinic and +77 / +136 % for ethanol /
IBO, whose costs are feedstock-dominated. It is not an unbounded-`z`
blow-up: 1.50 lies inside the training range of 1/η (0.52–1.80), so a
data-derived bound on `z` would not bind. The rest of the corn error is
corn-specific process structure no other feedstock shows: the DDGS / ethanol
coproduct credits (−60 / −173 % of the numerator, predicted ≈ 0) and corn
ethanol's purchased electricity and missing boiler (electricity −50 %,
capital of utilities +41 %). (ii) *Leave-one-product-out*: with the inner
calibration and early-stopping hold-outs grouped by product, every member
trains on 4 of the 6 remaining products and early-stops on one unseen product
whose validation loss is lowest in the first ~25 epochs (stops at 41–66 of
300 with patience 40): the members are nearly untrained, e.g. the
product-blind feedstock intensity of a held-out sugarcane product is 1/η =
1.0–1.2 instead of 1.80. On top come g (ethanol ĝ 0.81 vs 1.00), the
neutralization cost a held-out alcohol inherits from TAL / HP (+13–18 %),
materials, fermentation capital and IBO's ethanol coproduct credit (no other
product has it). The HPO-selected hyperparameters are not the cause (below).
Channel outputs stay inside the training range: 0–0.5 % of the predicted
(sample, channel, coefficient) entries of these folds lie outside the
training range ± 5 %, so bounding channel outputs would not bind either.

*The `fs_eta_le_1` check* (`b'_fs·w_sugar ≥ 1`, disabled by
`apply_constraint_report` because it fails for 33 % of the development
samples): lactic (3 configurations; 0.52–0.85, but the effective
1/η(y, t) = y·w_sugar·I_fs(y, t) ≥ 1.02 at every point: its feedstock sits
partly in d′ = 44–56, a basis split, not a physical violation); TAL corn,
glucose, SA-corn and SA-glucose (0.86–1.06; the effective 1/η falls to
0.68–0.72 at y = 0.267, t = 10 g/L and is < 1 at 76–80 % of their points:
the data imply less sugar than the product needs, the substrate-basis defect
D8, c′ = −17 … −22); HP on corn and glucose (8 configurations, 0.93–0.99,
effective ≥ 0.97, d′ 0.8–2.2: ≤ 3 %); IBO 0.98 and ethanol corn 0.999
(rounding). 13 % of all development points have an effective 1/η < 1.

*Safeguards* (§3, *Extrapolation safeguards*; `pnc` = sugar-availability
prior + neutralization gate + feedstock-input clipping), point MAPE %, mean ±
sd over seeds; paired difference to the baseline (the pipeline model without them) with t.

Leave-one-feedstock-out, leave-configurations-out and random: seeds 0–4;
leave-one-product-out: seeds 0–2 with the v1 product-grouped inner hold-outs
(`sel`, `pnc`) and seeds 0–4 with configuration-grouped ones (`pgrp`,
`pncg`):

| scheme (folds) | held out | baseline | `pnc` | `pnc` − baseline (t) |
|---|---|---|---|---|
| leave-one-feedstock-out | fold mean | 63.4 ± 14.5 | 34.4 ± 8.7 | −28.9 ± 12.2 (−5.3) |
| | corn | 127.5 ± 48.1 | 42.6 ± 5.1 | −84.8 (−4.1) |
| | corn stover | 40.0 ± 14.7 | 36.6 ± 6.0 | −3.3 (−0.5) |
| | glucose | 35.9 ± 13.1 | 7.7 ± 2.7 | −28.2 (−4.7) |
| | sugarcane | 50.2 ± 27.4 | 50.8 ± 29.4 | +0.7 (+0.1) |
| | configuration geo-mean | 28.6 ± 3.4 | 17.0 ± 0.4 | −11.7 ± 3.4 (−7.8) |
| leave-one-product-out, v1 inner grouping | fold mean | 81.9 ± 8.8 | 78.2 ± 8.2 | −3.7 ± 0.7 (−9.6) |
| | ethanol | 161.4 ± 84.1 | 148.0 ± 64.2 | −13.5 (−1.2) |
| | HP | 21.5 ± 0.5 | 20.1 ± 1.5 | −1.4 (−1.7) |
| | isobutanol | 193.3 ± 132.2 | 170.8 ± 115.0 | −22.5 (−2.2) |
| | lactic | 58.2 ± 5.3 | 56.8 ± 6.7 | −1.4 (−1.7) |
| | oxalic | 49.0 ± 33.8 | 61.5 ± 34.6 | +12.5 (+2.1) |
| | succinic | 36.9 ± 14.8 | 37.1 ± 13.3 | +0.3 (+0.2) |
| | TAL | 53.2 ± 4.8 | 53.1 ± 5.8 | −0.1 (−0.2) |
| leave-one-product-out, configuration-grouped | fold mean | 48.3 ± 7.7 | 47.1 ± 9.4 | −1.2 ± 3.2 (−0.8) |
| | ethanol | 65.2 ± 49.5 | 62.3 ± 41.6 | −2.9 (−0.7) |
| | HP | 20.1 ± 4.6 | 21.0 ± 8.2 | +0.9 (+0.4) |
| | isobutanol | 59.3 ± 8.6 | 53.9 ± 2.8 | −5.4 (−1.4) |
| | lactic | 51.1 ± 8.9 | 48.3 ± 8.6 | −2.9 (−2.2) |
| | oxalic | 51.2 ± 23.0 | 53.0 ± 27.3 | +1.7 (+0.8) |
| | succinic | 42.5 ± 15.8 | 43.2 ± 16.2 | +0.7 (+0.3) |
| | TAL | 48.9 ± 2.6 | 48.2 ± 3.7 | −0.7 (−1.0) |
| leave-configurations-out (5) | fold mean | 32.7 ± 4.7 | 29.6 ± 4.7 | −3.1 ± 2.2 (−3.1) |
| | configuration mean | 16.3 ± 1.8 | 15.7 ± 2.0 | −0.5 ± 1.0 (−1.2) |
| | configuration geo-mean | 9.2 ± 1.1 | 9.3 ± 1.1 | +0.0 ± 0.7 (+0.1) |
| random (5) | fold mean | 4.30 ± 0.09 | 4.24 ± 0.16 | −0.06 ± 0.15 (−0.9) |
| | configuration geo-mean | 3.56 ± 0.14 | 3.45 ± 0.22 | −0.12 ± 0.20 (−1.3) |

Mean |error| of the held-out feedstock channel (% of the true numerator,
configuration-balanced, seeds 0–2), baseline / clipping only / neutralization
gate only / `pnc`: leave-one-feedstock-out 23.2 / 12.2 / 23.4 / 3.5 (corn 46.8
/ 21.1 / 46.6 / 7.5, sugarcane 16.3 / 19.8 / 16.3 / 4.1, glucose 25.4 / 3.9 /
25.6 / 1.5), |log(1/η̂ ÷ 1/η)| 33 / 29 / 33 / 8 %; leave-one-product-out
12.1 (v1 grouping) and 6.5 (configuration-grouped) → 5.7 and 2.8 with `pnc`.
Neutralization channel, leave-one-product-out (configuration-grouped): 6.9 →
4.5 (gate only) and 4.4 (`pnc`). Clipping alone reaches the same
leave-one-feedstock-out MAPE (38.8 ± 9.6 %, −28.9 ± 9.5, t = −5.3; corn 51.8
± 25.1, sugarcane 53.5) but for partly wrong reasons: corn 1/η̂ 1.30 (truth
0.98) and sugarcane 0.97 (truth 1.76–1.80, −20 to −29 % of the numerator,
offsetting overpredicted electricity and materials); only with the prior is
the feedstock channel right (corn spread over seeds 3.8 vs 25.1 points). The
gate alone: leave-one-product-out (configuration-grouped) −1.2 ± 0.9 (t =
−2.3), ethanol −8.3; leave-one-feedstock-out −2.8 ± 2.7 (t = −1.8).

With two members (seeds 0–2): corn 145 → 60 % (−85, t = −3.1), ethanol 154
→ 145 (−10, t = −0.7), isobutanol 210 → 191 (−19, t = −2.1). The corrected
feedstock channel is visible per configuration (leave-one-feedstock-out
corn, 3 seeds): 1/η̂ 1.29 → 1.03, feedstock error of HP/TAL/succinic corn
+25…+38 → +4…+8 % of the numerator, their MAPE 25–49 → 6–15 %; corn stover
1.30 → 1.11 (truth 1.04–1.19), glucose 0.98, sugarcane 1.63 (truth
1.71–1.80); the neutralization error of held-out alcohols +13–18 → 0 %.
Oxalic (leave-one-product-out) gets worse because the baseline was right
for the wrong reason: its nearly untrained feedstock intensity (1/η̂ 1.00
vs 1.79, −20 % of the numerator) offset the overpredicted electricity,
materials and recovery (+24, +22 %, ĝ 0.87); with the prior the feedstock
error is −1 % and the offset is gone. Remaining failures are product- or
feedstock-specific process structure no training configuration shows
(coproduct credits of IBO / corn ethanol, corn ethanol's missing boiler, g
of new products), which no safeguard can supply.

Hyperparameters (default vs HPO-selected, paired, seeds 0–2, one member):
leave-one-product-out fold mean 81.3 ± 12.8 vs 81.9 ± 8.8 % (−0.7 ± 7.3,
t = −0.2; ethanol 141 vs 161, IBO 189 vs 193); with two members ethanol
+0.3 (t 0.0) and IBO +3.4 (t 0.8); leave-one-feedstock-out the defaults are
worse (+13.6, t = 6.8; corn 190 vs 129). The selected hyperparameters do
not explain the gap to round 3 (masked LOPO ethanol 83 ± 10 %): the inner
grouping does. With configuration-grouped inner hold-outs (`pgrp`, as in
round 3) the leave-one-product-out fold mean falls from 81.9 ± 8.8 to 49.7 ±
10.0 % (−32.2 ± 5.4, t = −10.3; ethanol 161 → 77, IBO 193 → 62, oxalic 49 →
35, succinic 37 → 50; seeds 0–2). (This diagnosis proposed
`INNER_GROUP['product'] = 'config_id'`; the v2 protocol, round 4a, goes
further: sample-level early stopping in every scheme.) Under
that grouping the safeguards change leave-one-product-out by −2.9 ± 1.9
(t = −2.7, seeds 0–2) and −1.2 ± 3.2 (t = −0.8, seeds 0–4): neutral to
slightly better in MAPE, with the feedstock and neutralization channels
right (above). Seeds 3–4 of `sel` / `pnc` ran after that change, so their
leave-one-product-out folds count as `pgrp` / `pncg`.

#### 4e. Combined validation of v2 (all changes, corrected protocol)

Setup: development set only (1306 samples, 37 configurations); the
pipeline's folds (`experiments.scheme_splits` with `results/plan.json`, seed 0:
random 5-fold, all 37 LOCO folds, 4 LOFO, 7 LOPO, 1 productivity fold) and the
dev proxy of the frozen test (4a, 6 draws). PI-NN with the **v1 settings**
(v1 protocol: CV hold-outs grouped by `INNER_GROUP_V1`, calibration 0.15;
proxy: the v1 final protocol; masking 'physics'; no r-structure; no
safeguards) vs **v2** (all adopted changes: `PROTOCOL`, 'physics_v2',
r-structure, safeguards), both with the v1-selected hyperparameters, water
gate and downstream features, 3 members, paired over 3 training seeds (seed
1000 s + fold) on identical folds; plain MLP (v2 protocol), ridge, HGB and kNN
(fit on the same rows as the v2 members) once, for reference. Scripts and raw
results: round-4 scratch directory `integration/` (`validate.py`,
`analyze.py`, `tables.md`). MPSP MAPE %, fold mean (pipeline headline),
mean ± sd over seeds; diff = v2 − v1, paired, t over seeds:

| scheme | folds | v1 | v2 | v2 − v1 (t) | plain | ridge | HGB | kNN |
|---|---|---|---|---|---|---|---|---|
| random | 5 | 4.10 ± 0.02 | **3.91 ± 0.03** | −0.20 ± 0.04 (−9.1) | 5.77 | 8.15 | 4.88 | 8.02 |
| leave-one-configuration-out | 37 | 16.06 ± 0.92 | **9.71 ± 0.25** | −6.35 ± 1.01 (−10.8) | 13.76 | 19.46 | 21.51 | 16.34 |
| leave-one-feedstock-out | 4 | 64.13 ± 1.43 | 34.89 ± 1.33 | −29.25 ± 2.10 (−24.1) | 71.57 | 141.80 | 73.69 | 34.63 |
| … configuration-balanced | 4 | 38.98 ± 1.36 | **21.96 ± 0.42** | −17.03 ± 1.67 (−17.6) | 40.64 | 95.54 | 32.92 | 23.77 |
| leave-one-product-out | 7 | 86.12 ± 9.47 | **50.08 ± 3.51** | −36.03 ± 8.36 (−7.5) | 51.29 | 88.25 | 107.69 | 52.57 |
| productivity | 1 | 18.30 ± 5.76 | **4.30 ± 0.07** | −13.99 ± 5.71 (−4.2) | 6.86 | 8.92 | 7.26 | 10.64 |
| LOCO, HP / TAL / succinic / lactic, family-balanced (matched OOD estimate) | 32 | 13.25 ± 1.27 | **6.23 ± 0.02** | −7.02 ± 1.28 (−9.5) | 9.79 | 11.65 | 8.85 | 13.93 |

Fold-level statistics (seeds averaged per fold, then paired over folds; the
seed t of the table measures training-seed noise on fixed folds only, not
whether a new configuration, feedstock or product benefits; script
`review_fix/foldstats.py` in the round-4 scratch directory), v2 − v1 in MAPE
points:

| scheme | units | mean | median | better / worse | t over units | Wilcoxon p | largest unit's share of Σ\|Δ\| |
|---|---|---|---|---|---|---|---|
| random | 5 folds | −0.20 | −0.19 | 5 / 0 | −5.7 | 0.06 | 0.28 |
| leave-one-configuration-out | 37 folds | −6.35 | −3.70 | 35 / 2 | −4.4 | 6·10⁻⁸ | 0.13 |
| … four OOD families | 32 configurations | −5.20 | −3.26 | 31 / 1 | −4.9 | 3·10⁻⁹ | 0.16 |
| leave-one-feedstock-out | 4 folds | −29.25 | −14.89 | 3 / 1 | −1.3 | (sign p 0.63) | 0.72 (corn) |
| … configuration-balanced | 4 folds | −17.03 | −13.31 | 3 / 1 | −1.2 | (sign p 0.63) | 0.58 (corn) |
| leave-one-product-out | 7 folds | −36.03 | −5.58 | 6 / 1 | −1.5 | 0.16 | 0.48 (isobutanol) |
| productivity (1 fold) | 36 held-out configurations | −10.43 | −8.04 | 36 / 0 | −9.2 | 3·10⁻¹¹ | 0.08 |
| proxy OOD | 6 draws | −4.73 | −1.50 | 6 / 0 | −1.9 | 0.03 | 0.53 (draw 0) |
| proxy ID | 6 draws | −10.50 | −10.72 | 6 / 0 | −4.3 | 0.03 | 0.31 |

The LOCO and productivity gains hold at the configuration level. The LOFO
(4 folds) and LOPO (7 folds) gains are **not statistically established at
the fold level**: one fold dominates each (corn, isobutanol), the LOPO
median is only −5.6 points, and cornstover (LOFO, +8.2) and succinic (LOPO,
+23) got worse. The proxy draws overlap (24 picks, 16 distinct
configurations; TAL_glucose and lactic_lowpH_SSCF_cornstover in 3 draws
each), so a t over draws overstates their independence. The same holds for
the single changes: 4d's LOFO −28.9 (seed t −5.3) is 3 of 4 feedstocks
better, median −15.8, t over folds −1.5 (corn 72 % of Σ|Δ|); 4b's LOPO
−10.5 is 3 better / 4 worse over the products, median +0.1; 4c's
productivity −0.9 is 28 of 36 configurations better (t −2.9).

Per held-out group (v1 → v2; plain / ridge / HGB / kNN): LOFO corn 145.6 →
50.2 (t −9.0; 165 / 351 / 83 / 52), glucose 29.9 → 5.1 (t −4.4), sugarcane
46.2 → 41.2 (t −0.6), cornstover 34.9 → 43.1 (t +1.0); LOPO ethanol 157.7 →
44.4 (t −2.4; seed sd 83 → 1.6), isobutanol 226.9 → 85.1 (t −3.5), HP 18.2 →
12.7 (t −3.6), lactic 56.6 → 53.5 (t −2.6), oxalic 65.3 → 59.3 (t −0.2), TAL
49.6 → 44.3 (t −1.9), **succinic 28.4 → 51.4 (t +5.0)**. LOCO per family
(v1 → v2): HP 5.3 → 2.9, TAL 17.6 → 8.7, succinic 7.2 → 2.9, lactic 22.9 →
10.4, ethanol 38.8 → 21.8, oxalic 68.3 → 38.8, **isobutanol 87.6 → 99.5**
(one configuration). The succinic LOPO regression comes with the protocol
(4d found succinic 37 → 50 % when the product-grouped early stopping was
replaced by configuration-grouped early stopping): with product-grouped
hold-outs the v1 members stopped within ~8 epochs (87 % of the LOPO members
had best epoch ≤ 20) and the nearly untrained members happened to be closer
for succinic; the v2 LOPO members' best epochs have median 123 (min 54).

Out-of-fold recovery g (pooled R², configuration-mean |ln ĝ/g|), v1 → v2:
LOCO 0.66 → 0.80 (6.4 → 3.3 %), productivity 0.78 → 0.97 (7.9 → 2.4 %), LOPO
−0.56 → −0.24 (23.0 → 20.9 %), LOFO 0.50 → 0.39 (7.5 → 6.7 %), random 0.973 both.
Finite-difference ∂MPSP/∂r > 0 rate of the held-out samples: v2 0 in every
scheme (v1 0.2–5.9 %, plain MLP 11–36 %). Best epochs, median (min; fraction
≤ 20) over the members: v1 LOCO 67 (1; 0.27), LOPO 8 (1; 0.87), productivity
23 (11; 0.44), LOFO 57 (8; 0.17), random 142 (83; 0); v2 121–139 (min 52; 0)
in every scheme. Seconds per 3-member fit, v1 → v2: LOCO 65 → 126, LOPO 26 →
100, LOFO 43 → 88, productivity 39 → 79, random 84 → 87.

Cross-fitted coverage of the LOCO out-of-fold scores (nominal 0.90; leave one
configuration out; seeds 0 / 1 / 2): v2 MPSP mean 0.94 / 0.94 / 0.89, 10 %
quantile over configurations 0.89 / 0.87 / 0.67, worst 0.0 / 0.21 / 0.19; g
mean 0.94 / 0.86 / 0.90; a′ 0.93. v1: MPSP 0.94 / 0.95 / 0.94 (worst 0.0–0.03).

Dev proxy of the frozen test (6 draws, MAPE %, mean ± sd; paired over draws):

| variant | OOD | OOD cfg | ID | ID cfg | g R² OOD |
|---|---|---|---|---|---|
| PI-NN v1 settings (v1 final protocol) | 11.81 ± 6.61 | 11.02 ± 5.68 | 14.28 ± 6.00 | 10.23 ± 4.15 | 0.58 ± 0.46 |
| **PI-NN v2** | **7.08 ± 0.96** | **6.56 ± 0.49** | **3.79 ± 0.26** | **3.79 ± 0.46** | 0.88 ± 0.02 |
| plain MLP (v2 protocol) | 10.30 ± 2.10 | 9.94 ± 1.87 | 6.13 ± 0.70 | 5.85 ± 0.59 | 0.86 ± 0.05 |
| ridge | 11.44 ± 2.18 | 12.14 ± 3.04 | 7.73 ± 0.98 | 7.66 ± 1.02 | 0.88 ± 0.01 |
| HGB | 9.22 ± 1.36 | 9.15 ± 1.94 | 4.79 ± 0.51 | 4.91 ± 0.70 | 0.86 ± 0.05 |
| kNN | 17.29 ± 2.32 | 15.71 ± 1.86 | 7.93 ± 0.70 | 8.41 ± 0.57 | 0.84 ± 0.02 |

v2 − v1: OOD −4.73 ± 6.16 (t −1.9), ID −10.50 ± 5.99 (t −4.3); v2 − HGB: OOD
−2.14 ± 1.18 (t −4.5), ID −1.01 ± 0.29 (t −8.7); v2 − ridge: OOD −4.36 (t
−7.5), ID −3.95 (t −12.7); v2 − plain: OOD −3.22 (t −5.1), ID −2.34 (t −8.5).
Configuration-balanced OOD per draw, v1 / v2: 20.7 / 6.4, 7.7 / 6.8, 15.2 /
7.2, 7.7 / 6.6, 6.8 / 6.7, 8.0 / 5.7: v2 is never worse and its draw-to-draw
spread is 12× smaller. On the proxy (unlike the v1 test) v2 beats every
baseline on both OOD and ID. Draws 0–5 also chose the protocol (4a), so this
table is not out-of-sample for it.

Confirmation on draws and folds that played no role in any round-4 decision
(after the review; script `review_fix/confirm.py`, tables
`review_fix/confirm_tables.txt`; same settings as above, development set
only). Dev proxy draws 6–11 (one training seed per draw, 100 + draw) and the
random 5-fold CV with fold seed 1 (the pipeline uses seed 0), v1 settings vs
v2 paired over 3 training seeds (1000 s + 500 + fold). No fresh LOCO, LOFO,
LOPO or productivity folds exist (they are deterministic) and every
development configuration was in the round-4 LOCO, so the fresh draws are
new combinations of held-out configurations, new ID samples and new seeds,
not unseen configurations (24 picks, 13 distinct;
lactic_lowpH_SHF_cornstover in 5 draws). MAPE %, mean ± sd over draws:

| variant | OOD | OOD cfg | ID | ID cfg | g R² OOD |
|---|---|---|---|---|---|
| PI-NN v1 settings (v1 final protocol) | 13.29 ± 6.44 | 13.02 ± 5.94 | 13.86 ± 6.47 | 11.15 ± 6.53 | 0.46 ± 0.63 |
| **PI-NN v2** | **6.88 ± 1.41** | **6.40 ± 0.87** | **3.99 ± 0.33** | **3.94 ± 0.55** | 0.88 ± 0.02 |
| plain MLP (v2 protocol) | 8.92 ± 1.80 | 8.61 ± 1.32 | 5.69 ± 0.56 | 5.45 ± 0.55 | 0.85 ± 0.03 |
| ridge | 9.76 ± 1.70 | 9.98 ± 1.72 | 8.00 ± 0.60 | 8.29 ± 0.68 | 0.86 ± 0.02 |
| HGB | 8.54 ± 1.11 | 8.53 ± 1.77 | 4.92 ± 0.28 | 5.32 ± 0.63 | 0.83 ± 0.06 |
| kNN | 16.57 ± 2.37 | 15.09 ± 1.84 | 7.77 ± 0.51 | 8.70 ± 0.24 | 0.81 ± 0.04 |

Paired over the fresh draws (better / worse draws, t over draws): v2 − v1
OOD −6.41 (6 / 0, t −2.6), ID −9.87 (6 / 0, t −3.7); v2 − HGB OOD −1.66
(5 / 1, t −2.8), ID −0.93 (6 / 0, t −10.0); v2 − plain OOD −2.04 (6 / 0, t
−5.1), ID −1.70 (6 / 0); v2 − ridge OOD −2.88, ID −4.01 (6 / 0 each).
Configuration-balanced OOD per draw, v1 / v2: 14.7 / 7.7, 7.6 / 6.4, 9.4 /
5.0, 7.2 / 6.3, 21.9 / 6.3, 17.3 / 6.6. LOCO-calibrated MPSP intervals (the
v2 LOCO scores of seed 0 without each draw's OOD configurations, q ≈ 12.2)
cover the fresh OOD proxy at 0.88–0.99 (mean 0.95; TAL_glucose 0.70 in one
draw) and the ID proxy at 0.97–1.00. Random CV, fold seed 1: v1 4.12 ±
0.10, v2 3.95 ± 0.03 % (−0.16, seed t −4.4; 4 of 5 folds better, t over
folds −2.9); plain 5.50, ridge 8.11, HGB 4.77, kNN 8.10. No ∂MPSP/∂r > 0
violation for v2 (v1 up to 6.1 %). So the proxy and random-CV conclusions of
4a / 4e hold on fresh draws and folds, with similar effect sizes.

#### 4f. Decisions and open issues

Together (4e) the changes improve LOCO at the configuration level (35 of 37
folds, t over folds −4.4), the productivity extrapolation (36 of 36
configurations) and the dev proxy of the frozen test (6 of 6 draws, OOD
11.8 → 7.1 %, ID 14.3 → 3.8 %; confirmed on fresh draws, 4e). The LOFO and
LOPO fold means fall a lot (64 → 35 %, 86 → 50 %), but the seed t-values
(−24, −7.5) measure seed noise only: over folds the gains are driven by one
fold each and are not statistically established (t −1.3 / −1.5). Held-out
groups that got worse: succinic in leave-one-product-out (28 → 51 %) and
cornstover in leave-one-feedstock-out (35 → 43 %). All v2 CV and proxy
figures are post-selection estimates on the folds and draws used to design
the changes (intro), so they are optimistic, most of all LOFO (designed from
the corn failure). Adopted (each is a protocol correction or a hard
physical guarantee; the accuracy benefit is stated with its fold-level
strength):

* the training protocol V3 for every fit and the out-of-fold LOCO calibration
  of the final intervals (4a; protocol corrections; the test split is never
  used to calibrate);
* the recovery blind to the reagent prices, preset 'physics_v2' (4b; hard
  constraint dg/dprice = 0, CV benefit in the right direction but within the
  noise);
* the hard productivity structure (4c; dMPSP/dr ≤ 0 by construction,
  productivity CV −0.9 pp, seed t −4.5, 28 of 36 configurations better);
* the sugar-availability prior with feedstock-input clipping (4d; the held-out
  feedstock channel 23 → 3.5 % of the numerator; LOFO −29 points, seed t
  −5.3, over the 4 folds t −1.5, corn-driven) and the neutralization gate
  (4d; hard constraint, verified on the data);
* equal ensemble sizes for CV / final / release (the scores divide by the
  member spread), LOCO of the PI-NN and plain MLP on every configuration, the
  ablations on 12 stratified configurations.

Not adopted: refit after early stopping (V2t / V2a / V2m: better than V1,
worse than V3; V3r: no gain at twice the cost) and the 'anneal' schedule (the
code stays, default off); `INNER_GROUP['product'] = 'config_id'` (superseded by
sample-level early stopping everywhere); feedstock, fermentation, scale or
titer-domain inputs for g (the titer-domain feature stays as an opt-in
`FeatureSpec(domain_features=True)`), w_g = 5; a learned β, product groups for
the electricity head, structuring other channels; re-enabling `fs_eta_bound`,
data-derived bounds on z or on the channel outputs, default hyperparameters,
distance-based shrinkage; the configuration-weighted ('group') quantile as the
default (0.81 vs 0.78 MPSP coverage on the proxy, within the noise; available
as `weighting='group'`); random-CV calibration for in-distribution queries
(stored for reporting only: it would need a rule such as the AD distance to
decide when a query is in distribution).

Open issues:

1. Interval coverage per configuration is uneven. The adopted calibration
   (LOCO out-of-fold scores of the same procedure, sample weighting) applied
   to the 4e v2 proxy models (draws 0–5, calibration = the v2 LOCO scores of
   seed 0 without each draw's OOD configurations, q ≈ 12.2) covers the OOD
   proxy at 0.86–0.99 per draw (mean 0.93) and the ID proxy at 0.95–0.99
   (conservative), but TAL_glucose / TAL_corn at 0.70 in 2 of 6 draws and
   lactic_lowpH_SHF_cornstover at 0.58 in one. (The 0.78–0.81 OOD-proxy
   coverage of 4a belongs to the 8-fold configuration-grouped calibration of
   each draw's training rows, q ≈ 8.6, not to the LOCO calibration.) The
   cross-fit restricted to the four OOD families, family-balanced, varies by
   seed: 0.96 / 0.97 / 0.83 (seeds 0 / 1 / 2; TAL 0.87 / 0.95 / 0.75, lactic
   0.98 / 0.91 / 0.64). Untested remedies: scores without the member spread
   (eps only, relative errors; the score tables now keep the raw |error| and
   spread) and Mondrian (per-family) calibration.
2. The final hyperparameters are re-selected under the new protocol in the v2
   run; the round-4 diagnostics (and 4e) used the v1-selected ones.
3. The 300-epoch limit binds for a few members under sample-level early
   stopping (1 of 18 on the proxy; 13 of 504 v2 members in 4e ran to 300
   epochs, one with its best epoch at 300; the report counts them per stage).
   400 epochs would cost little.
4. Early stopping on 'data' vs 'mpsp' under sample-level hold-outs is left to
   the hyperparameter search (dimension `early_stop`).
5. Electricity (with the coproducts) remains the largest error channel (about
   7 % of the numerator with the configuration in training, ~20 % in LOCO); its
   inputs need their own study (CHP / boiler inputs, coupling with natural gas).
6. Data / specification: the saccharification mode (SHF / SSCF, −5 % lactic g)
   is not in the regime specification; the D8 TAL substrate basis; TAL's ~10 %
   recovery noise floor from an unobserved `bp_*` parameter.
7. The g loss weighs absolute errors and the g prior is centred on the
   sample-weighted median (0.98, dominated by ethanol / oxalic / IBO); a
   log-space g loss or a configuration-balanced centre is untested.
8. The recovery of a new separation train is not predictable from the present
   descriptors (LOPO g errors −23 to +44 % for every input set): a
   user-supplied recovery or a richer separation vocabulary is needed.
9. Sugar-availability prior (4d): s = `p_saccharification_yield` is the
   glucan conversion in the data (adapters: enzymatic glucan-to-glucose,
   saccharification reactor R303; xylan is converted in pretreatment), but
   the prior applies it to every structural polysaccharide. Within
   configurations the data z = log(b′_fs w_sugar) follow the prior's
   log(1/s_avail) with slope 0.57 (ethanol corn stover, n 98), 0.72 / 0.53
   (lactic low-pH / neutral SHF), 1.08 (lactic low-pH SSCF), not 1; the
   level is right (median Δz +0.03…+0.07 for corn stover). Inside the
   training range of s (0.77–0.94) the residual, which sees s, absorbs it;
   s is clipped as a 'pretreatment' input but the prior uses the raw value,
   so far below the range the prior's full 1/s slope applies (corn stover at
   s = 0.6: b′_fs about +10 %). Fix (needs a paired dev check): s on the
   hexosans only, a separate pentosan conversion (default 0.9). Not
   relevant to the CV or the frozen test (s in range).
10. Reagent prices still reach the fermentation and downstream capital
   sub-networks (preset 'physics_v2' removed them from g only; a v1 legacy:
   the old 'neutralization' group held the prices). Review measurements
   (round-4 scratch `review_physics/`): the data show no capital
   response (within-configuration correlation of the capital a′ with the
   price −0.01 / +0.03 over 26 / 21 configurations), but a trained model
   moves MPSP through capital by a median 0.66 % of the numerator from the
   q05 to the q95 acid price (4× the change through the neutralization
   channel that buys the acid). Small for accuracy, spurious for price
   sweeps. Follow-up: 'reagent_prices' for the neutralization (and
   materials_other, which holds the pretreatment acid) channels only, with a
   paired LOCO / LOPO check (≥ 3 seeds) before adoption.
11. Scale response (final review, round 4; REPORT.md section 8 (ii-c)): each
   configuration's capacity is sampled only 0.81-1.19 x its baseline, and
   the release model reproduces a median 0.60 / 0.55 of the simulated MPSP
   decrease of a 2x / 3x plant (BioSTEAM at baseline parameters,
   `data_generation/checks/capacity_response.csv`). The user API flags a
   capacity outside the configuration's training range
   (`contextual_outside_range`), but the model is unchanged. Proposed for
   the next data and model round (needs a re-run): log-uniform capacity
   sampling over 0.3-3x the baseline in the campaigns, and/or a hard
   capital-scaling structure (the capital and capital-proportional fixed
   channels proportional to (F/F_ref)^(β−1), β learned in [0.5, 1]),
   validated against the simulated capacity response.
12. Interval calibration per family (final review, round 4; REPORT.md
   section 6, "Coverage for a KNOWN configuration"): the pooled LOCO
   quantiles over-cover HP, ethanol, isobutanol, oxalic and succinic acid
   and under-cover TAL's g (also for known configurations: 0.45 on the
   random-CV scores; the decarboxylation-conversion noise of item 6) and
   lactic's d′ (0.56). Candidates (need a re-run and a fresh calibration):
   Mondrian (per-family) conformal quantiles with a fallback for small
   families, or a score with an aleatoric term (a per-family floor ε, or a
   variance head for g).

Review of round 4 (fixes that change no trained model; the evidence above
is unchanged):

* Statistics and wording: fold-level statistics next to every seed t (4e);
  the post-selection caveat (intro, 4e, 4f) and the confirmation on fresh
  proxy draws 6–11 and random-CV fold seed 1 (4e); the report and READMEs
  call the v2 test an adaptive second use of the split (likely optimistic,
  v2 released whatever the v1 / v2 test comparison shows), qualify "the CV
  estimates the procedure" (HPO not nested), describe the coverage as
  marginal over the samples of new configurations (sample weighting),
  compare realized pooled / configuration-mean test coverage with the
  cross-fitted pooled / mean estimate, drop the unsupported "slightly
  conservative" claim for the release intervals, and note that the v1
  final-test baselines were fit without the 6 v1 calibration configurations.
* Intervals: the out-of-fold score tables keep the raw |error| and member
  spread, and `calibrate_conformal_oof` re-scores a′–d′ and g with the final
  ensemble's eps (§6); `<key>_level` of an out-of-fold calibration is
  labelled nominal (`level_kind`).
* Pipeline: an unusable plain-MLP out-of-fold table leaves the plain
  `final_dev` without intervals instead of stopping the run (the PI-NN's
  stays strict); plan.json is reused only while the requested budget, search
  size, HPO members and ablation folds are unchanged; `--ablation-loco-folds`
  is honoured in quick mode (default 4 there, 12 in full mode); the
  `--v1-protocol` report labels final_dev's hold-out calibration 'holdout';
  a quick report no longer compares itself with the real v1 tables; a
  full-mode smoke test (spawned pool, resume, re-plan) joins the suite.
* Features: a feedstock type outside the vocabulary (e.g. 'lignocellulosic')
  used to count as a known non-lignocellulosic type (s_avail = 0.05 floor, a
  20× feedstock-intensity prior); it is now unknown (prior not applied,
  imputed 1 and reported), after case / blank normalization. The prior is
  identical on all development rows.

#### v2 run vs. the round-4 estimates (written after the v2 run)

The v2 pipeline re-ran the HPO under the v2 protocol and selected different
hyperparameters (width 64, depth 2, lr 0.0045, `w_ch` 0.1, `w_phys` 1.0,
early stopping on MPSP; `results/tables/selected_hyperparameters.json`) from
those of v1 used in the round-4 estimates (width 64, depth 4, lr 0.0012,
`w_ch` 1.0, `w_phys` 0.03; `results/v1/tables/selected_hyperparameters.json`).
The released CV numbers therefore differ from 4e:

* fold means, v1 run → 4e (v2 settings, v1 hyperparameters) → v2 run: LOCO
  14.4 → (4e: 35 of 37 folds better) → 9.3 %; LOFO 66.3 → 35 → 23.9 %; LOPO
  75.1 → 50 → 57.2 %;
* LOPO per held-out product, v1 run → 4e → v2 run (%): ethanol 177.1 → 44.4
  → 59.4, isobutanol 148.2 → 85.1 → 66.2, **HP 21.8 → 12.7 → 77.8**,
  lactic 59.6 → 53.5 → 54.3, oxalic 47.5 → 59.3 → 41.1, **succinic 28.3 →
  51.4 → 59.7**, TAL 43.0 → 44.3 → 41.8. The median over products got
  worse (v1 47.5 → v2 59.4 %) although the mean fell (75.1 → 57.2 %, carried
  by ethanol and isobutanol).

**The HP fold is a regression of the v2 run**: 77.8 % against 21.1 (plain
MLP), 21.0 (ridge), 21.3 (gradient boosting) and 29.5 % (kNN) on the same
fold, and 79-82 % for all three PI-NN ablations, so the cause is shared by
every PI-NN variant of the v2 run, plausibly the re-selected hyperparameters
(the ablations use them too, and 4e with the v1 hyperparameters gave
12.7 %). Succinic got worse already in 4e (protocol, 4d) and again in the v2
run. Diagnosed in round 5 (below, 5d-5e): the cause is the HPO's choice
`clip_inputs = 0`, which clips the held-out HP's final-product
stoichiometric ratio (0.80) to the training minimum 1.00 (v2-selected hyperparameters with only
the clipping switched off: HP 12.9 %); the exact stoichiometric ratios are
now never clipped (`experiments.CLIP_EXEMPT`). The "All changes together" statement of REPORT.md
("Changes since v1": LOPO 86.1 → 50.1, 6 of 7 better) refers to 4e, not to
the released run.

### Development-set diagnostics, round 5: channel / price-intensity fidelity and the v2 LOPO regression (2026-10-17; never on test data)

Data: a frozen copy of the processed tables of the v2 run (data hash
51d0efe2b4f517a4, 1550 samples / 41 configurations), restricted to the
development set with `development_set` (1306 samples / 37 configurations);
no diagnostic or decision used a test sample, `models/final_dev` or
`models/release`. Settings: the v2 pipeline (`make_job` with the water gate,
masking 'physics_v2', r-structure 'physics', EXTRAPOLATION_SAFEGUARDS,
downstream FeatureSpec, PROTOCOL v2) and the v2-selected hyperparameters,
unless stated otherwise. With every new option off, the round-5 code
reproduces the v2 code's predictions exactly (max difference 0.0 on the HP
LOPO fold).

**5a. Fidelity metrics** (`evaluate.channel_fidelity`, out of fold,
configuration-balanced = mean over held-out configurations of their sample
means): `ch_<c>` = RMS over the sample's simulated points of the error of
channel c's numerator, relative to the RMS target total numerator (price
channels as price × intensity); `ch_rms` = root sum of squares over the 14
channels; per price key the predicted sensitivity dMPSP/dp = (q_a + q_b/y +
q_c/t + q_d/(yt))/ĝ against the simulated per-point intensity (points
`q__fs`, `q__el`, `q__ng`, exact at fixed design): `sens` = relative RMS
error, `bias` = median ratio − 1 (signed; absolute in averages), `resp` = RMS
of (dMPSP/dp predicted − simulated)·p/MPSP (the price response in % of
MPSP). Relative electricity / NG errors are reported only where the
simulated price share is material (≥ 2 % of MPSP: the net electricity of
exporters crosses zero); `*` = without the TAL and lactic configurations
whose feedstock targets are masked in option (2).

**5b. Data checks.**

* *Electricity without an electricity coproduct* (ethanol corn, three lactic
  configurations; 177 samples, 3667 points): the net electricity intensity
  is never negative (minimum 0.086 kWh/kg); a per-coefficient non-negative
  fit costs ≤ 0.09 % of MPSP (OLS 0.03 %). The power gate (§3) is therefore
  a valid hard constraint. Isobutanol lists electricity as a coproduct but
  its CHP only covers its own demand: q_el = 0 at every point (no gate
  applies; only the loss can fix its channel).
* *Feedstock intensity.* Within a sample F_feedstock/F_sugar_fed is constant
  (CV 0), so g·q_fs·y = (g/rec_i)·F_fs/F_sugar/basis: the non-b′ structure
  of the primed intensity comes from (i) the per-point recovery rec_i
  varying around the sample's median g (CV HP 2.5–5 %, TAL 6 %, lactic
  8–11 % with correlation 0.81–0.84 with y·t, succinic / ethanol / oxalic
  ≤ 0.6 %) and (ii) TAL's basis (D8: F_ferm/(y·F_sugar) has CV 0.13, up to
  1.64). A pure b′/y fit of g·q_fs is within ≤ 5 % RMS of every sample's
  intensity (lactic 6–9 %), and the product-blind b′ of each feedstock
  (median over the non-TAL, non-lactic samples: corn 1.486, corn stover
  2.232, glucose 1.124, sugarcane 6.783) matches each configuration's own
  pure-b′ fit within ±2–4 %, TAL included; only lactic (+5 … +20 %) and
  TAL_SA_sugarcane (−10 %) deviate. (TAL_SA_sugarcane's deviation is data
  defect D16, found later: its cane mill receives the juice recycle twice,
  so its data feed 15 % more sugar per kg cane than the cane supplies; its
  corrected b′ is 0.327, 4 % above TAL_sugarcane's, as for the other
  feedstocks.) TAL's c′ = −17 … −90 is mostly an
  ill-conditioned OLS split (4-coefficient fit residual 0.6–2.3 %, pure
  b′ 3.6–4.6 %). So lactic's feedstock targets are not defective data (D1 was fixed
  by regeneration, D9 is harmless): its d′ of 44–56 is the recovery shape
  plus the b′–d′ collinearity of its domain (y ≥ 0.3, t ≥ 40 g/L), which a
  product-blind feedstock head cannot represent either way.
* *In-sample infidelity of v2* (dev-trained, 2 members, 3 seeds, evaluated
  on its own training rows): every corn configuration's feedstock
  sensitivity is 19–27 % low (also ethanol and isobutanol, whose own
  targets are pure b′), isobutanol's electricity response is off by 94 % of
  MPSP and its feedstock channel by 72 % of the numerator, corn ethanol's
  electricity by 32 %, while MPSP fits at 2–8 %: with the coefficient-space
  channel term at w_ch = 0.1 the MPSP fit leaves the decomposition free
  (§4: a 25 %-of-MPSP electricity error costs 0.005 of loss).

**5c. Options** (paired over seeds 0–2 on the same folds, one member per
fit; LOCO = the 12 stratified LOCO folds of the pipeline ablations
(`stratified_subset`) for the screen, all 37 for the finalists below; LOPO
= 7 folds; LOFO = 4 folds; value (difference to v2, t over seeds); %):

* (1) channel / intensity term in surface space (`LossConfig.ch_space =
  'surface'`) with w_ch ∈ {1, 3, 10}; the total keeps coefficient space;
* (2) masking the feedstock-channel targets of TAL and lactic in the loss
  (`train.mask_channel_targets`);
* (3) the power gate (§3);
* (4) a pure-b′ feedstock intensity (`ModelConfig.fs_b_only`, motivated by
  5b), and combinations.

LOPO:

| variant | MPSP fold | MPSP cfg | ch_rms | ch_feedstock | ch_electricity | sens_fs* | resp_fs | resp_el | resp_ng |
|---|---|---|---|---|---|---|---|---|---|
| v2 | 59.1 | 66.5 | 41.0 | 10.4 | 13.1 | 22.9 | 12.3 | 16.7 | 12.7 |
| (1) w_ch 1 | 60.6 (+1.5, +0.9) | 66.8 (+0.3, +0.5) | 34.2 (−6.8, −23.6) | 2.9 (−7.5) | 10.3 (−2.8) | 19.5 (−3.4) | 9.3 (−3.1) | 14.2 (−2.5) | 13.5 (+0.8) |
| (1) w_ch 3 | 59.5 (+0.4, +0.3) | 62.4 (−4.1, −2.8) | 33.2 (−7.9, −12.7) | 2.4 (−8.0) | 10.5 (−2.6) | 17.6 (−5.3) | 8.5 (−3.8) | 13.3 (−3.4) | 12.2 (−0.6) |
| (1) w_ch 10 | 61.2 (+2.1, +4.6) | 65.4 (−1.0, −0.6) | 34.4 (−6.6, −6.1) | 2.3 (−8.1) | 10.3 (−2.8) | 20.5 (−2.4) | 9.2 (−3.1) | 13.1 (−3.5) | 14.1 (+1.4) |
| (2) mask | 60.9 (+1.8, +2.9) | 66.7 (+0.2, +0.2) | 40.7 (−0.3, −0.3) | 7.3 (−3.1) | 13.9 (+0.8) | 20.2 (−2.7) | 10.8 (−1.5) | 16.9 (+0.2) | 13.0 (+0.3) |
| (3) gate | 54.9 (−4.3, −5.5) | 62.5 (−3.9, −2.3) | 39.5 (−1.6, −2.5) | 10.5 (+0.2) | 12.9 (−0.2) | 21.1 (−1.8) | 11.8 (−0.5) | 16.9 (+0.2) | 12.3 (−0.4) |
| (4) b′ only | 62.3 (+3.2, +1.6) | 68.0 (+1.5, +1.0) | 42.9 (+1.9, +2.2) | 9.9 (−0.5) | 13.5 (+0.4) | 19.2 (−3.7) | 11.7 (−0.6) | 16.0 (−0.7) | 14.8 (+2.1) |
| (1, w 3) + (3) | 54.8 (−4.3, −3.5) | 62.6 (−3.8, −7.9) | 33.3 (−7.7, −10.7) | 2.6 (−7.8) | 11.2 (−1.9) | 18.7 (−4.2) | 8.8 (−3.5) | 15.4 (−1.3) | 12.8 (+0.0) |
| (1, w 3) + (4) | 59.2 (+0.1, +0.3) | 64.2 (−2.3, −3.6) | 37.5 (−3.6, −17.3) | 2.0 (−8.4) | 10.9 (−2.2) | 17.6 (−5.3) | 8.3 (−4.0) | 13.7 (−3.0) | 13.9 (+1.1) |
| (1, w 3) + (3) + (4) | 56.7 (−2.4, −2.4) | 65.0 (−1.5, −1.4) | 36.9 (−4.1, −2.9) | 2.0 (−8.4) | 11.8 (−1.3) | 18.2 (−4.7) | 8.7 (−3.7) | 15.8 (−0.8) | 13.8 (+1.1) |
| (1, w 3) + (3) + (2) | 54.9 (−4.2, −4.1) | 64.3 (−2.2, −2.1) | 35.6 (−5.4, −4.2) | 2.3 (−8.1) | 11.8 (−1.3) | 19.4 (−3.5) | 8.9 (−3.4) | 15.9 (−0.8) | 13.5 (+0.8) |

LOFO:

| variant | MPSP fold | MPSP cfg | ch_rms | ch_feedstock | ch_electricity | sens_fs* | resp_fs | resp_el | resp_ng |
|---|---|---|---|---|---|---|---|---|---|
| v2 | 26.9 | 22.0 | 39.0 | 11.8 | 18.2 | 17.5 | 9.2 | 16.5 | 7.0 |
| (1) w_ch 1 | 27.5 (+0.6, +0.2) | 22.7 (+0.7, +0.4) | 37.9 (−1.1, −1.1) | 9.1 (−2.8) | 17.6 (−0.5) | 13.2 (−4.3) | 7.5 (−1.6) | 15.7 (−0.8) | 6.1 (−0.9) |
| (1) w_ch 3 | 26.4 (−0.5, −0.2) | 22.5 (+0.5, +0.3) | 36.6 (−2.4, −2.7) | 9.9 (−2.0) | 17.7 (−0.4) | 14.0 (−3.5) | 7.6 (−1.5) | 15.9 (−0.6) | 5.7 (−1.3) |
| (1) w_ch 10 | 28.5 (+1.5, +0.8) | 23.6 (+1.6, +1.0) | 36.3 (−2.7, −8.7) | 9.2 (−2.6) | 17.9 (−0.2) | 13.8 (−3.7) | 7.7 (−1.5) | 16.1 (−0.3) | 5.3 (−1.7) |
| (2) mask | 27.9 (+1.0, +0.6) | 22.3 (+0.3, +0.3) | 37.0 (−2.1, −2.8) | 6.2 (−5.6) | 18.2 (+0.0) | 10.9 (−6.6) | 6.5 (−2.7) | 16.6 (+0.2) | 6.6 (−0.4) |
| (3) gate | 33.3 (+6.4, +7.1) | 24.1 (+2.1, +7.6) | 37.2 (−1.9, −1.8) | 9.3 (−2.5) | 14.8 (−3.3) | 14.3 (−3.2) | 8.0 (−1.2) | 14.9 (−1.6) | 6.4 (−0.5) |
| (4) b′ only | 25.7 (−1.2, −0.6) | 21.5 (−0.5, −0.4) | 37.3 (−1.7, −2.9) | 6.6 (−5.2) | 18.2 (+0.1) | 12.4 (−5.1) | 6.3 (−2.9) | 16.3 (−0.1) | 7.2 (+0.2) |
| (1, w 3) + (3) | 35.0 (+8.1, +9.7) | 25.3 (+3.3, +6.4) | 35.8 (−3.2, −7.6) | 10.2 (−1.6) | 14.9 (−3.3) | 14.9 (−2.6) | 8.4 (−0.8) | 14.7 (−1.8) | 5.5 (−1.5) |
| (1, w 3) + (4) | 25.7 (−1.2, −0.5) | 22.0 (−0.0, −0.0) | 34.8 (−4.3, −4.1) | 2.7 (−9.2) | 17.4 (−0.8) | 5.3 (−12.2) | 4.0 (−5.2) | 15.6 (−0.8) | 5.8 (−1.1) |
| (1, w 3) + (3) + (4) | 33.0 (+6.0, +1.9) | 24.6 (+2.6, +1.4) | 33.9 (−5.2, −4.0) | 3.2 (−8.6) | 14.8 (−3.4) | 6.2 (−11.3) | 4.5 (−4.6) | 14.8 (−1.7) | 5.8 (−1.2) |
| (1, w 3) + (3) + (2) | 32.9 (+6.0, +2.1) | 24.4 (+2.4, +1.4) | 33.8 (−5.2, −6.9) | 4.2 (−7.6) | 14.9 (−3.2) | 7.1 (−10.4) | 4.9 (−4.3) | 15.0 (−1.5) | 5.6 (−1.4) |

LOCO, all 37 folds (finalists):

| variant | MPSP fold = cfg | MPSP geo-mean | ch_rms | ch_feedstock | ch_electricity | sens_fs* | bias_fs* | resp_fs | resp_el | resp_ng |
|---|---|---|---|---|---|---|---|---|---|---|
| v2 | 9.0 | 5.2 | 25.7 | 9.9 | 10.8 | 15.2 | 9.1 | 7.8 | 10.4 | 4.6 |
| (1, w 3) + (3) | 8.5 (−0.5, −0.8) | 5.4 (ratio 1.04: 1.06 / 1.00 / 1.06) | 16.2 (−9.5, −12.3) | 2.3 (−7.6) | 6.9 (−3.9) | 3.9 (−11.3) | 1.5 (−7.6) | 2.0 (−5.7) | 7.8 (−2.6) | 3.3 (−1.2) |
| (1, w 3) + (3) + (2) **adopted** | 8.6 (−0.4, −0.6) | 5.6 (ratio 1.08: 1.13 / 1.01 / 1.09) | 16.5 (−9.2, −9.7) | 2.2 (−7.7) | 7.1 (−3.7) | 3.8 (−11.4) | 1.5 (−7.7) | 2.1 (−5.6) | 8.0 (−2.4) | 3.4 (−1.2) |
| (1, w 3) + (3) + (4) | 9.1 (+0.1, +0.2) | – | 17.4 (−8.3, −9.9) | 2.1 (−7.7) | 7.2 (−3.6) | 3.4 (−11.8) | 1.5 (−7.6) | 2.4 (−5.4) | 8.1 (−2.3) | 3.5 (−1.1) |

In-sample (all development rows, 2 members, 3 seeds): v2 → adopted MPSP
3.4 → 3.3 %, ch_rms 21.2 → 5.7 %, feedstock channel 8.9 → 1.6 %,
electricity 9.5 → 2.2 %, sens_fs* 13.4 → 2.5 %, resp_el 9.6 → 2.9 %, resp_ng
4.1 → 1.8 %.

Findings:

* *(1) is the fix.* Comparing channels where their economic size is
  measured, at the sample's own points, removes the coefficient-space
  failure modes: the corn feedstock sensitivity bias of every corn
  configuration (LOCO, per configuration median: v2 −34 … −18 %; adopted
  −1.2 … +4.4 %), isobutanol's feedstock response (LOCO 95 → 7 % of MPSP)
  and most of the electricity error. w_ch = 3 is chosen on fidelity:
  w_ch 1 → 3 improves every scheme; w_ch 10 improves LOCO slightly more but
  worsens the LOPO feedstock sensitivity (20.5 vs 17.6 %) and NG response,
  i.e. fidelity saturates at 3 where it matters for extrapolation. MPSP:
  fold and configuration-balanced means are unchanged or better in every
  scheme; the typical in-family configuration is slightly worse on LOCO
  (geometric mean 5.2 → 5.4–5.6 %, median +0.25 points; TAL / succinic /
  ethanol +0.5–2 points, HP ±0.5), but this statistic moved
  non-monotonically with w_ch on the 12-fold screen (ratios 1.04 / 1.06 /
  0.95 for w_ch 1 / 3 / 10) and between seeds (0.90–1.13), i.e. it is at the
  noise level of single-member fits. To check in the next pipeline run (5
  members).
* *(3) is valid physics and fixes a qualitative failure, at one
  configuration's cost.* Without the gate 41–100 % of the held-out samples
  of the non-exporting configurations get a negative (export) electricity
  sensitivity (v2: 98–100 % in every scheme; predicted p_el·dMPSP/dp_el /
  MPSP −0.7 … −42 % against +2 … +10 % simulated); with it 0 % by
  construction and lactic's response error falls from 8–14 to 1–5 % of
  MPSP. But in the corn LOFO fold ethanol corn's MPSP rises (v2 33 %, gated
  92–106 %; LOCO 14 → 21 %): the ungated model predicted power export
  (electricity −44 % of the numerator) against a missed DDGS / corn-oil
  credit (coproducts +61 %) and capex_utilities +30 %; the gate removes the
  compensating wrong-sign term and exposes the coproduct error. Without
  ethanol corn the adopted variant's LOFO configuration-balanced MPSP
  change is +0.6 points (with it +2.4; fold mean +6.0, carried by the corn
  fold). Isobutanol (electricity listed as coproduct, net electricity 0,
  ungated, the only configuration of its product) trades electricity
  against other channels between variants (LOCO MPSP 59 → 37 %, electricity
  response 78 → 109 % of MPSP), which dominates the averaged electricity
  metrics in LOCO and LOPO.
* *(2) helps new feedstocks.* With (1) + (3), masking TAL's and lactic's
  feedstock targets reduces the LOFO feedstock-channel error 10.2 → 4.2 %
  and the feedstock sensitivity error 14.9 → 7.1 % (MPSP −2.1); it is
  neutral on LOCO and costs LOPO a little (configuration-balanced +1.7, t
  2.1; ch_rms +2.3). Alone, in coefficient space, it is weaker than (1)
  (LOCO feedstock channel 12.8 → 8.7 % vs 2.2 %). The masked
  configurations' own feedstock sensitivity becomes product-blind (LOCO
  bias lactic −6 … +29 %, TAL +0.2 … +16 %), the price of keeping their
  unrepresentable shapes out of the shared feedstock heads.
* *(4) is not adopted.* With (1) it gives the best new-feedstock
  sensitivities (LOFO sens_fs* 5.3 %) but the recovery-shaped feedstock
  terms the data do have (5b) move into other channels: LOPO ch_rms +4.3
  (t 10.3) vs (1) alone, full LOCO MPSP +0.6 and ch_rms +1.2 vs (1) + (3).
  It stays a library option (`ModelConfig.fs_b_only`, CLI `--fs-b-only`).
* The channel-supervision ablation of the v2 run (`pinn-no_channels`, LOCO
  −0.2 pp: "the intermediate data add little") measured the
  coefficient-space term at w_ch = 0.1; in surface space the decomposition is
  actually learned: in-sample channel errors of 1–2 % of the numerator
  (feedstock 1.6, electricity 2.2, NG 1.4, coproducts 1.1 %) instead of
  4–10 %, at unchanged MPSP. The MPSP gain from the intermediate data
  remains small; their value is the faithful breakdown and price response.

**Proposals of the study (pipeline, `experiments.py`; library defaults
unchanged; each with a CLI flag restoring v2, all at once `--v2-model`;
decided in 5f after a re-check on the corrected data):** (1) with w_ch = 3
as `CHANNEL_LOSS`, a fixed setting chosen on fidelity, not an HPO dimension
(`hpo_candidates(fixed=...)` keeps the other candidate dimensions;
`--v2-channel-loss`); (3) `POWER_GATE` (`--no-power-gate`); (2)
`DEFECTIVE_TARGETS = {'feedstock': ('TAL', 'lactic')}` (`--no-target-mask`);
(4) opt-in only. The CV stage now stores the PI-NN's out-of-fold channels
and intensities and reports the 5a metrics per scheme and per held-out
configuration (`tables/cv_fidelity_summary.csv`,
`cv_fidelity_per_config.csv`, `cv_channel_fidelity.csv.gz`, REPORT section
5). Open: the coproduct credit of an unseen starch feedstock (ethanol corn
LOFO), isobutanol's electricity (net zero in the data, no structural
anchor), and the recovery-shaped part of the feedstock intensity (a
product-specific recovery surface, not a constant g, would be needed).

**5d. The v2 leave-one-product-out regression** (all seven LOPO folds,
2 members, seeds 0–2, v2 settings; point MAPE of the held-out product, mean
over seeds ± seed sd):

| hyperparameters | ethanol | HP | IBO | lactic | oxalic | succinic | TAL | fold mean | median |
|---|---|---|---|---|---|---|---|---|---|
| v2-selected (candidate 7) | 61.6 ± 5.4 | **82.0 ± 3.8** | 65.6 ± 10.3 | 53.5 ± 1.6 | 50.9 ± 4.4 | 57.5 ± 2.4 | 42.6 ± 0.9 | 59.1 ± 1.1 | 57.5 |
| v2-selected, `clip_inputs` 0 → −1 | 61.6 ± 5.6 | **12.9 ± 0.9** | 75.2 ± 11.6 | 48.4 ± 0.5 | 85.9 ± 7.6 | 58.3 ± 2.7 | 38.9 ± 1.1 | 54.5 ± 1.4 | 58.3 |
| v1-selected (candidate 3) | 45.8 ± 5.1 | 12.8 ± 0.6 | 89.6 ± 10.5 | 51.7 ± 3.0 | 94.4 ± 55.3 | 41.9 ± 5.7 | 45.3 ± 1.0 | 54.5 ± 8.0 | 45.8 |
| package defaults (candidate 0) | 55.0 ± 4.6 | 16.5 ± 2.9 | 83.7 ± 21.2 | 49.7 ± 1.1 | 61.3 ± 15.8 | 54.3 ± 2.8 | 41.4 ± 0.9 | 51.7 ± 5.2 | 54.3 |

The v2 run's HP fold (77.8 %) is reproduced (82.0 %), and changing only
`clip_inputs` from 0 to −1 in the v2-selected hyperparameters removes it
(12.9 %; seed sd 0.9): the HP regression is the input clipping chosen by
the v2 HPO, not the protocol or the other dimensions. Mechanism: with HP
held out, only three continuous encoder inputs of HP lie outside the
training range, the final-product stoichiometric mass ratio and its MW-ratio
twin (HP → acrylic acid 0.80; in training 1.00 for every product sold
without a mass-changing conversion and 1.19 for TAL_SA → sorbic acid, so
HP lies 4.1 training sd below the minimum) and the final product's
mass-specific heat of vaporization (1.3 sd). Clipping sets HP's ratio to
1.00 while its upgrading and final-product descriptors say dehydration to
acrylic acid: in training an upgrading conversion only occurs with the
ratio 1.19, so the combination is new, and the sub-networks fed by the
'final_product' group (materials_other, capex_downstream, the natural-gas
intensity, g) extrapolate badly: signed contributions to HP's relative
numerator error materials_other +22.0, capex_downstream +15.3, natural gas
+11.0 points (unclipped: +5.0, −4.5, +4.3), and ĝ 0.56 vs 0.66 (unclipped
0.69); MPSP bias +82 % (unclipped +5 %). The same clipping helps oxalic,
whose product and final-product descriptors lie 2–9 sd outside the range
(pKa 7.0 sd, electron-based maximum yield 8–9 sd, final-product mass
fraction 7.5 sd): 85.9 → 50.9 % with clipping, mostly through g (ĝ 0.80 →
0.95, truth 1.00). Succinic is not affected by the hyperparameters' clipping
(ĝ 0.70–0.71 vs 0.98 under all four rows; its separation train resembles
TAL's and lactic's, whose recoveries are 0.64–0.73, open issue 8 of 4f).

**5e. HPO criterion and search space.** The 16 HPO candidates of the v2
run (same seeded sequence) on the seven LOPO folds (2 members, seed 0, as
the HPO), against their LOCO criterion of the v2 run (same data hash, same
5 stratified LOCO folds; `hpo_candidates.csv`):

* the two criteria are unrelated across candidates (Spearman 0.05); the
  LOCO top three (candidates 7, 15, 1) all clip the inputs, and every one of
  the six candidates with `clip_inputs = 0` has a broken HP fold (24–113 %),
  every one of the ten without clipping a sound one (10–20 %);
* a criterion (1 − w)·LOCO + w·LOPO selects candidate 3 (the v1 choice) for
  every w ≥ 0.25: LOCO criterion 10.1 vs 8.8 % for candidate 7 (+1.3
  points, about one standard error of 1.2–1.4), LOPO criterion 46.1 vs
  59.4 % (−13.3; HP 13 vs 86 %);
* post-selection optimism: candidate 3's selection-seed LOPO fold mean is
  45.8 %, its mean over seeds 0–2 54.5 ± 8.0 % (oxalic 43 % on seed 0, 94 ±
  55 % over seeds): a LOPO-weighted selection inherits the seed noise of
  single-seed LOPO folds, so its LOPO CV estimate is optimistic by several
  points (here ~9) unless re-estimated on fresh seeds;
* clipping margin: `clip_inputs` = 2 or 5 sd does not resolve the
  HP / oxalic trade-off (HP 35.8 / 12.9 %, oxalic 79.1 / 87.0 %; edge
  clipping gives 82.0 / 50.9 %, none 12.9 / 85.9 %): oxalic's benefit needs
  clipping at the range edge;
* exempting the exact stoichiometric ratio (and its MW-ratio twin) from
  clipping (`ModelConfig.clip_exempt`) with the v2-selected
  hyperparameters otherwise unchanged: HP 82.0 → 14.8 % (seed sd 2.6), and
  every other LOPO fold is bit-identical to v2 (oxalic keeps 50.9 %: the
  ratio lies inside the training range for every other held-out group, so
  the clip bounds of the training rows do not change); LOPO fold mean
  59.1 → 49.5 %, median 57.5 → 53.5 %. LOCO is unaffected for the same
  reason.

Cost of the LOPO term: the HPO trains 7 more folds per candidate (training
rows 7836 vs 5224 for the 5 LOCO folds), i.e. 2.5× the HPO time; at the v2
run's size 0.27 → 0.66 h on 4 workers, the whole run 1.89 → 2.29 h (+21 %,
inside the 3 h budget of v2).

**Proposals of the study** (decided in 5f): (i) search space / feature handling: the exact
stoichiometric ratios are never clipped (`experiments.CLIP_EXEMPT`, every
network of the pipeline; `--no-clip-exempt`): it removes the mechanism of
the v2 regression at no cost elsewhere, whatever the HPO selects (adopted);
(ii) criterion: LOPO enters with weight 0.5 (`--hpo-lopo-weight 0.5`; not
adopted, 5f), chosen a priori from the tool's documented uses
(new configurations and new products, README 4.3), not tuned on these
numbers: the LOCO criterion alone is blind to product-extrapolation failures
(its ranking is unrelated to LOPO's), at +21 % run time and, on this search,
about one standard error of the LOCO criterion. Caveats: the 16-candidate
LOPO values are single-seed (seed sd of a LOPO fold mean 1–8 points, of one
product up to 55), so the selection is noisy and the selected candidate's
LOPO CV estimate is post-selection (here ~9 points optimistic on the
selection seed); the next run's REPORT must label its LOPO numbers as such
(like LOCO's, §6 intro), ideally re-estimated on fresh seeds. LOFO is not
added (4 folds, dominated by corn and by the feedstock structures of 4d;
its failure modes are handled by those hard structures).

**5f. Integration on the corrected data, and decisions** (2026-10-18;
development rows of the regenerated tables only: data hash
c345efda41c4fdca, `development_set` = 1353 samples / 37 configurations; no
test sample, `models/final_dev` or `models/release` used). The patch of the
study was applied to the code that already held the data corrections (D15,
H1), and the proposals were re-checked through the pipeline's own job path
(`experiments.make_job` with the `mk_extra` / `lk_extra` / `target_mask`
context that `experiments.main` builds): v2 model settings against the
proposals, v2-selected hyperparameters (w_ch = 3 under the surface loss),
one member per fit, seeds 0–2 paired on the same folds, all 37 LOCO, 7 LOPO
and 4 LOFO folds; metrics as in 5a (`*` = without TAL and lactic; value
(difference, seed t); %). "final" = (1) + (3) + clip exemption + the TAL
mask only (the decision below); "+ lactic mask" = the study's proposal.

| scheme | variant | MPSP fold mean | MPSP cfg | MPSP geo-mean | ch_rms | ch_feedstock | sens_fs* | bias_fs* | resp_fs | resp_el | wrong-sign el. |
|---|---|---|---|---|---|---|---|---|---|---|---|
| LOCO | v2 | 9.7 | 9.7 | 5.6 | 25.2 | 9.9 | 15.1 | 9.7 | 7.8 | 9.8 | 84.9 |
| LOCO | final | 9.3 (−0.4, −0.9) | 9.3 | 5.7 (+0.1, +1.9) | 15.7 (−9.5, −92) | 2.3 | 4.2 (−10.9) | 1.7 | 2.3 | 7.7 (−2.2) | 0 |
| LOCO | + lactic mask | 9.1 (−0.6, −1.3) | 9.1 | 5.7 (+0.1, +0.4) | 15.8 | 2.1 | 3.7 | 1.5 | 2.2 | 8.0 | 0 |
| LOCO | no mask | 9.4 (−0.3) | 9.4 | 5.8 | 15.8 | 2.3 | 4.5 | 1.9 | 2.4 | 7.8 | 0 |
| LOCO | final without the gate | 10.2 (+0.4, +0.8) | 10.2 | 5.6 (+0.0) | 15.2 | 2.1 | 4.0 | 1.6 | 2.2 | 5.1 | 38.4 |
| LOPO | v2 | 57.2 | 55.0 | 52.0 | 38.6 | 10.7 | 19.4 | 13.0 | 11.0 | 14.8 | 100 |
| LOPO | final | 46.5 (−10.7, −4.4) | 36.6 (−18.5, −7.5) | 28.1 | 31.5 (−7.1, −11.7) | 2.4 | 18.1 | 16.5 (+3.4) | 8.0 | 14.0 | 0 |
| LOPO | + lactic mask | 47.1 (−10.1, −3.0) | 38.7 (−16.4, −5.8) | 29.6 | 33.6 | 2.1 | 18.5 | 17.0 | 8.1 | 15.3 | 0 |
| LOPO | no mask | 44.5 (−12.7) | 35.8 | 27.4 | 31.8 | 2.3 | 18.4 | 17.1 | 8.3 | 14.1 | 0 |
| LOPO | final without the gate | 52.0 (−5.2) | 38.9 | 29.5 | 32.3 | 2.4 | 17.3 | 15.5 | 7.8 | 11.0 | 42.0 |
| LOFO | v2 | 30.3 | 23.5 | 17.6 | 38.9 | 11.0 | 16.6 | 8.2 | 8.9 | 16.1 | 100 |
| LOFO | final | 34.1 (+3.8, +2.0) | 24.7 (+1.1, +1.2) | 16.6 (−1.0) | 35.0 (−3.9, −12.6) | 8.5 | 12.7 (−3.9, −5.1) | 6.0 | 7.5 | 14.4 | 0 |
| LOFO | + lactic mask | 35.1 (+4.8, +1.0) | 25.2 | 17.3 | 33.6 | 4.1 | 6.9 (−9.7, −1.8) | 5.2 | 4.8 | 14.5 | 0 |
| LOFO | no mask | 33.6 (+3.3) | 24.6 | 16.9 | 34.6 | 7.9 | 11.9 | 5.8 | 7.2 | 14.4 | 0 |
| LOFO | final without the gate | 29.6 (−0.7) | 23.4 | 16.5 | 36.1 | 9.0 | 13.3 | 6.6 | 7.6 | 15.8 | 98.3 |

Folds (seeds averaged per fold, final vs v2): LOCO 14 better / 23 worse,
median +0.12 points, t over folds −0.40; LOPO 6 / 1 (median −5.2, t −1.45;
HP 54.9 → 13.2, isobutanol 70.5 → 39.7, lactic 55.7 → 49.3, oxalic 52.1 →
47.3, succinic 61.1 → 55.9, TAL 46.4 → 43.3, ethanol 59.7 → 76.7); LOFO 2 /
2 (corn 47.3 → 63.3, corn stover 31.7 → 35.8, glucose 9.5 → 7.2, sugarcane
32.7 → 30.1). Per configuration (LOCO): the corn feedstock bias is −32 …
−16 % under v2 (every corn configuration: HP 3, TAL 2, succinic 2,
ethanol, isobutanol) and −0.6 … +1.0 % with the final settings (isobutanol
+8.7 %); isobutanol 70.5 → 39.7 %, oxalic 52.1 → 47.3 %; the losers are the
three ethanol configurations (corn 18.0 → 27.9, corn stover 27.4 → 38.0,
sugarcane 17.8 → 20.2 %), consistent over seeds (corn stover, seeds 0 / 1
/ 2: 25.8 / 28.0 / 28.2 → 39.2 / 28.9 / 45.7 %), and the same holds out
of family (LOFO: ethanol corn 48 → 111 %, ethanol corn stover 22 → 35 %;
LOPO ethanol 60 → 77 %). Attribution, with the final settings minus the
power gate (same seeds and folds): the ethanol corn stover loss comes from
the surface channel loss (LOCO 27.4 → 51.9 % without the gate; the gate
brings it back to 38.0 %), the ethanol corn loss from the gate (LOCO 19.3
→ 27.9 %, LOFO corn fold 48.8 → 63.3 %, all of the LOFO deterioration:
without the gate LOFO is 29.6 % vs v2 30.3 %), and LOPO ethanol mostly from
the surface loss (59.7 → 73.5 %, gate +3). The gate improves the LOCO
fold mean (10.2 → 9.3 %) and LOPO (52.0 → 46.5 %, seed t −2.7; isobutanol
72.1 → 39.7 %, lactic 55.6 → 49.3 %), and without it 38 % (LOCO), 42 %
(LOPO) and 98 % (LOFO) of the non-exporting held-out samples are again
predicted to export power.

The study's findings replicate on the corrected rows: the surface channel
loss removes most of the decomposition error in every scheme (ch_rms −9.5
LOCO, −7.1 LOPO, −3.9 LOFO points; seed t −12 to −92), the corn
feedstock-sensitivity bias disappears, the power gate removes every
wrong-sign electricity prediction (85–100 % of the held-out samples of the
non-exporting configurations under v2), and the clip exemption fixes HP's
LOPO fold (54.9 → 13.2 % here). MPSP: LOCO and LOPO means better (LOPO −10.7
points, seed t −4.4; 6 of 7 products), the typical in-family configuration
unchanged (geometric mean +0.1, median fold +0.12), LOFO worse (+3.8, seed t
+2.0), carried by ethanol corn and ethanol corn stover. The lactic mask's
new-feedstock gain (sens_fs* 12.7 → 6.9 % in LOFO) comes from the corn fold
alone and is seed-dependent (corn fold, seeds 0 / 1 / 2: 27.8 / 64.8 / 6.8 %
without, 13.3 / 6.1 / 16.1 % with it; the other three folds unchanged or
slightly worse), while it costs LOPO 2.1 points configuration-balanced
(seed t 2.6) and makes lactic's own feedstock sensitivity product-blind.
The TAL mask alone is neutral on every metric (vs. no mask: LOCO +0.0,
LOPO +2.0 / LOFO +0.5 points of MPSP, within seed noise).

HPO criterion with the clip exemption (study data, 16 candidates, seed 0, 2
members; the HP and TAL LOPO folds of the six `clip_inputs = 0` candidates
re-fit with `clip_exempt`, the other folds and candidates unchanged because
the exemption only acts where the ratio lies outside the training range):
every candidate's HP fold becomes sound (10–24 %), the Spearman correlation
of the LOCO and LOPO criteria rises from 0.05 to 0.31, and w = 0.5 still
selects candidate 3 (the v1 choice; LOCO 10.1, LOPO 46.1 %) over the
LOCO-only choice, candidate 7 (LOCO 8.8, LOPO 49.7 % with the exemption; w
≤ 0.25 keeps candidate 7). Over seeds 0–2 that ranking reverses (5d / 5e:
candidate 3 LOPO 54.5 ± 8.0 %, candidate 7 with the exemption 49.5 %), so
after the exemption the LOPO term selects on single-seed noise, at +1.3
points of the LOCO criterion and +21 % run time.

**Decisions (v3 pipeline; library defaults unchanged; each with its CLI
flag, `--v2-model` for (1)–(4)):**

1. *Adopt* the surface channel loss with w_ch = 3 (`CHANNEL_LOSS`,
   `--v2-channel-loss`): the decomposition and price-sensitivity gains are
   large and replicated on both data versions; MPSP is better on average and
   unchanged for the typical in-family configuration; the ethanol losses are
   reported (REPORT "Changes v1 -> v2 -> v3", and the per-configuration
   fidelity table) and stay open.
2. *Adopt* the power gate (`POWER_GATE`, `--no-power-gate`): a hard physical
   constraint verified on the data, which removes a qualitative failure; its
   ethanol-corn cost is the exposure of a compensating error (missed DDGS /
   corn-oil credit), not a new one.
3. *Adopt* the clip exemption (`CLIP_EXEMPT`, `--no-clip-exempt`): removes
   the mechanism of the v2 HP regression, bit-identical elsewhere.
4. *Adopt* the mask of the TAL feedstock-channel targets only
   (`DEFECTIVE_TARGETS = {'feedstock': ('TAL',)}`, `--no-target-mask`): a
   correction of defective targets (D8), neutral within noise. *Do not
   adopt* the lactic mask: lactic's targets are not defective, and its gain
   is one fold and seed-dependent at a LOPO cost.
5. *Do not adopt* the LOPO term of the HPO criterion (`HPO_LOPO_WEIGHT = 0`;
   opt-in `--hpo-lopo-weight`): the exemption removes its motivation, and
   the remaining selection effect is seed noise at +21 % run time and a
   post-selection LOPO estimate.
6. *Do not adopt* the pure-b′ feedstock intensity (opt-in `--fs-b-only`),
   w_ch 1 or 10, or wider clip margins (5c–5e).

`experiments.CODE_VERSION` 4 → 5 (new fit-job semantics: loss space, gate,
clip exemption, target masks, the channel outputs stored for the fidelity
tables). Open: ethanol (its coproduct credits and its turbogenerator /
boiler structure have no product neighbour), isobutanol's electricity (net
zero in the data), the recovery-shaped part of the feedstock intensity
(lactic, a product-specific recovery surface would be needed), succinic's g
under LOPO, and the D8 normalisation / acetate co-feed input.

**Final review of the v3 run (2026-10-18).** Fixed without a refit (serving
path, report wording and report-only statistics; the report regenerated
from the cached fits): the user guides describe the v3 checkpoint; the CV
summary carries the paired fold-level comparison with the baselines (for a
new product kNN − PI-NN +7.9 ± 9.7 pp, PI-NN better on 4 / 7 folds; for a
new feedstock −1.8 ± 10.1 pp, 2 / 4: no advantage over a
nearest-configuration lookup is established); the g-interval coverage per
family (shipped quantile on the LOCO out-of-fold scores: TAL 0.45, lactic
0.62, the other families 1.00); undefined relative price errors masked
below a 2 % share; the user API encodes a missing agent without
neutralization as 'none', accepts boolean strings (pandas-3 'str' columns),
treats a config_id with a composition of its own as a new feedstock, uses a
configuration's baseline titer and prices by name, and warns for
configurations the inputs cannot distinguish (lactic SHF vs. SSCF).
Open, each needs a re-run (a further use of the test split), not adopted
here:

1. *Feedstock intensity of a new composition.* The product-blind feedstock
   head also carries the lactic corn-stover recovery shape (feedstock d′
   43-51, 5b). For the README's wheat straw in the HP_neutral_cornstover
   process it predicts int__fs (a, b, c, d) = (0.247, 1.617, −3.80, 20.06)
   against (0.191, 2.097, 0.920, 1.149) for corn stover: b′ 23 % lower
   although the straw has 5.5 % more sugar, d′ 17 times larger. Proposed:
   keep b′ of the feedstock intensity product-blind (sugar availability)
   and predict its a′, c′, d′ (the recovery shape, F_feedstock / F_sugar is
   constant within a sample) from a product / recovery-side sub-network, or
   use the pure-b′ intensity (`--fs-b-only`) for compositions outside the
   training range; to be judged on LOFO channel fidelity and LOPO, since
   the pure-b′ variant cost LOPO channel error +4 points (5c).
2. *Saccharification regime.* SHF and SSCF lactic configurations differ in
   the data (low pH b′ 0.265 vs. 0.118, baseline MPSP 1.425 vs. 1.308) but
   not in the inputs; a `saccharification` regime field (the lactic module
   records it) needs an adapter field and a data rebuild.
3. *g interval calibration by family.* Pooled g coverage is at nominal, but
   per family it is 0.45-1.00; a group-weighted or family-stratified
   conformal quantile would be an option (it changes the shipped quantiles).
4. *Fixed reagent prices* (succinic_neutral H2SO4, TAL NaOH) recorded as
   `p_acid_price` / `p_base_price` at the next data rebuild (§1, minor items).
5. *Monotonicity priors per product.* dMPSP/dy ≤ 0 and dMPSP/dt ≤ 0 are off
   for every product because isobutanol and oxalic acid violate them in
   the data (§1, "Product-specific cost structure"); enabling them for the
   other products changes training.
6. *Composition features.* For a custom feedstock the energy channels
   follow the composition as a proxy for the feedstock (diluting glucose
   with water lowers the predicted electricity and natural-gas costs, −5.6 %
   MPSP at w_sugar = 0.5, README section 4.2); dry-basis composition
   features plus an explicit water fraction feeding only the channels it
   affects, and a sugar-throughput scale input instead of the wet capacity.
7. *Calibration on known-product folds.* 290 of the 1353 LOCO out-of-fold
   scores that calibrate the shipped intervals come from the two
   new-product folds (isobutanol, oxalic acid); calibrating on the 35
   known-product folds only would match the intervals' stated target (a new
   configuration of a known product; coverage of the shipped quantiles on
   those folds: MPSP 0.90, g 0.87).

### Validation-design checks (final review, round 3; report only)

Added to the report from the saved outputs, without any refit or new
choice (REPORT.md sections 2, 5 and 6; `experiments.clustered_and_twin_checks`):

* **Input twins.** Configurations whose baseline rows encode to the same
  model input (`report_checks.input_twins`, features of `models/final_dev`):
  the two lactic SHF / SSCF pairs (0 of 151 features differ; no other pair
  is closer than 5.09 in the standardized features). One of them crosses
  the development / test boundary: the held-out
  `lactic_neutral_SSCF_cornstover` is input-identical to the development
  configuration `lactic_neutral_SHF_cornstover`, so its OOD result
  measures the unmodelled SHF-to-SSCF difference. The OOD test is also
  reported without it (PI-NN 7.1 % vs. 6.6 % with it; 22 of the 42
  in-domain OOD samples were the twin's, the in-domain fraction without it
  is 0.26), and the LOCO folds whose twin is in training (the low-pH pair)
  are marked; the one lactic fold without a twin in training has 18.1 %.
  The split rule, fixed before modelling, is unchanged; a check of
  identical encodings across a new split belongs before the split is
  frozen.
* **New-product LOCO folds.** `IBO_EtOH_corn` and `oxalic_neutral_sugarcane`
  are their product's only configuration, so their LOCO folds are
  new-product folds; LOCO over the other 35 folds: PI-NN 6.8 %, plain MLP
  8.8 %, kNN 15.5 %.
* **Family clustering.** The 37 LOCO folds cluster in 7 families (15 HP);
  the fold-level p-values treat them as independent. At the family level
  (fold differences averaged per family) the PI-NN is better than the
  plain MLP in 7 of 7 families and than kNN in 6 of 7 (Wilcoxon p 0.016,
  0.031).
* **Data's own constraint violations** next to the methods' (sections 5,
  6), and the measured vs. nominal y and t per configuration (section 1).

## 7. Processed-data interface (contract between `dataset.py` and the model)

`dataset.build(raw_dirs) → ProcessedData` writes/reads three tables (CSV.gz,
in `coeff_predictor/data/processed/`):

**`samples`** — one row per sample that passes quality filters:
* ids: `sample_id`, `config_id`, `biorefinery`, `sample_kind`
  (`baseline|productivity|random`), `seed`;
* specs (from the config record of the sample's own raw file, else the
  ranked one — non-smoke files first; smoke samples superseded by a
  production file are dropped): `feedstock`, `feedstock_type`,
  `fermentation_product`, `final_product`, `neutralization`,
  `neutralizing_agent`, `n_acid_groups`, `pKa1`,
  `final_product_stoich_mass_ratio`, `upgrading`, `separation_variant`,
  `separation_principle`, `aerobic`, `fermentation_mode`, `T_fermentation`,
  `product_water_solubility`, the downstream fields of §5 (patch S1),
  `theoretical_yield`, `inflection_yield` (the
  sample's own; config value in `inflection_yield_config`),
  `y_min, y_max, t_min, t_max`, `baseline_productivity`, `r` (g/L/h),
  `r_mult`;
* descriptors: `fp_*` (fermentation product), `pp_*` (final product) from
  the config record; `fs_*` (feedstock) = median over the sample's ok points
  (harness v2 describes the feedstock stream at every point: sample-specific
  composition and price), else the `sample_start` record (feedstock after
  applying the sample's parameters), else the config record;
* parameters: every `p_*` (canonical) and `bp_*` (biorefinery-specific) value;
  prices actually used: `feedstock_price`, `electricity_price` (median over
  points), `natural_gas_price` (v2: effective price of the purchased gas =
  NG cost / priced NG flow, §2; v1: adapter price);
* targets (primed basis, i.e. ×g): `a, b, c, d`, `g`, `R2`, `n_points`,
  `fit_mape` (MAPE of the 4-coefficient fit itself = irreducible error),
  channel coefficients `ch__<channel>__<k>` for k in a,b,c,d and channels of
  §2, intensity coefficients `int__<fs|el|ng>__<k>` (price-explicit channels:
  `ch__feedstock__k = feedstock_price · int__fs__k`, etc.), `phi`, `F_fixed`
  (median over points of the fixed channel, USD/yr), `closure_rel` (RMS
  closure / RMS MPSP), `channels_ok` (φ identified, closure small, NG pools
  admissible, `F ≥ 0`, measured labor ≤ FOC);
* diagnostics `diag_*`, among them `diag_ng_pool` (`split` = v2 pools,
  else the v1 pool), `diag_ng_price_source`, `diag_ng_mixed_pool`,
  `diag_ng_unpriced_flow_share`, `diag_fixed_source`
  (`labor|labor_foc|foc_intercept`), `diag_labor_cost`,
  `diag_labor_multiplier`, `diag_capex_weights`, `diag_isbl_factor`,
  `diag_foc_alloc`, `diag_F_foc_intercept` (cross-check).

**`points`** — one row per feasible simulated TY point of retained samples:
`sample_id, y, t, MPSP, recovery` (+ optional extra columns).

**`configs`** — one row per configuration with specs and descriptors.

`report.json` additionally holds `rejections_by_biorefinery` and
`per_biorefinery` (`dataset.biorefinery_summary`: kept/rejected, median
R², fit MAPE, closure, φ and its max CV over the productivity samples of a
configuration, `F_fixed`/labor, `channels_ok` fraction, NG price / priced
kg per kg / source, capex weighting).

Harness compatibility: v1 files (no NG pools, labor, OSBL split or
per-point `fs_*`) are processed with the legacy fallbacks above;
`synthetic.generate(..., harness_version=1)` produces such files for the
tests.

Units: coefficient `a'` USD/kg, `b'` USD/kg·(g/g), `c'` USD/kg·(g/L),
`d'` USD/kg·(g/g)(g/L); intensities analogously per unit price
(`int__fs` in kg feedstock/kg final ×(basis), `int__el` in kWh/kg, `int__ng`
in kg NG/kg).
