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
| D4 | ethanol_sugarcane, adapter < 2 | The MPSP was net of the module's RIN credit, recorded as `coproduct_sales` (0.2593 USD per kg of denatured product, checked to be constant before patching). Patched: `x_RIN_credit = coproduct_sales/h`, `MPSP += coproduct_sales/(h·F_final_product)`, `x_MPSP_denatured`/`x_MESP_…` raised by the credit, `coproduct_sales = 0`. |
| D6 | HP, adapter < 2 | `base_cost = x_lime_fermentation · p_base_price`. The TEA charges `lime_fresh`, which is linked to `fermentation_lime`, at the sampled price. v1 read `fermentation_lime`, whose price was frozen. |
| D7 | HP, adapter < 2 | `x_inflection_yield = (1 − p_cell_mass_yield)·theoretical_yield`, `x_above_inflection = y > x_inflection_yield`. Basis checked: the config inflection is 0.755 = 1 − 0.245, and `x_cell_mass_conversion` equals `p_cell_mass_yield` at low y. |
| D2 | ethanol_corn, adapter < 2 | The molar yeast growth destroyed sugar mass. This cannot be patched, so the samples are rejected (`unpatchable_D2`) and the configuration was regenerated. |
| M1 | TAL, succinic low-pH, adapter < 2 | `neutralizing_agent` = 'NaOH' (TAL pH control) or 'Ca(OH)2' (succinic partial pH control). `neutralization` stays False. |
| M3 | TAL, adapter < 2 | `bp_fermentation_fermentation_y_lipolytica_yield` → `p_cell_mass_yield`, in the sample values and in the config's parameter list. |
| M12 | TAL, adapter < 2 | `product_water_solubility` = 9.8 g/L water at 25 °C, interpolated from the TAL module's experimental data. v1 used 5.41 from the fitted model. |
| M23 | HP non-neutral, adapter < 2 | `separation_principle` loses its 'acidulation+' prefix. |

Resumed samples append duplicate rows (M11). `read_raw` already keeps only
the points of each sample's chosen, complete attempt. Points above the
sample's own inflection are dropped (`x_inflection_yield`, `x_above_inflection`;
§7). Defects that need regeneration or are not corrected here (HP acid
regeneration D5, TAL substrate basis D8, lactic D1/D9, isobutanol D10, oxalic
D13/D14, and the minor items not listed) are tracked in the review report.

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
   w.r.t. the productivity input;
4. non-negative feedstock and natural-gas intensities over the domain.

**Collocation points**: unlabeled inputs created by perturbing training inputs
(productivity log-uniform over a range wider than the data, prices over wide
ranges); only physics residuals 1–4 are evaluated there (PINN-style), which
regularizes extrapolation.

Options (`LossConfig`; defaults = the loss above):

* `coef_space='surface'` compares the total and every channel/intensity in
  *surface* space instead of coefficient space: at each simulated point of
  the sample, `(N̂ − N)/N` for the total and `Σ_c ((m̂_c − m_c)/N)²` for the
  channels (price channels as price × intensity), with `N` = target total
  numerator `g·MPSP_fit`. This is well conditioned (the per-sample OLS
  estimates of a', b', c', d' are strongly correlated, §6), but in the
  development CV it made in-family configurations worse (§6), so it is off;
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
* regime: neutralization flag and agent, aerobic flag, fermentation mode,
  temperature, productivity (log r and 1/r), inflection/theoretical yield ratio,
  separation principle;
* contextual (optional, baseline-imputed with missing-indicators): electricity
  price, natural gas price, IRR, income tax, operating days, boiler/
  turbogenerator efficiency, etc.

Descriptors of known feedstocks/products are stored in a registry so users can
specify them by name (composition and baseline capacity per feedstock; yields,
acid properties per product; per product pair the stoichiometric ratio, the
fermentation regime where it is unique in training — otherwise the user is
asked to specify it — and the baseline value of every contextual `p_*`
parameter, or "not applicable" which keeps its missing indicator at 1); new
chemicals can be described with thermosteam or by supplying descriptor values.
Omitted categoricals of a query are encoded as the unknown level (trained by
categorical dropout), and every imputed / defaulted field is reported per row.
Input units and bases are listed in `features.FIELD_UNITS` (feedstock prices
per kg WET feedstock). The applicability domain (kNN + range) uses encoder
features and productivity only: price-explicit channels are exactly linear in
price, so prices outside the training range are reported separately.

## 6. Training & validation protocol

* Deep ensemble (M members, different seeds/bootstraps) → mean ± epistemic
  spread; intervals calibrated by split-conformal scaling on a held-out
  calibration set of whole groups (configurations by default; in CV the
  inner calibration / early-stopping splits group like the outer scheme:
  feedstock, fermentation product or configuration). One score per
  calibration sample (MPSP: max over its points); if there are too few
  samples for the level the quantile is +inf (no interval is reported);
  `ens.conformal` stores n, the guaranteed level, the grouping actually used
  and `coverage_guaranteed` (needs ≥ 10 calibration groups).
* Splits are **grouped** (never split TY points of one sample; MC samples of a
  configuration are correlated):
  - random sample-level K-fold (interpolation);
  - leave-one-configuration-out;
  - leave-one-feedstock-out;
  - leave-one-product-out (hardest; generalization to new products);
  - productivity extrapolation (train on the central productivity range).
* Hyperparameters chosen on inner validation only; a final test set
  (untouched until the end) reported once. It is **frozen** in
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
  (full mode), selection by the fold mean of the **configuration-balanced**
  held-out MPSP MAPE (`mpsp_mape_cfg`: mean over held-out configurations of
  their point MAPE; `compute_metrics` also reports its geometric mean
  `mpsp_mape_cfg_gmean`). Reasons (development diagnostics below): the
  point-weighted fold MAPE is dominated by the fold's largest configuration,
  and single-member results differ between seeds by up to ~50 % (relative).
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
  fit), physics-constraint violation rates, interval coverage.

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
  `product_water_solubility`, `theoretical_yield`, `inflection_yield` (the
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
