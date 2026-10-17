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
| S1 | every known configuration, `downstream_spec_version` < 1 (all recorded data) | Not a defect fix: adds the structured final-product / downstream-process specification of `downstream.TABLE` (§5) to the config record (the table wins; a differing recorded `coproducts` / `neutralization_degree` is kept as `<field>_recorded`). Versioned by `downstream_spec_version`; `downstream.apply_to_frame` does the same for processed tables built before S1. |

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
| neutralization | neutralization flag / agent / degree, `n_acid_groups`, `pKa1`, `phys_{base,acid,gypsum}_per_product`, `p_{base,acid}_price`, `ds_sep__acidulation`, `ds_sep__gypsum_filtration` |
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
sub-network of its groups (`CHANNEL_INPUT_PRESETS['physics']`):

| output | groups |
|---|---|
| feedstock intensity | feedstock, pretreatment |
| capex_upstream | feedstock, pretreatment, scale, financial |
| fixed, closure | feedstock, scale, financial |
| capex_fermentation | fermentation, productivity, neutralization, feedstock, scale, financial |
| neutralization | neutralization |
| capex_downstream | product, final_product, upgrading, separation, water_route, neutralization, scale, financial |
| utilities_other | product, final_product, upgrading, separation, water_route, fermentation, utilities, scale |
| natural-gas intensity | product, final_product, upgrading, separation, water_route, feedstock, utilities, scale |
| materials_other | feedstock, pretreatment, fermentation, product, final_product, upgrading, separation, water_route, scale |
| capex_wastewater | feedstock, fermentation, water_route, coproducts, scale, financial |
| capex_utilities | feedstock, water_route, upgrading, coproducts, utilities, scale, financial |
| electricity intensity | feedstock, fermentation, productivity, water_route, upgrading, coproducts, utilities, scale |
| coproducts | feedstock, fermentation, coproducts, scale |
| g | product, final_product, upgrading, separation, water_route, neutralization |

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

Known caveat: the succinic_neutral flowsheet acidulates (R401 + S405), but
its recorded `acid_cost` is 0 (module behaviour; not corrected).

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
  fit), physics-constraint violation rates, interval coverage. Additions
  (2026-10-17):
  - MAPE split by sample kind (`*_spec`, `*_random`, also configuration-
    balanced): *spec-determined* samples (`sample_kind` baseline /
    productivity: every biorefinery-specific `bp_*` parameter at baseline,
    so the inputs determine the target) vs *random* samples (unobserved
    `bp_*` draws: irreducible noise for any model);
  - `evaluate.noise_floor`: within each configuration (≥ 10 random
    samples), a 5-fold ridge on the observable inputs (`p_*`, prices, log r,
    1/r) predicts a held-out random sample's coefficient vector; its surface
    MAPE estimates the error the unobserved `bp_*` cause for a model that
    knows the configuration exactly (also reported: + `bp_*` inputs, and the
    configuration mean);
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
