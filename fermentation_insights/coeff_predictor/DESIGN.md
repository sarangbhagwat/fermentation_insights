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
labor-type fixed costs (intercept of FOC regressed on `I`; constant within a
sample), `FOC_i − F` the capital-related O&M (maintenance, insurance), and
`φ_K` the capital charge per USD of FCI (annualized capital recovery at the
IRR incl. taxes and working capital) from the no-intercept least squares of
`R − ΣE − FOC` on `I` (financial parameters are constant within a sample).
`φ = φ_K + dFOC/dI` is reported. `ε_i` is the closure residual (kept as its
own channel so that channels sum exactly to MPSP); it absorbs DCF weights
≠ 1 (e.g. startup fractions) instead of letting a free intercept turn them
into a negative or productivity-dependent `F`. Without FOC the legacy free
fit `R − ΣE = φ·I + F` is used. Samples with a negative `F` are flagged
(`channels_ok = False`).

Channels (contribution to MPSP = channel USD/yr ÷ (F_final·h)):

| channel | content | price-explicit? |
|---|---|---|
| `feedstock` | wet feedstock purchase | yes: `p_fs × feedstock intensity` |
| `electricity` | net purchased (+) / exported (−) electricity | yes: `p_el × net electricity intensity` |
| `natural_gas` | purchased natural gas | yes: `p_ng × NG intensity` |
| `neutralization` | base for pH control + acid for acidulation | no |
| `materials_other` | all other purchased materials (enzymes, nutrients, catalysts, solvents, boiler chemicals …) | no |
| `utilities_other` | non-electric utility agents (net, after boiler offsets) | no |
| `capex_upstream` | capital charge (φ_K·FCI + FOC − F) × weight share of feedstock + conversion areas | no |
| `capex_fermentation` | capital charge × weight share of fermentation | no |
| `capex_downstream` | capital charge × weight share of separation + upgrading | no |
| `capex_wastewater` | capital charge × weight share of wastewater | no |
| `capex_utilities` | capital charge × weight share of CHP + facilities (+ uncovered installed cost) | no |
| `fixed` | F (labor-type fixed costs) | no |
| `coproducts` | − coproduct sales | no |
| `closure` | ε (should be small) | no |

(Area weights are the areas' TDC contributions `I_isbl·(1 + f) + I_osbl`
when the harness records each area's OSBL installed cost
(`<area>__installed_cost_osbl`; `f` = warehouse + site development +
additional piping, recorded or estimated from FCI), else installed costs.
Natural gas: when the harness records the NG cost per pool
(`ng_cost_utility`, `ng_cost_material`, at the prices BioSTEAM charges, e.g.
boiler 'Fuel' vs. dryer 'Natural gas'), each part is subtracted from its own
pool, `natural_gas_price` is the sample's cost-weighted price and `int__ng`
the NG flow equivalent at that price (exact `p × q`). Legacy data (one
recorded price) use a pool decided once per configuration (CHP-area
material-cost evidence, else the sign of the remainders); samples whose
chosen pool leaves a negative remainder — mispriced or split NG — get
`channels_ok = False`. Known limitation of existing data: TAL/succinic/corn
ethanol dryer NG is priced differently from the recorded price until the
adapters record per-pool costs.)

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
* Hyperparameters chosen on inner validation only; a final test set of whole
  configurations (untouched until the end; requires ≥ 2 configurations)
  reported once.
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
* descriptors: `fp_*` (fermentation product), `pp_*` (final product), `fs_*`
  (feedstock) from the config record; `fs_*` in the `sample_start` record
  (feedstock described after applying the sample's parameters) override;
* parameters: every `p_*` (canonical) and `bp_*` (biorefinery-specific) value;
  prices actually used: `feedstock_price`, `electricity_price`,
  `natural_gas_price` (median over points);
* targets (primed basis, i.e. ×g): `a, b, c, d`, `g`, `R2`, `n_points`,
  `fit_mape` (MAPE of the 4-coefficient fit itself = irreducible error),
  channel coefficients `ch__<channel>__<k>` for k in a,b,c,d and channels of
  §2, intensity coefficients `int__<fs|el|ng>__<k>` (price-explicit channels:
  `ch__feedstock__k = feedstock_price · int__fs__k`, etc.), `phi`, `F_fixed`,
  `closure_rel` (RMS closure / RMS MPSP), `channels_ok` (φ identified,
  closure small, NG pool admissible, `F ≥ 0`).

**`points`** — one row per feasible simulated TY point of retained samples:
`sample_id, y, t, MPSP, recovery` (+ optional extra columns).

**`configs`** — one row per configuration with specs and descriptors.

Units: coefficient `a'` USD/kg, `b'` USD/kg·(g/g), `c'` USD/kg·(g/L),
`d'` USD/kg·(g/g)(g/L); intensities analogously per unit price
(`int__fs` in kg feedstock/kg final ×(basis), `int__el` in kWh/kg, `int__ng`
in kg NG/kg).
