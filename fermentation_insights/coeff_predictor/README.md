# coeff_predictor: physics-informed prediction of the unifying-equation coefficients

`coeff_predictor` predicts the coefficients of the unifying equation of
Bhagwat et al. ([*Nat. Commun.* 2026](https://doi.org/10.1038/s41467-026-75285-1))
for a fermentation-based biorefinery from a short specification of its
**feedstock**, **fermentation product**, **final product** and
**fermentation regime**, without building or running a process model. It is
a deep ensemble trained on BioSTEAM simulations of 7 Bioindustrial-Park
biorefineries (41 configurations; uncertain parameters and productivities
sampled, titer-yield grids simulated). The simulations record intermediate
process- and system-level data, and the model is built around that
structure. The design is documented in [DESIGN.md](DESIGN.md); § below
refers to its sections.

## 1. What it predicts, and why

    MPSP(y, t) = (a' + b'/y + c'/t + d'/(y·t)) / g

MPSP is the minimum product selling price [USD/kg final product], `y` the
yield [g fermentation product / g glucose-eq. fermentable sugar] and `t` the
titer [g/L]. The coefficients are in the *primed* basis of the publication
(a' = g·A, where A, B, C, D are the MPSP-level coefficients). Since g·MPSP is
the cost per kg of **fermentation product**:

| coefficient | unit | meaning |
|---|---|---|
| a′ | USD/kg | titer- and yield-independent cost per kg fermentation product (e.g. costs proportional to the product flow, such as purification and upgrading) |
| b′ | USD/kg·(g/g) | yield-scaled cost (∝ 1/y): the sugar, i.e. feedstock, needed per kg product; feedstock purchase and, at the modules' fixed feedstock capacity, the costs sized by the feedstock flow |
| c′ | USD/kg·(g/L) | titer-scaled cost (∝ 1/t): the broth volume and water handled per kg product; fermenters, water removal (evaporation, distillation), wastewater, utilities |
| d′ | USD/kg·(g/g)(g/L) | interaction of yield and titer (∝ 1/(y·t)) |
| g | kg final / kg fermentation product | recovery: final product sold per kg fermentation product leaving fermentation, incl. the stoichiometric mass change of upgrading (≤ the stoichiometric ratio `g_max`) |

With the coefficients, the MPSP over the whole titer-yield landscape, the
titer and yield targets for a price, and the relative impact of yield vs.
titer improvements follow in closed form (section 4.4). The training data are
per-sample least-squares fits of exactly this equation to simulated MPSPs
(yields up to the inflection yield, as in the publication).

## 2. Physics-informed design in brief

* **Cost channels (§2).** The revenue required at the MPSP is decomposed
  exactly into annual cost channels measured in the simulations (feedstock,
  electricity, natural gas, neutralization, other materials and utilities,
  capital by process area, fixed costs, coproduct credits, closure). Every
  channel is fitted with the same basis `[1, 1/y, 1/t, 1/(y·t)]`, so the
  channel coefficients sum exactly to a′, b′, c′, d′. The network predicts
  the channels (intermediate supervision), and the totals are their sum (a
  hard identity).
* **Price-explicit intensities (§2, §3).** Feedstock, electricity and natural
  gas enter as price × predicted intensity (kg feedstock, kWh, kg gas per kg
  product). The prices are not network inputs, so predictions are exactly
  linear in them. The feedstock intensity is structured as
  `b'_fs = (1/η)/w_sugar` (η = sugar reaching fermentation); in the pipeline
  1/η = (1/s_avail)·exp(Δz), where the composition's sugar availability
  s_avail (free sugars and starch fully, structural polysaccharides by the
  saccharification yield) is a physics prior and the network predicts only
  the residual Δz.
* **Hard productivity structure (§3).** The fermentation capital and the
  electricity intensity are θ0(x) + θ1(x)·r_ref/r with θ1 ≥ 0, and no other
  output sees the productivity, so MPSP never increases with r.
* **Equation layer and hard bounds (§3, §4).** The equation turns the totals
  into MPSP at the simulated points, and the MPSP error is part of the loss.
  The recovery is bounded, `g = g_max·sigmoid(·)`.
* **Soft constraints at collocation points (§4).** MPSP non-increasing in
  yield, titer and productivity, positive MPSP and non-negative intensities are
  penalized at labelled points and at unlabelled, perturbed inputs (only
  constraints the training data satisfy are enabled).
* **Channel-wise input masking (§3).** Each output has its own sub-network
  that sees only the features of the process areas that can affect it. For
  example, the feedstock intensity and upstream capital see no product
  information, and the recovery sees no feedstock information.
* **Downstream specification and water-removal gate (§3, §5).** The final
  product and its recovery route are described with controlled
  vocabularies: purity, phase, water-removal route, separation and upgrading
  operations, coproducts, neutralization degree, acid regeneration. When the
  broth itself is the product (route `none`), the gate sets the
  titer-dependent (c′, d′) parts of the broth-water channels to zero.
* **Extrapolation safeguards (§3).** The feedstock-composition inputs are
  clipped to their training range (a new feedstock's minor mass fractions
  cannot saturate the network), and the neutralization cost is exactly zero
  for products with nothing to neutralize (no acid group, no neutralizing
  agent).
* **Uncertainty (§6).** Deep ensemble, split-conformal intervals calibrated
  on the out-of-fold scores of the development leave-one-configuration-out
  CV of the same training procedure (so they target new configurations), and
  a kNN applicability domain.

Defaults: the library defaults are `FeatureSpec(downstream=False)` and
`ModelConfig(water_gate=False, channel_inputs=None, r_structure=None,
fs_prior=False, clip_groups=None, neutralization_gate=False)`. The experiment
pipeline (`experiments.py`, v2), and so the release checkpoint, uses the
downstream features, the gate, masking with the preset `'physics_v2'` (the
recovery does not see the reagent prices), the productivity structure and the
extrapolation safeguards (`--no-downstream`, `--no-water-gate`,
`--no-masking`, `--v1-g-inputs`, `--no-r-structure`,
`--no-extrapolation-safeguards` switch them off; `--v1-protocol` restores the
v1 training protocol); DESIGN.md §3 and §6 have the development-set evidence
for these choices.

## 3. Installation

The predictor needs Python ≥ 3.12, PyTorch (CPU is enough), NumPy, pandas,
SciPy and scikit-learn. matplotlib is needed for the experiment figures and
pytest for the tests. The package is imported from the repository, so add
the repository root to `PYTHONPATH`:

```bash
conda create -n fi_ml python=3.12 && conda activate fi_ml
pip install torch numpy pandas scipy scikit-learn matplotlib pytest
pip install thermosteam==0.54.2      # optional: descriptors of new chemicals
export PYTHONPATH=/path/to/fermentation_insights
python -m pytest fermentation_insights/coeff_predictor/tests   # ~1 min, synthetic data only
```

Tested with torch 2.6.0, numpy 2.5.4, pandas 3.0.6, scipy 1.18.1,
scikit-learn 1.9.1, matplotlib 3.11.2 and thermosteam 0.54.2. Neither
BioSTEAM nor thermosteam is needed for predictions. The descriptors of the
known chemicals are stored in `registry/chemicals.json`, and thermosteam is
only used to describe a new chemical (`descriptors.chemical_descriptors`).

A checkpoint is a directory written by `train.Ensemble.save`. The release
checkpoint is `coeff_predictor/models/release`, written by `experiments.py`
(section 8).

## 4. Quickstart

The examples below are implemented in
[`examples/quickstart.py`](examples/quickstart.py), which runs them end to end
against a checkpoint:

```bash
python fermentation_insights/coeff_predictor/examples/quickstart.py [models/release] [--no-thermosteam]
```

### 4.1 A known configuration, by name

```python
from fermentation_insights.coeff_predictor.predict import CoefficientPredictor
P = CoefficientPredictor.load('fermentation_insights/coeff_predictor/models/release')

known = dict(
    config_id='HP_neutral_cornstover',          # fills the downstream fields (downstream.TABLE)
    feedstock='cornstover', fermentation_product='HP', final_product='AcrylicAcid',
    neutralization=True, neutralizing_agent='Ca(OH)2',   # regime fields that differ between HP designs
    r=0.28,                                     # productivity [g/L/h]
    feedstock_price=0.0517,                     # USD per kg WET corn stover
)
p = P.predict(known)                    # one row per spec (dict, list of dicts or DataFrame)
p[['a', 'b', 'c', 'd', 'g']]
P.unprimed(known)                       # A, B, C, D = a'/g, ...
```

Names fill everything the training data define for them: the feedstock
composition and baseline capacity, the yields and acid properties of the
product, the stoichiometric ratio, the regime fields that are unique for the
product pair in training, and the baseline contextual parameters. The output
lists what was filled (`filled_from_registry`), defaulted
(`contextual_defaults`) or imputed (`imputed`), and UserWarnings repeat it. A
regime field that varies between the configurations of a product pair should
be given; otherwise it is encoded as unknown, and a warning names it. The values of the known configurations are in
their adapters' `configurations()` (`data_generation/adapters/`), and
`P.known_downstream(config_id)` returns their downstream fields.

### 4.2 A new feedstock (custom composition and price)

A new feedstock needs no name. Give its type, its complete wet-basis
composition (`fs_w_*`; components left out are imputed), its price and,
optionally, its capacity:

```python
from fermentation_insights.coeff_predictor.features import W_SUGAR   # 'fs_w_fermentable_sugar_eq'
dry = dict(glucan=0.354, xylan=0.197, arabinan=0.024, galactan=0.008, mannan=0.003,
           lignin=0.170, ash=0.065, protein=0.030, extractives=0.120)       # wheat straw, dry basis
moisture = 0.15
comp = {d['name']: 0. for d in P.spec.numeric if d['name'].startswith('fs_w_')}
comp.update({f'fs_w_{k}': v * (1 - moisture) for k, v in dry.items()}, fs_w_water=moisture)
comp['fs_w_other'] = 1 - sum(v for k, v in comp.items() if k not in ('fs_w_other', W_SUGAR))
gain = dict(glucan=180.156/162.141, galactan=180.156/162.141, mannan=180.156/162.141,
            xylan=150.13/132.115, arabinan=150.13/132.115)             # hydrolysis mass gain
comp[W_SUGAR] = sum(comp[f'fs_w_{k}'] * f for k, f in gain.items())  # glucose-eq. sugar, kg/kg wet
straw = {**known, 'feedstock_type': 'lignocellulose', 'feedstock_price': 0.0595,
         'p_feedstock_capacity': 1.0e5, **comp}
straw.pop('feedstock')
P.predict(straw)
```

Other gains, as in `data_generation/harness.describe_feedstock`: 1 for
monomeric sugars, 2·180.16/342.30 for sucrose and 180.16/162.14 for starch.

### 4.3 A new fermentation product

A new product needs its yields, stoichiometric ratio, acid properties,
regime, chemical descriptors and, for checkpoints with the downstream
specification, the downstream fields:

```python
from fermentation_insights.coeff_predictor.descriptors import chemical_descriptors
desc = chemical_descriptors('Itaconic acid')       # needs thermosteam; or a manual dict:
# desc = dict(MW=130.099, Tb=601.0, Tm=444.15, Hvap_298=1.1407e5, Psat_ferm=1.98e-3,
#             LHV=1.8517e6, n_C=5, n_H=6, n_O=4, n_N=0, gamma_inf=7.78, log10_alpha_inf=-5.61)
itaconic = dict(
    feedstock='glucose', feedstock_price=0.28, r=0.6,
    fermentation_product='ItaconicAcid', final_product='ItaconicAcid',
    ferm_product_chemical='ItaconicAcid', final_product_chemical='ItaconicAcid',
    theoretical_yield=0.722, inflection_yield=0.58, final_product_stoich_mass_ratio=1.0, upgrading=False,
    n_acid_groups=2, pKa1=3.85, product_water_solubility=83.,
    neutralization=False, neutralizing_agent='none', aerobic=True, fermentation_mode='fed-batch',
    T_fermentation=310.15,
    final_product_mass_fraction=0.99, final_product_is_broth=False, final_product_phase='solid',
    water_removal_route='crystallize',
    separation_ops='centrifugation+evaporation+crystallization+filtration+drying',
    upgrading_ops='none', coproducts='electricity', neutralization_degree=0., acid_regeneration='none',
)
p_new = P.predict(itaconic, descriptors={'ItaconicAcid': desc})
```

The descriptor keys are those of `registry/chemicals.json`. `dor` and
`Ymax_*` are completed from the atom counts and MW.
`features.register_chemical(ID, desc)` (or `features.describe_new_chemical`)
registers a chemical for the session instead of passing `descriptors=`. A
value outside a vocabulary raises ValueError (`P.downstream_vocabulary`), and
`P.missing_fields(spec)` lists the required fields that are still missing.

### 4.4 MPSP surface and relative impact of yield vs. titer

```python
import numpy as np
r = P.predict(known).iloc[0]
Y, T = np.meshgrid(np.linspace(0.2, 0.755, 5), [10, 20, 50, 100, 200])   # y up to the inflection yield
MPSP = (r.a + r.b / Y + r.c / T + r.d / (Y * T)) / r.g                   # USD/kg final product
# MPSP decrease from +0.01 g/g yield relative to +1 g/L titer (publication, MPSP_YT_fit_only.py);
# > 1: the yield improvement is worth more. g cancels.
RI = (0.01 * (r.b + r.d / T) / Y**2) / ((r.c + r.d / Y) / T**2)
m = P.mpsp(known, y=Y.ravel(), t=T.ravel())     # ensemble: MPSP, MPSP_std, MPSP_lo, MPSP_hi per (y, t)
```

`P.mpsp` averages the members' MPSPs (not the equation at the mean
coefficients) and adds the calibrated intervals. It raises ValueError for
y ≤ 0 or t ≤ 0, returns NaN above the theoretical yield, and warns above the
inflection yield or outside the training titer range.

### 4.5 Uncertainty and applicability domain

* `<k>_std` is the ensemble spread, and `<k>_lo`/`<k>_hi` the split-conformal
  interval (level `1 − P.ensemble.conformal['alpha']`, default 90 %) for
  k = a, b, c, d, g. `mpsp()` returns `MPSP_lo`/`MPSP_hi` (the MPSP interval
  is simultaneous over a sample's titer-yield points). Intervals are omitted
  when the calibration set was too small for the level.
* The release (and `final_dev`) checkpoints have no calibration hold-out:
  their quantiles come from the out-of-fold scores of the development
  leave-one-configuration-out CV of the same training procedure and ensemble
  size (`P.ensemble.conformal['source']`, `['n']`, `['n_calib_groups']`).
  `P.ensemble.conformal['coverage_guaranteed']` is therefore False (the
  scores come from the CV fold ensembles, not from the released one), and
  `P.ensemble.conformal['crossfit']` holds a cross-fitted estimate of the
  coverage (per key: mean, pooled, worst and 10 % quantile over the
  configurations). The coverage is marginal over the samples of new
  configurations (configurations weighted by their numbers of samples), not
  guaranteed per configuration: it varies widely between configurations, so a new configuration of a family whose
  errors are larger than usual can be under-covered, while for a query close
  to the training configurations the intervals are conservative. REPORT.md
  (section 6) compares the estimate with the realized test coverage; no test
  sample was used to calibrate. `P.ensemble.meta['conformal_random_cv']` holds,
  for information only, the quantiles calibrated on the random-CV scores
  (in-distribution queries).
* `ad_ratio` is the mean distance to the 5 nearest training rows, divided by
  the 99 % quantile of the training rows' own distances, in the space of the
  encoder features and productivity. `in_domain` = `ad_ratio ≤ 1` and no
  feature outside its training range (`outside_range` lists those). Outside
  the domain, an `ApplicabilityWarning` is issued; treat the prediction as an
  extrapolation. Prices are not part of the domain (the price channels are
  exactly linear), and prices outside the training range are listed in
  `price_outside_range`.

### 4.6 Channel breakdown

`ch__<channel>__<k>` (k = a, b, c, d; primed basis) are the predicted cost
channels. They sum to the totals, and `ch__<channel>__a + …/y + …/t +
…/(y·t)` divided by g is that channel's share of the MPSP. For the
price-explicit channels, `int__<fs|el|ng>__<k>` are the intensities and
`ch__feedstock__k = feedstock_price · int__fs__k`, so the effect of a
feedstock, electricity or natural-gas price change is exact. The channel
names are in `P.ensemble.channel_names()`.

```python
names = P.ensemble.channel_names()
breakdown = {c: [r[f'ch__{c}__{k}'] for k in 'abcd'] for c in names}
```

### 4.7 A known recovery

If the recovery of the downstream process is known, pass it as `recovery=`
(a scalar, or one value per row; None/NaN = predicted) or as a `recovery`
field of the spec (the argument takes precedence):

```python
P.predict(itaconic, descriptors={'ItaconicAcid': desc}, recovery=0.85)   # g = 0.85, g_source = 'user'
P.mpsp(itaconic, y=0.5, t=80., descriptors={'ItaconicAcid': desc}, recovery=0.85)
```

The recovery replaces the predicted g in `g`, in the MPSP and in
`unprimed()`. The predicted g is kept as `g_predicted`, and a′, b′, c′, d′,
their channels and their intervals are unchanged. It must satisfy
0 < g ≤ g_max = `final_product_stoich_mass_ratio`, else ValueError.

*Approximation.* g·MPSP is the cost per kg of fermentation product. At a
fixed process design this cost does not depend on how much of the product is
recovered downstream. To first order, a different recovery only changes the
amount of final product over which the same costs are spread:
MPSP = N(y, t)/g. Three effects are not captured: costs that scale with the
recovered product mass (purification, drying, upgrading reagents, storage);
the design changes, and their costs, that a different recovery would need;
and coproduct credits from unrecovered product. The MPSP intervals then
reflect the coefficient uncertainty only (their calibration included the
error of the predicted g).

### Output columns of `predict()`

`a, b, c, d, g` (ensemble means; `g` is the user recovery where one was
given), `g_predicted`, `g_source`, `<k>_std`, `<k>_lo`, `<k>_hi`,
`ch__<channel>__<k>`, `int__<fs|el|ng>__<k>`, `ad_ratio`, `in_domain`,
`outside_range`, `price_outside_range`, `imputed`, `contextual_defaults` and
`filled_from_registry`. `mpsp()` returns `row, y, t, MPSP, MPSP_std,
MPSP_lo, MPSP_hi, g_source`.

## 5. Field reference

Generated from `features.FIELD_UNITS` (which merges `downstream.FIELDS`;
`P.field_units` at run time). Multi-valued fields are '+'-joined, and
`'none'` is an explicit empty value. Categorical values that were not seen
in training are encoded as the unknown level. Every omitted, defaulted or
imputed field is reported per row.

| field | unit / allowed values | if omitted |
|---|---|---|
| **Productivity, prices, recovery** | | |
| `r` | g/L/h (fermentation productivity, > 0) | **required** |
| `feedstock_price` | USD/kg wet (as-received) feedstock | **required** (aliases `p_feedstock_price`, `fs_price`) |
| `electricity_price` | USD/kWh | training baseline (`imputed`) |
| `natural_gas_price` | USD/kg natural gas | training baseline (`imputed`) |
| `recovery` | kg final product / kg fermentation product (optional known recovery g, replaces the predicted g; 0 < g <= final_product_stoich_mass_ratio) | predicted g (`recovery=` argument takes precedence) |
| **Names and chemical IDs** | | |
| `feedstock` | name; a known one fills composition, type and baseline capacity | no registry fill |
| `fermentation_product` | name; a known one fills yields and acid properties | no registry fill |
| `final_product` | name; a known product pair fills the stoichiometric ratio, the regime where unique and the contextual baselines | no registry fill |
| `config_id` | name of a known configuration (downstream.TABLE); fills the downstream fields | downstream fields from the product pair (if unique) or imputed |
| `ferm_product_chemical` | chemical ID of the fermentation-product descriptors (registry/chemicals.json, register_chemical or descriptors=) | product name (registry), else descriptors imputed |
| `final_product_chemical` | chemical ID of the final-product descriptors | product pair (registry), else descriptors imputed |
| **Feedstock** | | |
| `feedstock_type` | pure_sugar, starch, sugar_juice or lignocellulose | feedstock name, else unknown level |
| `fs_w_*` | kg/kg wet feedstock (mass fractions, wet basis) | feedstock name, else training median (`imputed`) |
| `fs_w_fermentable_sugar_eq` | kg glucose-eq. fermentable sugar / kg wet feedstock, in (0, 1.12] (> 1 only for anhydro sugars, e.g. pure starch 1.11) | **required** unless the feedstock is known by name |
| `p_feedstock_capacity` | kg/h wet feedstock | feedstock name (baseline), else training baseline (`contextual_defaults`) |
| **Fermentation product and regime** | | |
| `theoretical_yield` | g fermentation product / g glucose-eq. sugar | **required** unless the fermentation product is known by name |
| `inflection_yield` | g fermentation product / g glucose-eq. sugar | **required** unless the fermentation product is known by name |
| `final_product_stoich_mass_ratio` | kg final product / kg fermentation product (stoichiometric, > 0) | **required** unless the product pair is known by name |
| `n_acid_groups` | carboxylic acid groups per molecule of fermentation product (0 if not an acid) | fermentation-product name, else training median |
| `pKa1` | first pKa of the fermentation product (NaN if not an acid) | fermentation-product name, else training median |
| `product_water_solubility` | g/L water at 25 C (1000 if miscible) | fermentation-product name, else training median |
| `upgrading` | bool: chemical conversion of the fermentation product after fermentation | product pair, else training median |
| `neutralization` | bool: pH control by base addition in fermentation | product pair (if unique), else training median |
| `neutralizing_agent` | base dosed in fermentation, e.g. 'Ca(OH)2', 'NaOH'; 'none' if none | product pair (if unique), else unknown level |
| `aerobic` | bool | product pair (if unique), else training median |
| `fermentation_mode` | batch, fed-batch or continuous | product pair (if unique), else unknown level |
| `T_fermentation` | K | product pair (if unique in training), else training median |
| `separation_principle` | free text; used only by checkpoints trained without the downstream fields | product pair (if unique), else unknown level |
| **Final product and downstream process (checkpoints with `FeatureSpec(downstream=True)`)** | | |
| `final_product_mass_fraction` | kg product / kg final-product stream at baseline, in (0, 1]; None for a broth product (see final_product_is_broth) | `config_id` or product pair (if unique), else imputed (none for a broth product) |
| `final_product_is_broth` | bool: the final product is the fermentation broth (w_f ~ t/1000, no water removal) | `config_id` or product pair (if unique), else imputed |
| `final_product_phase` | one of liquid, solid, aqueous_solution | `config_id` or product pair (if unique), else imputed |
| `water_removal_route` | one of evaporate_water, crystallize, distill_product, none | `config_id` or product pair (if unique), else imputed |
| `separation_ops` | '+'-joined subset of: centrifugation, filtration, acidulation, gypsum_filtration, adsorption, ion_exchange, membrane, flash, evaporation, distillation, LLE, decantation, molecular_sieves, crystallization, drying, esterification, hydrolysis | `config_id` or product pair (if unique), else imputed |
| `upgrading_ops` | '+'-joined subset of: dehydration, hydrogenation, ring_opening_hydrolysis, esterification, hydrolysis, decarboxylation, oxidation, amination, polymerization or 'none' | `config_id` or product pair (if unique), else imputed |
| `coproducts` | '+'-joined subset of: electricity, DDGS, corn_oil, ethanol, lignin, gypsum, CO2, fertilizer or 'none' | `config_id` or product pair (if unique), else imputed |
| `neutralization_degree` | mol base-eq. / mol product acid groups converted to salt in fermentation, in [0, 1] | `config_id` or product pair (if unique), else imputed |
| `acid_regeneration` | acid regenerating the free acid downstream: 'H2SO4' or 'none' | `config_id` or product pair (if unique), else imputed |
| **Contextual parameters** | | |
| `p_operating_days` | d/yr | product-pair baseline, else training baseline (`contextual_defaults`) |
| `p_IRR` | fraction (internal rate of return) | product-pair baseline, else training baseline (`contextual_defaults`) |
| `p_income_tax` | fraction | product-pair baseline, else training baseline (`contextual_defaults`) |
| `p_boiler_efficiency` | fraction | product-pair baseline, else training baseline (`contextual_defaults`) |
| `p_turbogenerator_efficiency` | fraction | product-pair baseline, else training baseline (`contextual_defaults`) |
| `p_cell_mass_yield` | g cell mass / g sugar (biorefinery convention) | product-pair baseline, else training baseline (`contextual_defaults`) |
| `p_base_price` | USD/kg base | product-pair baseline, else training baseline (`contextual_defaults`) |
| `p_acid_price` | USD/kg acid | product-pair baseline, else training baseline (`contextual_defaults`) |
| `p_enzyme_loading` | biorefinery-specific units | product-pair baseline, else training baseline (`contextual_defaults`) |
| `p_saccharification_yield` | fraction | product-pair baseline, else training baseline (`contextual_defaults`) |
| **Query and output** | | |
| `y` | g fermentation product / g glucose-eq. sugar | argument of `mpsp()` |
| `t` | g/L | argument of `mpsp()` |
| `MPSP` | USD/kg final product | output |

Names of the planned configurations: the feedstocks `glucose`, `corn`,
`sugarcane` and `cornstover`, and the product pairs (`fermentation_product`
`|` `final_product`) `HP|AcrylicAcid`, `TAL|TAL`, `TAL|KSA` (potassium
sorbate), `SuccinicAcid|SuccinicAcid`, `LacticAcid|LacticAcid`,
`Ethanol|Ethanol`, `OxalicAcid|OxalicAcid` and `Isobutanol|Isobutanol`. A
checkpoint knows the names in its training data (`P.spec.spec_registry`),
and an unrecognized name raises a warning with close matches.

## 6. Validation

`experiments.py` writes every validation number, table and figure to
**[`results/REPORT.md`](results/REPORT.md)** (with `results/tables/*.csv` and
`results/figures/*.png`; written by the full run). A `--quick` debug run
writes to `results/quick/`; its numbers are not results. The protocol (§6):

* **Frozen final test set.** [`results/split.json`](results/split.json) was
  committed before any modelling and is derived only from the planned
  configuration list and a seed. It holds out four whole configurations
  (`HP_corn`, `TAL_sugarcane`, `succinic_glucose`,
  `lactic_neutral_SSCF_cornstover`; *test_ood*) and 10 % of the samples of
  the other configurations (*test_id*). No test sample is used to fit, tune,
  select or calibrate anything. The v1 test results did, however, direct
  which failure modes round 4 investigated, so in the v2 report the test
  split is a re-used hold-out and its numbers are likely optimistic (see
  *Two runs* below).
* **One training procedure everywhere.** Every network fit (hyperparameter
  search, every CV scheme, `final_dev`, release) trains on all its samples,
  each ensemble member early-stopping on its own random 15 % sample-level
  hold-out, without a calibration hold-out. The CV numbers therefore
  estimate the training procedure that made the final models, with the
  selected hyperparameters (the hyperparameter search is not nested in the
  CV, so they are slightly optimistic for search + fit). The intervals of the
  final models are calibrated from the leave-one-configuration-out
  out-of-fold scores (section 4.5).
* **Hyperparameters** are selected by leave-configurations-out CV on the
  development set only, using the configuration-balanced MPSP error.
* **Generalization CV** (development set; grouped, so the titer-yield points
  of a sample and the samples of a held-out group never straddle the split).
  How to read the schemes:
  - *random (by sample)*: interpolation. The configuration is in training,
    with other parameter draws. This is the best case.
  - *leave-one-configuration-out*: a new configuration of known feedstocks,
    products and regimes (e.g. a known product on another feedstock or with
    another separation design). This is the typical use case.
  - *leave-one-feedstock-out*: a feedstock never seen (there are only four,
    so this is coarse).
  - *leave-one-product-out*: a new product family, never seen. This is the
    hardest case and the relevant one for a new product (section 4.3).
  - *productivity extrapolation*: trained on the central productivity range,
    tested on both tails.
* **Metrics.** The headline metric is the MAPE of the MPSP surface on the
  held-out samples' simulated titer-yield points, compared with the fit MAPE
  of the 4-coefficient equation itself (irreducible). Also reported (REPORT.md
  and `tables/cv_headline.csv`): the
  configuration-balanced MAPE (each configuration counts once), MAPE split
  into spec-determined and random samples, the noise floor, per-coefficient
  MAE/R², a conditioning-aware coefficient error (`coef_gram_rel`), interval
  coverage and constraint violations. The basis is collinear on the
  titer-yield domain, so b′ and d′ are individually ill-determined (an R² of
  d′ near or below 0 can come with an accurate surface). Judge the surface
  errors and `coef_gram_rel`, not single-coefficient R².
* **Baselines and ablations.** Nearest configuration (kNN), ridge, gradient
  boosting, plain MLP, and the PI-NN without channels, without price-explicit
  channels and without physics residuals or collocation.
* **Release model** (`models/release`). It is trained on all data and has no
  independent test. The CV and final-test estimates (`models/final_dev`)
  measure its accuracy.
* **Two runs.** The first full run (v1, commit fd3f386; report and headline
  tables in [`results/v1/`](results/v1/README.md)) was the first use of the
  test split. Its results led to the round-4 changes (training protocol,
  interval calibration, recovery inputs, productivity structure,
  extrapolation safeguards), which were chosen on the development set only
  (DESIGN.md §6, round 4). The current report (v2) uses the same frozen test
  split a second time and shows v1 and v2 side by side. Because the v1 test
  failures set the round-4 agenda, the v2 test numbers (and the v1 → v2
  test improvement) are likely optimistic and are not an independent
  confirmation; the development CV is the primary evidence, and it too is
  optimistic for the leave-one-feedstock-out and leave-one-product-out
  schemes, whose folds were used to design the changes (4 and 7 folds; at
  the fold level those gains are not statistically established, DESIGN.md
  round 4f). v2 was fixed before its test evaluation and is released
  whatever the v1 / v2 test comparison shows.

## 7. Known limitations

* **Training coverage.** 7 product families (3-HP → acrylic acid, TAL and
  TAL → potassium sorbate, succinic acid, lactic acid, ethanol, oxalic acid,
  isobutanol), 4 feedstocks, and the scale, financial and design assumptions
  of the Bioindustrial-Park modules. Yields are covered up to
  each configuration's inflection yield, and titers and productivities over
  the simulated ranges.
* **New-product extrapolation is uncertain.** Products without a close
  relative in training (in development CV: ethanol, isobutanol, oxalic,
  lactic) have errors several times those of new configurations of a known
  product. The errors come from the recovery g and from the electricity,
  natural-gas and downstream channels, which legitimately depend on the
  product, and from process structure no training configuration has
  (coproduct credits of isobutanol and corn ethanol, corn ethanol's missing
  boiler; DESIGN.md §6, rounds 1–4). The recovery of a new separation train
  is not predictable from the present descriptors (leave-one-product-out g
  errors of −23 to +44 % for every input set tried). Check
  `in_domain`/`ad_ratio`, and supply a known recovery (section 4.7) where
  possible.
* **Per-configuration bias.** The remaining error is mostly a systematic
  bias per configuration, not noise. The ensemble spread does not capture
  it, and the conformal intervals cover it only on average over
  configurations (`coverage_guaranteed` is False; section 4.5).
* **Noise floor.** Random samples differ in biorefinery-specific parameters
  (`bp_*`) that are not model inputs. Even a model that knows the
  configuration exactly cannot predict those samples better than the noise
  floor reported in REPORT.md (`evaluate.noise_floor`).
* **Primed coefficients, not process detail.** The model predicts costs at
  the biorefineries' own designs. A user recovery is a first-order
  approximation (section 4.7). Price linearity is exact only for the feedstock,
  electricity and natural-gas prices.
* **Data corrections and module quirks** (DESIGN.md §1, §3, §5).
  - Recorded data of older harness and adapter versions are corrected by
    versioned patches (`data_patches.py`: D3 boiler natural gas, D4
    sugarcane RIN credit, D6 HP lime price, D7 HP inflection yield, M1/M3/
    M12/M23 metadata, S1 downstream fields). The D2 corn-ethanol samples
    were regenerated.
  - Defects that are not corrected are tracked: HP acid regeneration (D5),
    TAL substrate basis (D8, a negative feedstock c′), lactic (D1, D9),
    isobutanol (D10), oxalic (D13, D14).
  - Known module behaviour: the succinic_neutral flowsheet acidulates but
    records an acid cost of 0, and TAL's final-product purity varies with
    titer.
  - Some cost-allocation conventions differ between modules and cannot be
    explained by the features: labor per module, the upstream capital of
    the cellulosic and sugarcane ethanol modules, and the broth-volume term
    of the HP and lactic neutralization channel.

## 8. Reproduction

From simulations to the release model (paths below are the launchers'
defaults; override them with the environment variables named in each
script):

1. **Bioindustrial-Park.** Base commit
   `b3988834106b0a31b79f9c7f01a1f4c730c58102` plus one patch series per
   biorefinery ([`data_generation/bip_patches/`](data_generation/bip_patches/README.md);
   API ports to biosteam 2.54.0 / thermosteam 0.54.2 / pandas 3 and bug
   fixes, no process assumptions changed):

   ```bash
   FI=/path/to/fermentation_insights        # this repository
   git clone https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park bip && cd bip
   for b in HP TAL succinic lactic ethanol oxalic isobutanol; do
     git worktree add -b fi-$b ../bip_$b b3988834106b0a31b79f9c7f01a1f4c730c58102
     (cd ../bip_$b && git am $FI/fermentation_insights/coeff_predictor/data_generation/bip_patches/fi-$b.patch)
   done
   ```

2. **One conda environment per biorefinery** (`fi_hp`, `fi_tal`,
   `fi_succinic`, `fi_lactic`, `fi_ethanol`, `fi_oxalic`, `fi_isobutanol`).
   Each module sets global BioSTEAM state, and the dependencies differ. Use
   Python 3.12 with biosteam 2.54.0, thermosteam 0.54.2, numpy 2.5.4,
   pandas 3.0.6, scipy, chaospy, numba, flexsolve, matplotlib and openpyxl,
   plus any biorefinery-specific requirements of Bioindustrial-Park. Isobutanol uses Python 3.14
   with nskinetics 0.6.0 and tellurium 2.2.13.1 (libroadrunner; its launcher
   sets `LD_LIBRARY_PATH`). A separate ML environment (`fi_ml`, section 3) runs
   everything after the simulations.
3. **Campaigns.** The adapters (`data_generation/adapters/`, contract in
   [`data_generation/ADAPTERS.md`](data_generation/ADAPTERS.md)) connect each
   biorefinery to the shared harness. To run or resume all campaigns:

   ```bash
   SEED=1 bash fermentation_insights/coeff_predictor/data_generation/scripts/run_all.sh
   ```

   Each launcher (`run_<name>.sh`; `ENV`/`ENV_PY`, `BIP`, `FI`,
   `OUT`/`OUT_DIR`, `LOG_DIR`, `N_RANDOM`, `N_PROD`, `N_WORKERS`) starts
   background workers. Each worker runs every configuration in its own
   Python process, writes one JSON-lines file per configuration, and skips
   completed samples when re-run. `stop_campaigns.sh` stops them. Compute
   (an estimate from the per-point times in `run_all.sh`: oxalic 0.2 s,
   ethanol 0.3 s, HP 0.5–3 s, TAL 1–3 s, isobutanol 2–4 s, lactic ~8 s,
   succinic ~14 s per titer-yield point, ~25 points per sample, with
   run_all.sh's sample counts): roughly 30 CPU-hours in total, run
   concurrently.
4. **Dataset** (optional; the experiments build it themselves). This step
   reads the raw files, applies the data patches, fits the coefficients and
   channels, and filters (minutes):

   ```bash
   python -m fermentation_insights.coeff_predictor.dataset --raw /path/to/raw --out /path/to/processed
   ```

5. **Experiments** (`fi_ml`). The pipeline runs data → frozen split →
   hyperparameter search → CV → final test + release model → checks →
   figures → REPORT.md. Every fit is cached, so an interrupted run resumes.

   ```bash
   python -m fermentation_insights.coeff_predictor.experiments --raw /path/to/raw --threads 4 --budget-hours 3 --members 5
   python -m fermentation_insights.coeff_predictor.experiments --quick     # debug run, minutes
   ```

   The search size is set so that the whole run fits `--budget-hours` on
   `--threads` single-threaded workers. The estimate is made from the
   measured cost of the first candidate; `--plan-only` writes the estimate
   to `results/plan.json` and stops. A later run reuses `plan.json` only
   while the data, split, code version, protocol, workers, members and the
   requested `--budget-hours` / `--hpo-candidates` / `--hpo-members` /
   `--ablation-loco-folds` are unchanged (otherwise it re-plans; `--replan`
   forces it). A full run can take longer than a terminal session or a
   2-hour job limit: start it detached (e.g. `setsid nohup ... &`) or simply
   relaunch the same command, which resumes from the cache. Outputs: `results/` (REPORT.md, tables,
   figures, log), `models/final_dev` and `models/release`. The CV must run
   before the final stage: the final intervals are calibrated from its
   leave-one-configuration-out out-of-fold scores.

## Package layout

| module | content |
|---|---|
| `predict.py` | `CoefficientPredictor` (user API) |
| `features.py`, `downstream.py`, `descriptors.py` | specifications → features, downstream vocabularies and table, chemical descriptors (`registry/chemicals.json`) |
| `model.py`, `train.py` | network (channels, price-explicit heads, bounds, gate, masking), loss, ensemble, conformal intervals, checkpoints |
| `physics.py`, `dataset.py`, `data_patches.py` | channel decomposition and coefficient fits, processed tables, data corrections |
| `evaluate.py`, `experiments.py` | metrics, baselines, the validation pipeline |
| `data_generation/` | harness, adapters, Bioindustrial-Park patches, campaign launchers |
| `synthetic.py`, `tests/` | synthetic data and the test suite |

## Citation

Bhagwat, S. S.; Rao, C. V.; Zhao, H.; Singh, V.; Guest, J. S. A unifying
equation for fermentation sustainability across the titer-rate-yield
landscape. *Nature Communications* 2026.
https://doi.org/10.1038/s41467-026-75285-1
