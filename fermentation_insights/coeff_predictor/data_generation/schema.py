# -*- coding: utf-8 -*-
"""
Canonical data schema shared by all biorefinery adapters.

Every simulated TRY point is written as one row of a long-format table. Columns
fall into the groups below; adapters may add extra columns (prefixed ``x_``)
but must fill every canonical column they can (NaN otherwise).

Units
-----
yield ``y``           g fermentation product / g fermentable sugar (glucose-eq.)
titer ``t``           g fermentation product / L broth
productivity ``r``    g fermentation product / L / h
MPSP                  USD / kg pure final product
recovery ``g``        kg final product / kg fermentation product in broth
flows                 kg/h (mass), m3/h (volume), kmol/h (molar)
duties                GJ/h; electricity in MW; costs in USD/h or USD
"""

#: Canonical process areas used for the additive decomposition of MPSP.
#: Each adapter maps its own UnitGroup names onto these.
CANONICAL_AREAS = (
    'feedstock',       # feedstock acquisition, handling, storage
    'conversion',      # pretreatment, saccharification, juicing, liquefaction
    'fermentation',    # seed train, fermenters, neutralization in fermentation
    'separation',      # product recovery & purification (incl. acidulation)
    'upgrading',       # catalytic/chemical upgrading to the final product
    'wastewater',      # wastewater treatment & resource recovery
    'chp',             # boiler & turbogenerator (combined heat and power)
    'facilities',      # cooling tower, chilled water, process water, storage,
                       # heat exchanger network, other facilities
)

#: Per-area quantities recorded at every TRY point (column = f'{area}__{q}').
AREA_QUANTITIES = (
    'installed_cost',          # USD
    'material_cost',           # USD/h (purchased inlet streams)
    'heating_duty',            # GJ/h
    'cooling_duty',            # GJ/h
    'electricity_consumption', # MW
    'electricity_production',  # MW
    'utility_cost',            # USD/h (heat + power utilities)
    'installed_cost_osbl',     # USD, part of installed_cost from TEA OSBL units
)

#: Additional per-point columns written by harness v2 (generic, adapter-independent):
#: natural gas split by BioSTEAM cost pool -- 'utility' (stream utility priced
#: by ``bst.stream_prices``, part of TEA utility_cost) and 'material' (priced
#: feed stream, part of TEA material_cost):
#:   ng_flow_utility, ng_cost_utility, ng_flow_material, ng_cost_material
#:   (kg/h, USD/h);
#: labor_cost, labor_burden (TEA attributes, USD/yr and fraction, if present);
#: fs_* feedstock descriptors of the feedstock stream at that point (sample-
#: specific composition, e.g. when carbohydrate content is sampled).
#: v3: natural-gas inlet utilities include BioSTEAM's boiler 'Fuel' utility
#: (``BoilerTurbogenerator.define_utility('Fuel', ...)``, priced at
#: ``bst.stream_prices['Fuel']`` = ``BT.natural_gas_price``) and any inlet
#: utility stream >= 50 wt% CH4, priced at ``bst.stream_prices[name]`` when the
#: stream's own price is 0. v2 missed 'Fuel' and booked the boiler gas as an
#: unpriced material feed (review defect D3; v2 rows are corrected by
#: ``coeff_predictor/data_patches.py``). Adapter-level fixes are versioned by
#: the 'adapter_version' field of the configuration (config record).
HARNESS_VERSION = 3

#: Specification / design-decision columns (constant within a sample).
SPEC_COLUMNS = (
    'biorefinery',          # e.g. 'HP'
    'config_id',            # unique configuration label within the biorefinery
    'feedstock',            # e.g. 'glucose', 'corn', 'sugarcane', 'cornstover'
    'fermentation_product', # thermosteam chemical ID-like label, e.g. 'HP'
    'final_product',        # e.g. 'AcrylicAcid'
    'neutralization',       # bool: pH control by base addition in fermentation
    'separation_variant',   # label of separation design choice ('' if default)
    'aerobic',              # bool
    'theoretical_yield',    # g fp / g sugar (glucose-eq.) stoichiometric max
    'inflection_yield',     # g/g; above this the cell-mass/byproduct regime changes
    'y_min', 'y_max', 't_min', 't_max',  # TRY domain used for fitting [g/g, g/L]
    'baseline_productivity',# g/L/h
)

#: Core per-point outputs.
POINT_COLUMNS = (
    'sample_id', 'point_id', 'y', 't', 'r',
    'status',               # 'ok' | 'infeasible' | 'failed'
    'error',                # short error string
    'sim_seconds',
    'MPSP', 'recovery',
    'GWP', 'FEC',           # optional, NaN if not available
    # TEA (annual basis unless stated)
    'AOC', 'TCI', 'FCI', 'installed_cost', 'material_cost', 'utility_cost',
    'coproduct_sales', 'FOC', 'VOC', 'operating_hours', 'IRR', 'income_tax',
    # production / mass & energy balances (hourly basis)
    'F_final_product',      # kg/h pure final product
    'F_ferm_product',       # kg/h fermentation product in broth (acid-eq. mass)
    'F_sugar_fed',          # kg/h fermentable sugars (glucose-eq. mass) to fermentation
    'F_feedstock',          # kg/h wet feedstock
    'feedstock_price',      # USD/kg wet
    'V_broth',              # m3/h broth leaving fermentation
    'F_water_broth',        # kg/h water in broth
    'V_fermenters',         # m3 total fermenter volume (all vessels)
    'tau_fermentation',     # h
    'base_cost',            # USD/h base dosed in fermentation at its purchase price:
                            # neutralization or pH control (TAL NaOH, succinic low-pH
                            # lime) and HP's medium lime (bought also when not
                            # neutralizing); 0 if none. See data_generation/ADAPTERS.md
    'acid_cost',            # USD/h acid for acidulation/regeneration (0 if none)
    'heating_duty', 'cooling_duty',            # GJ/h, system totals (process units)
    'electricity_consumption', 'electricity_production',  # MW
    'natural_gas_flow',     # kg/h purchased natural gas (0 if none)
    'electricity_price', 'natural_gas_price',  # USD/kWh, USD/kg
)

ALL_AREA_COLUMNS = tuple(f'{a}__{q}' for a in CANONICAL_AREAS for q in AREA_QUANTITIES)

#: Canonical names for contextual/technological parameters shared across
#: biorefineries. Adapters map their uncertain-parameter names onto these when a
#: clear equivalent exists (columns prefixed ``p_``); all other sampled
#: parameters are stored with prefix ``bp_`` (biorefinery-specific).
CANONICAL_PARAMETERS = (
    'p_feedstock_price',          # USD/kg wet feedstock
    'p_feedstock_capacity',       # kg/h wet feedstock (or dry, see adapter doc)
    'p_electricity_price',        # USD/kWh
    'p_natural_gas_price',        # USD/kg
    'p_operating_days',           # d/y
    'p_IRR',                      # fraction
    'p_income_tax',               # fraction
    'p_boiler_efficiency',        # fraction
    'p_turbogenerator_efficiency',# fraction
    'p_cell_mass_yield',          # g cell mass / g sugar (or fraction of sugar)
    'p_base_price',               # USD/kg neutralization base
    'p_acid_price',               # USD/kg acid (e.g., sulfuric)
    'p_enzyme_loading',           # adapter units
    'p_saccharification_yield',   # fraction
)
