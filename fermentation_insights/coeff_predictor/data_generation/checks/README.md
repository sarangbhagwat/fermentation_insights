# Simulation checks (`sim_checks.py`, `scripts/run_sim_checks.sh`)

Report-only re-simulations of every configuration, all uncertain parameters
at their baselines and the configuration's baseline productivity. Nothing
here feeds the training data. The tables of this directory are those of the
current data (harness v4, Bioindustrial-Park fi-TAL with the D16 fix,
2026-10-18); `v3/` keeps the tables of the checks that found D16 and D17
(harness v3 state, fi-TAL before the fix) and that REPORT.md section 8 of the
v3 run (`results/v3/REPORT.md`) quotes.

| file | check | content |
|---|---|---|
| `capacity_response.csv` | `capacity` | MPSP at the baseline yield and titer with `p_feedstock_capacity` x 0.5 / 0.85 / 1 / 1.15 / 2 / 3 (`harness_version` of the simulation state); read by `report_checks.capacity_response` (REPORT.md section 8) |
| `ic_biogas_points.csv`, `ic_biogas_bias.csv` | `ic_biogas` (D17) | the baseline point and a 3 x 3 TY grid of every configuration with BioSTEAM's `InternalCirculationRx` (lumped method), simulated with BioSTEAM's original method (`as_is`, as harness <= 3), in the current harness's state (`harness`) and with the corrected method (`fixed`); per configuration the MPSP change at the baseline and the lowest yield, the refitted coefficients and the biogas / IC mass gain; the 14 configurations without such a reactor are listed as not affected (`affected` = False; the criterion of `data_patches.D17_CONFIGS`) |
| `duplicate_inlets_bias.csv` | `duplicate_inlets` (D16) | units whose inlet list holds the same stream twice: none in the 41 configurations (TAL_SA_sugarcane's mill `U201` had one before the fi-TAL fix: `v3/duplicate_inlets_*.csv`) |
| `unit_bounds.csv` | `unit_bounds` (D18) | design results outside the unit class's own `_bounds` at the baseline point, with the unit's installed cost and share of the installed equipment cost (one `note` row per configuration without any) |
| `column_trays_points.csv`, `column_trays_bias.csv`, `column_trays_capacity.csv` | `unit_bounds` (D18) | configurations with a column wider than 24 ft: the baseline point and the 3 x 3 TY grid as simulated (`as_is`) and with that column's trays priced as N = ceil((D / 24 ft)^2) parallel stacks of the same area (`fixed`, a sensitivity), the MPSP change and refitted coefficients per configuration, and the capacity response (0.85 / 1 / 1.15 / 2x) in both modes |

Run of 2026-10-18 (harness v4):
`N_WORKERS=4 CHECKS="capacity ic_biogas duplicate_inlets" RAW=/home/user/work/data/sim_checks_v4 bash scripts/run_sim_checks.sh`.

Results of that run:

* `ic_biogas`: 27 configurations affected (all HP and TAL, ethanol_sugarcane,
  IBO_EtOH_corn, oxalic_neutral_sugarcane), 14 not. In the `harness` state
  the biogas is half the `as_is` biogas (0.50-0.52) and equals `fixed` (MPSP within
  0.01 % at every point: `dMPSP_harness_vs_fixed_max_abs_pct`; the IC mass
  gain at the baseline is −47 to +800 kg/h instead of +4.6 kg/h to
  +18.8 t/h). The `as_is` and `fixed` columns reproduce `v3/` (baseline MPSP changes
  within 0.001 percentage points) except for TAL_SA_sugarcane, whose system no longer has
  the D16 duplicate inlet (baseline +1.42 % instead of +1.48 %; its IC
  biogas, 18 960 kg/h as is, now equals TAL_sugarcane's). One grid point
  each of TAL_glucose (`harness`) and TAL_SA_glucose (`fixed`) did not
  converge (the wastewater recycle at the maximum-titer edge, DESIGN.md H1).
* `duplicate_inlets`: no unit with a repeated inlet in the 41
  configurations.
* `capacity`: 216 points of 36 configurations, all converged; the MPSPs of
  the configurations without an IC reactor and of oxalic equal `v3/`,
  those of HP, TAL and IBO_EtOH_corn are higher by the D17 (and, for
  TAL_SA_sugarcane, D16) correction (+0.5 to +16 %).

Run of 2026-10-18 (harness v4), check `unit_bounds` (D18):
`N_WORKERS=3 CHECKS=unit_bounds RAW=/home/user/work/data/sim_checks_bounds bash scripts/run_sim_checks.sh`.

* `unit_bounds`: every HP configuration has a purification column wider
  than BioSTEAM's 24 ft bound at the baseline (D408 27.5-40.2 ft; D401
  29.3-40.3 ft in the 8 hexanol designs; 4-16 % of the installed equipment
  cost), and D409 is 350 ft tall (bound 170 ft). No column of the other 25
  configurations exceeds an upper bound; the lactic `D401` and some ethanol /
  isobutanol columns are below the lower height / weight bounds (small
  columns, costed by the same correlations at the low end).
* `column_trays_*` (the 16 HP configurations): with the trays of a column
  wider than 24 ft priced as parallel stacks of the same area, the baseline
  MPSP falls by 3.9-4.1 % (HP_hexanol_sugarcane, HP_neutral_hexanol_sugarcane),
  0.9-1.0 % (HP_sugarcane, HP_neutral_sugarcane), 0.5-0.8 % (the other
  hexanol designs) and at most 0.13 % (the rest); 0 at the lowest grid
  yield. The as-is FCI grows x1.12-1.17 for x1.15 capacity, and the MPSP of
  the two hexanol-sugarcane configurations rises with the capacity (+0.7 /
  +0.8 % at 1.15x, +20 / +21 % at 2x); with parallel stacks every HP
  configuration's MPSP falls with the capacity (−0.6 to −1.5 % at 1.15x,
  −2 to −6 % at 2x).
