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
