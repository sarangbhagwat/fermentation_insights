# Bioindustrial-Park patches

Each biorefinery's data were generated on its own Bioindustrial-Park branch
(`fi-<name>`), created from Bioindustrial-Park master commit
`b3988834106b0a31b79f9c7f01a1f4c730c58102`, with minimal changes needed to run
on biosteam 2.54.0 / thermosteam 0.54.2 / pandas 3 (API ports) and a few bug
fixes that the adapters' verification uncovered. No process assumptions were
changed. To reproduce:

```bash
git clone https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park
cd Bioindustrial-Park
git checkout -b fi-HP b3988834106b0a31b79f9c7f01a1f4c730c58102
git am /path/to/fi-HP.patch
```

| branch | commits | summary |
|---|---|---|
| fi-HP | 2 | StirredTankReactor import; BatchBioreactor `V_max`; AdsorptionColumnTSA (removed from biosteam) vendored from biosteam 2.49.2; pandas 3 indexing; product-tank stream linking off (thermosteam linked-stream copy bug); shared cellulosic/corn fixes |
| fi-TAL | 4 | API ports via `TAL/_bst_compat.py`; restores pump flow linking, feed-scaled fire-water tank, autorename default, upgrading-reactor volume basis; CaSO4 locked as solid; sugarcane KSA system feeds the juice recycle to the mill once (D16, 2026-10-18) |
| fi-succinic | 1 | StirredTankReactor import; `autorename=False`; thermosteam `trim_cache` workaround; shared cellulosic/corn fixes; `load_yield` cell-mass conversion fix above inflection |
| fi-lactic | 7 | titer specifications moved onto units (nested-system specs skipped in biosteam 2.54); solver accuracy and SHF fixes; parameter fix for acetic acid yield; pandas 3 |
| fi-ethanol | 2 | cellulosic `size_batch` argument order (fermenter sizing); corn SSF `_N_ins = 2`; corn SSF mass-conserving yeast growth (`Glucose -> Yeast` by weight, as `bst.Fermentation`) |
| fi-oxalic | 1 | pandas 3; `Model(metrics=...)` → `indicators`; negative-flow guard in a polishing filter; shared HP/cellulosic fixes |
| fi-isobutanol | 2 | nskinetics pin update; beer-column Rmin guard fix |

Defects found in the final data review (DESIGN.md section 1, D16 and D17;
`../ADAPTERS.md`, "Module and BioSTEAM defects corrected by the regeneration
of 2026-10-18"); the data of the v3 models carry them, the current data do
not:

* `biorefineries/TAL/systems/system_SA_solubility_exploit_ethanol_sugarcane.py:116`
  appended the juice recycle `M201-0` to the crushing mill `U201` a second
  time (D16, TAL_SA_sugarcane). Fixed by the fourth fi-TAL commit (the line
  is commented out, as in the other TAL systems); `fi-TAL.patch` was
  re-exported and applies to the base commit above (`git am`, tree identical
  to the branch). The THF variant
  `system_SA_THF_Ethanol_solubility_exploit_ethanol_sugarcane.py`, which no
  configuration uses, keeps the upstream line.
* BioSTEAM's `InternalCirculationRx` (lumped method) counts its biogas twice
  (D17; every HP and TAL system, ethanol sugarcane, isobutanol, oxalic). It
  is in biosteam, not in Bioindustrial-Park, so no branch patches it: the
  harness applies the runtime correction `data_generation/biosteam_fixes.py`
  since harness v4.
