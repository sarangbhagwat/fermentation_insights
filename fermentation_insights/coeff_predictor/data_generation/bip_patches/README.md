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
| fi-TAL | 3 | API ports via `TAL/_bst_compat.py`; restores pump flow linking, feed-scaled fire-water tank, autorename default, upgrading-reactor volume basis; CaSO4 locked as solid |
| fi-succinic | 1 | StirredTankReactor import; `autorename=False`; thermosteam `trim_cache` workaround; shared cellulosic/corn fixes; `load_yield` cell-mass conversion fix above inflection |
| fi-lactic | 7 | titer specifications moved onto units (nested-system specs skipped in biosteam 2.54); solver accuracy and SHF fixes; parameter fix for acetic acid yield; pandas 3 |
| fi-ethanol | 2 | cellulosic `size_batch` argument order (fermenter sizing); corn SSF `_N_ins = 2`; corn SSF mass-conserving yeast growth (`Glucose -> Yeast` by weight, as `bst.Fermentation`) |
| fi-oxalic | 1 | pandas 3; `Model(metrics=...)` → `indicators`; negative-flow guard in a polishing filter; shared HP/cellulosic fixes |
| fi-isobutanol | 2 | nskinetics pin update; beer-column Rmin guard fix |
