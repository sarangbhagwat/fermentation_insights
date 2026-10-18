# v1 results (first full run, 2026-10-17 06:40)

Copy of the report and headline tables of the first complete run (commit
fd3f386; full artifacts, figures and models in that commit). This run was
the first use of the frozen test split. The training protocol of the final
ensembles was revised afterwards (round 4, development set only; see
DESIGN.md), and the pipeline was re-run, so the test split was used a second
time (v2, `results/v2/`), after the data corrections and round 5 a third
time (v3, `results/v3/`) and, after the D16 / D17 simulation corrections, a
fourth time (v4). The current `results/REPORT.md` compares the four runs
(v1, v2, v3, v4).
