#!/usr/bin/env bash
# Simulation checks (report-only; data_generation/sim_checks.py): the
# capacity response at the baseline TRY point, the InternalCirculationRx
# biogas double count (known data defect D17: BioSTEAM's original method, the
# current harness's state and the corrected method) and, with
# CHECKS=duplicate_inlets, units that receive the same stream twice (D16),
# for every configuration of every biorefinery, each configuration in its
# own Python process (the HP systems are module-level singletons), N_WORKERS
# at a time. Systems are simulated in the current harness's state (its
# BioSTEAM corrections applied, as in the campaigns). A configuration whose
# JSON already exists in RAW is skipped, so give each set of checks its own
# RAW. Then merges the per-configuration JSON files into CSV_DIR (default:
# data_generation/checks, see checks/README.md; capacity_response.csv is read
# by report_checks.capacity_response for REPORT.md section 8).
#
# Usage:  N_WORKERS=3 CHECKS="capacity ic_biogas duplicate_inlets" RAW=/home/user/work/data/sim_checks_v4 bash run_sim_checks.sh
#         N_WORKERS=3 bash run_sim_checks.sh                  # capacity + ic_biogas
#         N_WORKERS=3 CHECKS=duplicate_inlets RAW=/home/user/work/data/sim_checks_dup bash run_sim_checks.sh
#         ADAPTERS="HP TAL" CHECKS="ic_biogas" bash run_sim_checks.sh
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FI=${FI:-$(cd "$HERE/../../../.." && pwd)}
BIP_BASE=${BIP_BASE:-/home/user/work}
ENV_BASE=${ENV_BASE:-/home/user/miniforge3/envs}
RAW=${RAW:-/home/user/work/data/sim_checks}
CSV_DIR=${CSV_DIR:-$HERE/../checks}
N_WORKERS=${N_WORKERS:-2}
ADAPTERS=${ADAPTERS:-"HP TAL ethanol succinic lactic oxalic isobutanol"}
CHECKS=${CHECKS:-"capacity ic_biogas"}
ML_PY=${ML_PY:-$ENV_BASE/fi_ml/bin/python}
MOD=fermentation_insights.coeff_predictor.data_generation.sim_checks
mkdir -p "$RAW"
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1

env_of() { case "$1" in HP) echo fi_hp ;; TAL) echo fi_tal ;; *) echo "fi_$1" ;; esac; }

jobs=()
for a in $ADAPTERS; do
    py="$ENV_BASE/$(env_of "$a")/bin/python"
    bip="$BIP_BASE/bip_$a"
    while read -r cid; do
        [ -n "$cid" ] && jobs+=("$a $py $bip $cid")
    done < <(PYTHONPATH="$bip:$FI" "$py" -m fermentation_insights.coeff_predictor.data_generation.run \
                 --adapter "$a" --out /dev/null --list-configs 2>/dev/null)
done
echo "${#jobs[@]} configurations"

printf '%s\n' "${jobs[@]}" | xargs -P "$N_WORKERS" -L 1 bash -c '
    a=$0; py=$1; bip=$2; cid=$3
    out="'"$RAW"'/${cid}.json"
    [ -s "$out" ] && { echo "skip $cid"; exit 0; }
    PYTHONPATH="$bip:'"$FI"'" nice -n 5 "$py" -m '"$MOD"' --adapter "$a" --config "$cid" \
        --out "$out" --checks '"$CHECKS"' > "'"$RAW"'/${cid}.log" 2>&1
    echo "done $cid (exit $?)"
'
PYTHONPATH="$FI" "$ML_PY" -m $MOD --collect "$RAW" --csv-dir "$CSV_DIR"
