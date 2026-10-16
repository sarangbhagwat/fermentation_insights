#!/usr/bin/env bash
# Data-generation campaign for the corn -> isobutanol (+ ethanol) biorefinery
# (adapter: adapters/isobutanol.py; TRY-imposed fermentation).
#
# Starts N_WORKERS background workers on disjoint sets of configurations
# (round-robin over `--list-configs`), each `nice -n 5`, logging to
# $LOG_DIR/campaign_isobutanol_<k>.log and writing
# $OUT/isobutanol_w<k>.jsonl. Re-running resumes (completed samples are
# skipped). The biorefinery currently has ONE configuration, so N_WORKERS > 1
# leaves the extra workers idle (they exit immediately).
#
# Usage:  bash run_isobutanol.sh            # defaults below
#         N_RANDOM=48 SEED=2 bash run_isobutanol.sh
set -euo pipefail

ENV=${ENV:-/home/user/miniforge3/envs/fi_isobutanol}
BIP=${BIP:-/home/user/work/bip_isobutanol}
FI=${FI:-/home/user/fermentation_insights}
OUT=${OUT:-/home/user/work/data/raw/isobutanol}
LOG_DIR=${LOG_DIR:-/home/user/work/logs}
N_RANDOM=${N_RANDOM:-24}
N_PROD=${N_PROD:-4}
SEED=${SEED:-1}
N_WORKERS=${N_WORKERS:-1}
EXTRA_ARGS=${EXTRA_ARGS:-}

mkdir -p "$OUT" "$LOG_DIR"
export PYTHONPATH="$BIP:$FI${PYTHONPATH:+:$PYTHONPATH}"
# libroadrunner (tellurium) pip wheel needs the env's libpython on the loader path
export LD_LIBRARY_PATH="$ENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1

MOD=fermentation_insights.coeff_predictor.data_generation.run
mapfile -t CONFIGS < <("$ENV/bin/python" -m $MOD --adapter isobutanol --out "$OUT" --list-configs)

for ((k = 0; k < N_WORKERS; k++)); do
    mine=()
    for ((i = k; i < ${#CONFIGS[@]}; i += N_WORKERS)); do mine+=("${CONFIGS[$i]}"); done
    if [ ${#mine[@]} -eq 0 ]; then
        echo "worker $k: no configurations"; continue
    fi
    nohup nice -n 5 "$ENV/bin/python" -m $MOD --adapter isobutanol --out "$OUT" \
        --configs "${mine[@]}" --n-random "$N_RANDOM" --n-productivity "$N_PROD" \
        --seed "$SEED" --tag "_w$k" $EXTRA_ARGS \
        > "$LOG_DIR/campaign_isobutanol_$k.log" 2>&1 &
    echo "worker $k (pid $!): ${mine[*]} -> $LOG_DIR/campaign_isobutanol_$k.log"
done
