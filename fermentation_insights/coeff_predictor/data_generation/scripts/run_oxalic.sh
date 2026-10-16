#!/usr/bin/env bash
# Data-generation campaign for the oxalic acid biorefinery
# (adapter: adapters/oxalic.py; sugarcane, lime-neutralized fermentation,
# clarified-broth product).
#
# Starts N_WORKERS background workers on disjoint sets of configurations
# (round-robin over `--list-configs`), each `nice -n 5`, logging to
# $LOG_DIR/campaign_oxalic_<k>.log. Each worker runs its configurations
# sequentially, each in its own python process, writing
# $OUT/oxalic_<config_id>.jsonl. Re-running resumes (completed samples are
# skipped).
#
# Usage:  bash run_oxalic.sh              # defaults below
#         N_WORKERS=2 N_RANDOM=48 SEED=2 bash run_oxalic.sh
set -uo pipefail

ENV=${ENV:-/home/user/miniforge3/envs/fi_oxalic}
BIP=${BIP:-/home/user/work/bip_oxalic}
FI=${FI:-/home/user/fermentation_insights}
OUT=${OUT:-/home/user/work/data/raw/oxalic}
LOG_DIR=${LOG_DIR:-/home/user/work/logs}
N_RANDOM=${N_RANDOM:-24}
N_PROD=${N_PROD:-4}
SEED=${SEED:-1}
N_WORKERS=${N_WORKERS:-1}
EXTRA_ARGS=${EXTRA_ARGS:-}

mkdir -p "$OUT" "$LOG_DIR"
export PYTHONPATH="$BIP:$FI${PYTHONPATH:+:$PYTHONPATH}"
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1

MOD=fermentation_insights.coeff_predictor.data_generation.run
mapfile -t CONFIGS < <("$ENV/bin/python" -m $MOD --adapter oxalic --out "$OUT" --list-configs 2>/dev/null | grep -v '^TEA breakdown')

for ((k = 0; k < N_WORKERS; k++)); do
    mine=()
    for ((i = k; i < ${#CONFIGS[@]}; i += N_WORKERS)); do mine+=("${CONFIGS[$i]}"); done
    if [ ${#mine[@]} -eq 0 ]; then
        echo "worker $k: no configurations"; continue
    fi
    LOG="$LOG_DIR/campaign_oxalic_$k.log"
    nohup nice -n 5 bash -c '
        MOD=$1; ENV=$2; OUT=$3; N_RANDOM=$4; N_PROD=$5; SEED=$6; EXTRA_ARGS=$7; shift 7
        for cid in "$@"; do
            echo "=== $(date -Is) worker config $cid"
            "$ENV/bin/python" -m $MOD --adapter oxalic --out "$OUT" \
                --configs "$cid" --n-random "$N_RANDOM" --n-productivity "$N_PROD" \
                --seed "$SEED" --tag "_$cid" $EXTRA_ARGS
            echo "=== $(date -Is) finished $cid (exit $?)"
        done
    ' _ "$MOD" "$ENV" "$OUT" "$N_RANDOM" "$N_PROD" "$SEED" "$EXTRA_ARGS" "${mine[@]}" \
        > "$LOG" 2>&1 &
    echo "worker $k (pid $!): ${mine[*]} -> $LOG"
done
