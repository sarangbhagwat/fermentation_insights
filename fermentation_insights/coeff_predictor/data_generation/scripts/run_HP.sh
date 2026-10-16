#!/usr/bin/env bash
# Data-generation campaign for the 3-HP -> acrylic acid biorefineries
# (adapter: adapters/HP.py; 16 configurations = {HP (low pH), HP_neutral,
# HP_hexanol, HP_neutral_hexanol} x {glucose, corn, sugarcane, cornstover}).
#
# Starts N_WORKERS background workers on disjoint sets of configurations
# (round-robin over `--list-configs`), each `nice -n 5`, logging to
# $LOG_DIR/campaign_HP_<k>.log. Each worker runs its configurations
# sequentially, EACH IN ITS OWN PYTHON PROCESS (the HP systems are
# module-level singletons; only one system module can be loaded per process),
# writing $OUT/<config_id>.jsonl (e.g. HP_neutral_cornstover.jsonl). Re-running
# resumes (completed samples are skipped).
#
# Usage:  bash run_HP.sh              # defaults below
#         N_WORKERS=2 N_RANDOM=48 SEED=2 bash run_HP.sh
set -uo pipefail

N_RANDOM=${N_RANDOM:-24}
N_PROD=${N_PROD:-4}
SEED=${SEED:-1}
N_WORKERS=${N_WORKERS:-1}
ENV=${ENV:-/home/user/miniforge3/envs/fi_hp}
BIP=${BIP:-/home/user/work/bip_HP}
FI=${FI:-/home/user/fermentation_insights}
OUT=${OUT:-/home/user/work/data/raw/HP}
LOG_DIR=${LOG_DIR:-/home/user/work/logs}
EXTRA_ARGS=${EXTRA_ARGS:-}

mkdir -p "$OUT" "$LOG_DIR"
export PYTHONPATH="$BIP:$FI${PYTHONPATH:+:$PYTHONPATH}"
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1

MOD=fermentation_insights.coeff_predictor.data_generation.run
# CONFIGS_OVERRIDE="cfg1 cfg2 ..." restricts the run to those configurations
# (e.g. to add workers on configurations no running worker will touch);
# WORKER_OFFSET shifts the worker/log index so existing logs are not reused.
if [ -n "${CONFIGS_OVERRIDE:-}" ]; then
    read -ra CONFIGS <<< "$CONFIGS_OVERRIDE"
else
    mapfile -t CONFIGS < <("$ENV/bin/python" -m $MOD --adapter HP --out "$OUT" --list-configs)
fi

for ((k = 0; k < N_WORKERS; k++)); do
    mine=()
    for ((i = k; i < ${#CONFIGS[@]}; i += N_WORKERS)); do mine+=("${CONFIGS[$i]}"); done
    if [ ${#mine[@]} -eq 0 ]; then
        echo "worker $k: no configurations"; continue
    fi
    LOG="$LOG_DIR/campaign_HP_$((k + ${WORKER_OFFSET:-0})).log"
    nohup nice -n 5 bash -c '
        MOD=$1; ENV=$2; OUT=$3; N_RANDOM=$4; N_PROD=$5; SEED=$6; EXTRA_ARGS=$7; shift 7
        for cid in "$@"; do
            echo "=== $(date -Is) worker config $cid"
            "$ENV/bin/python" -m $MOD --adapter HP --out "$OUT" \
                --configs "$cid" --n-random "$N_RANDOM" --n-productivity "$N_PROD" \
                --seed "$SEED" --tag "_${cid#HP_}" $EXTRA_ARGS
            rc=$?
            echo "=== $(date -Is) finished $cid (exit $rc)"
        done
    ' _ "$MOD" "$ENV" "$OUT" "$N_RANDOM" "$N_PROD" "$SEED" "$EXTRA_ARGS" "${mine[@]}" \
        > "$LOG" 2>&1 &
    echo "worker $k (pid $!): ${mine[*]} -> $LOG"
done
