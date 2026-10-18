#!/usr/bin/env bash
# Run a list of configurations of any biorefineries with the campaign
# settings of run_all.sh, N_WORKERS configurations at a time (each in its own
# Python process, as the per-biorefinery launchers do; the queue keeps every
# worker busy, so long configurations of one biorefinery do not serialize
# behind each other). Used for the regeneration of 2026-10-18 (harness v4,
# D16 / D17; DESIGN.md section 1).
#
# Usage:  N_WORKERS=4 bash run_queue.sh JOBS_FILE
#   JOBS_FILE: one "<adapter> <config_id>" per line, longest first (the queue
#   starts them in that order); '#' lines are skipped.
# Output, file names (tags), seeds and sample plans are those of the
# launchers called by run_all.sh, so a queue run resumes / extends the same
# raw files: $RAW_BASE/<adapter>/<file>.jsonl. Per-configuration logs go to
# $LOG_DIR/queue_<config_id>.log.
set -uo pipefail
JOBS=${1:?jobs file}
N_WORKERS=${N_WORKERS:-4}
SEED=${SEED:-1}
RAW_BASE=${RAW_BASE:-/home/user/work/data/raw}
BIP_BASE=${BIP_BASE:-/home/user/work}
ENV_BASE=${ENV_BASE:-/home/user/miniforge3/envs}
FI=${FI:-/home/user/fermentation_insights}
LOG_DIR=${LOG_DIR:-/home/user/work/logs}
mkdir -p "$LOG_DIR"
export RAW_BASE BIP_BASE ENV_BASE FI LOG_DIR SEED

job() {
    local a=$1 cid=$2 env n_random n_prod tag extra
    # settings of run_all.sh and of run_<adapter>.sh
    case "$a" in
        oxalic)     env=fi_oxalic;     n_random=200; n_prod=10; tag="_$cid"; extra="--prod-range 0.3 8" ;;
        isobutanol) env=fi_isobutanol; n_random=100; n_prod=8;  tag="_w0";  extra="" ;;
        ethanol)    env=fi_ethanol;    n_random=100; n_prod=8;  tag="_$cid"; extra="" ;;
        HP)         env=fi_hp;         n_random=24;  n_prod=4;  tag="_${cid#HP_}"; extra="" ;;
        TAL)        env=fi_tal;        n_random=24;  n_prod=4;  tag="_$cid"
                    extra="--n-y 5 --n-t 5 --point-timeout 600" ;;
        lactic)     env=fi_lactic;     n_random=30;  n_prod=5;  tag="_$cid"; extra="" ;;
        succinic)   env=fi_succinic;   n_random=14;  n_prod=3;  tag="_${cid#succinic_}"; extra="--n-y 4 --n-t 4" ;;
        *) echo "unknown adapter $a"; return 1 ;;
    esac
    local log="$LOG_DIR/queue_${cid}.log"
    echo "=== $(date -Is) start $a $cid"
    (
        cd "$FI"
        export PYTHONPATH="$BIP_BASE/bip_$a:$FI" MPLBACKEND=Agg PYTHONWARNINGS=ignore
        export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
        [ "$a" = isobutanol ] && export LD_LIBRARY_PATH="$ENV_BASE/$env/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
        nice -n 5 "$ENV_BASE/$env/bin/python" -m fermentation_insights.coeff_predictor.data_generation.run \
            --adapter "$a" --out "$RAW_BASE/$a" --configs "$cid" --n-random "$n_random" \
            --n-productivity "$n_prod" --seed "$SEED" --tag "$tag" $extra
    ) > "$log" 2>&1
    echo "=== $(date -Is) finished $a $cid (exit $?)"
}
export -f job
grep -v '^\s*#' "$JOBS" | grep -v '^\s*$' | xargs -P "$N_WORKERS" -L 1 bash -c 'job "$0" "$1"'
echo "=== $(date -Is) queue done"
