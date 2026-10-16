#!/usr/bin/env bash
# Launch (or resume, after an interruption) every data-generation campaign
# used to train coeff_predictor. Each launcher skips samples already
# completed in its JSON-lines output, so re-running this script is safe.
#
# Settings balance the number of samples per biorefinery against the cost of
# one titer-yield (TY) point (median seconds/point on a 4-core machine):
#   oxalic 0.2 s, ethanol 0.3 s, TAL ~1-3 s, HP 0.5-3 s, isobutanol 2-4 s,
#   lactic ~8 s, succinic ~14 s.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SEED=${SEED:-1}
run() { echo "--- $*"; env "$@"; }
run N_RANDOM=200 N_PROD=10 SEED=$SEED EXTRA_ARGS="--prod-range 0.3 8" bash "$HERE/run_oxalic.sh"
run N_RANDOM=100 N_PROD=8  SEED=$SEED bash "$HERE/run_isobutanol.sh"
run N_RANDOM=100 N_PROD=8  SEED=$SEED bash "$HERE/run_ethanol.sh"
run N_RANDOM=24  N_PROD=4  SEED=$SEED bash "$HERE/run_HP.sh"
run N_RANDOM=24  N_PROD=4  SEED=$SEED bash "$HERE/run_TAL.sh"
run N_RANDOM=30  N_PROD=5  SEED=$SEED N_WORKERS=2 bash "$HERE/run_lactic.sh"
run N_RANDOM=14  N_PROD=3  SEED=$SEED N_WORKERS=2 EXTRA_ARGS="--n-y 4 --n-t 4" bash "$HERE/run_succinic.sh"
