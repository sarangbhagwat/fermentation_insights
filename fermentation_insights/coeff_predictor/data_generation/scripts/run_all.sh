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
# completed <name> <target>: true if the campaign's raw output already holds
# <target> completed samples (skipping avoids reloading every configuration
# after an interruption; the launchers would only skip them one by one).
completed() {
    local n
    n=$(python3 - "$1" <<'PY'
import glob, sys
ids = set()
for f in glob.glob(f'/home/user/work/data/raw/{sys.argv[1]}/*.jsonl'):
    for l in open(f):
        if '"sample_end"' in l:
            i = l.find('"sample_id": "'); ids.add(l[i + 14:l.find('"', i + 14)])
print(len(ids))
PY
)
    if [ "$n" -ge "$2" ]; then echo "--- $1 complete ($n samples)"; return 0; fi
    return 1
}
completed oxalic 211 || run N_RANDOM=200 N_PROD=10 SEED=$SEED EXTRA_ARGS="--prod-range 0.3 8" bash "$HERE/run_oxalic.sh"
completed isobutanol 109 || run N_RANDOM=100 N_PROD=8  SEED=$SEED bash "$HERE/run_isobutanol.sh"
completed ethanol 327 || run N_RANDOM=100 N_PROD=8  SEED=$SEED bash "$HERE/run_ethanol.sh"
completed HP 464 || run N_RANDOM=24  N_PROD=4  SEED=$SEED bash "$HERE/run_HP.sh"
completed TAL 232 || run N_RANDOM=24  N_PROD=4  SEED=$SEED bash "$HERE/run_TAL.sh"
completed lactic 144 || run N_RANDOM=30  N_PROD=5  SEED=$SEED N_WORKERS=2 bash "$HERE/run_lactic.sh"
completed succinic 144 || run N_RANDOM=14  N_PROD=3  SEED=$SEED N_WORKERS=2 EXTRA_ARGS="--n-y 4 --n-t 4" bash "$HERE/run_succinic.sh"
