#!/usr/bin/env bash
# Stop data-generation campaigns. Kills the launcher loops first (otherwise a
# loop would start its next configuration as soon as its worker dies), then
# the python workers. Usage: bash stop_campaigns.sh [adapter-name ...]
# (no argument = all adapters). Completed samples are kept; re-running a
# launcher resumes.
set -u
names=("$@"); [ ${#names[@]} -eq 0 ] && names=("")
for n in "${names[@]}"; do
    pat="scripts/run_${n}"
    pgrep -f "$pat" | xargs -r kill 2>/dev/null
done
sleep 1
for n in "${names[@]}"; do
    pat="data_generation.run --adapter ${n}"
    pgrep -f "$pat" | xargs -r kill 2>/dev/null
done
sleep 3
echo "remaining workers: $(pgrep -fc 'data_generation.run')"
