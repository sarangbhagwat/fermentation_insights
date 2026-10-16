#!/usr/bin/env bash
# Data-generation campaign for the TAL biorefineries (TAL and TAL -> potassium
# sorbate; 8 configurations). Configurations are distributed round-robin over
# N_WORKERS background workers; each worker runs its configurations one after
# another, each in its OWN Python process (TAL systems are module-level and
# cannot be switched within a process). Each configuration is written to its
# own JSONL file OUT_DIR/TAL_<config_id>.jsonl, so runs are resumable
# (completed samples are skipped) independently of N_WORKERS.
#
# Usage: bash run_TAL.sh            (defaults below; override via env vars,
#        N_WORKERS=2 bash run_TAL.sh  e.g., N_WORKERS=2)
set -euo pipefail

N_RANDOM=${N_RANDOM:-24}
N_PROD=${N_PROD:-4}
SEED=${SEED:-1}
N_WORKERS=${N_WORKERS:-1}
N_Y=${N_Y:-5}
N_T=${N_T:-5}
POINT_TIMEOUT=${POINT_TIMEOUT:-600}
EXTRA_ARGS=${EXTRA_ARGS:-}   # e.g. "--sample-limit 2" for smoke tests
TAG_SUFFIX=${TAG_SUFFIX:-}   # e.g. "_smoke"

ENV_PY=${ENV_PY:-/home/user/miniforge3/envs/fi_tal/bin/python}
BIP=${BIP:-/home/user/work/bip_TAL}
FI=${FI:-/home/user/fermentation_insights}
OUT_DIR=${OUT_DIR:-/home/user/work/data/raw/TAL}
LOG_DIR=${LOG_DIR:-/home/user/work/logs}

export PYTHONPATH="${BIP}:${FI}"
export PYTHONWARNINGS=ignore
mkdir -p "${OUT_DIR}" "${LOG_DIR}"

mapfile -t CONFIGS < <(cd "${FI}" && "${ENV_PY}" -m fermentation_insights.coeff_predictor.data_generation.run \
    --adapter TAL --out "${OUT_DIR}" --list-configs 2>/dev/null)
if [ "${#CONFIGS[@]}" -eq 0 ]; then echo "no TAL configurations found" >&2; exit 1; fi

for ((k=0; k<N_WORKERS; k++)); do
    MY_CONFIGS=()
    for ((i=k; i<${#CONFIGS[@]}; i+=N_WORKERS)); do MY_CONFIGS+=("${CONFIGS[$i]}"); done
    [ "${#MY_CONFIGS[@]}" -eq 0 ] && continue
    LOG="${LOG_DIR}/campaign_TAL_${k}.log"
    (
        cd "${FI}"
        for CID in "${MY_CONFIGS[@]}"; do
            echo "=== $(date -Is) worker ${k}: ${CID}"
            nice -n 5 "${ENV_PY}" -m fermentation_insights.coeff_predictor.data_generation.run \
                --adapter TAL --out "${OUT_DIR}" --configs "${CID}" \
                --n-random "${N_RANDOM}" --n-productivity "${N_PROD}" --seed "${SEED}" \
                --n-y "${N_Y}" --n-t "${N_T}" --point-timeout "${POINT_TIMEOUT}" \
                --tag "_${CID}${TAG_SUFFIX}" ${EXTRA_ARGS} \
                || echo "!!! $(date -Is) worker ${k}: ${CID} exited with status $?"
        done
        echo "=== $(date -Is) worker ${k} done"
    ) > "${LOG}" 2>&1 &
    echo "worker ${k} (pid $!): ${MY_CONFIGS[*]} -> ${LOG}"
done
