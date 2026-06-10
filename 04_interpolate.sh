#!/usr/bin/env bash
set -euo pipefail

WORK="${1:?Usage: $0 <work_dir>}"
PREFIX="converted"

interp() {
    local vtu="$1"
    local nc="$2"
    echo "[$(date)] Interpolating ${vtu} ..."
    python3 -m ginterp.interp \
        "${vtu}" \
        --spherical \
        --radii 1.208,2.208 \
        --dims 360,181,129 \
        --output "${nc}"
}

interp "${WORK}/${PREFIX}.vtu"                  "${WORK}/${PREFIX}.nc"
interp "${WORK}/${PREFIX}_srts_filtered.vtu"    "${WORK}/${PREFIX}_srts_filtered.nc"
interp "${WORK}/${PREFIX}_llnl_filtered.vtu"    "${WORK}/${PREFIX}_llnl_filtered.nc"

echo "[$(date)] All done."
