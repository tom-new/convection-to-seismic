#!/bin/bash
#PBS -N interpolate
#PBS -P xd2
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -l ncpus=1
#PBS -l mem=128GB
#PBS -l storage=scratch/xd2+gdata/fp50
#PBS -l wd
#PBS -j oe

# Step 4: interpolate VTU outputs onto a regular lon/lat/depth grid (NetCDF).
#
# Pass via qsub -v:
#   NAME  run identifier (required)

set -euo pipefail

: "${NAME:?must pass -v NAME=<run-id>}"

WORK=/scratch/xd2/sg8812/kat-conversion
PREFIX="${NAME}_converted"

module use /g/data/fp50/modules
module load firedrake/main-20260417
export PYTHONPATH=/scratch/xd2/sg8812/g-drift:/scratch/xd2/sg8812/local/lib/python3.11/site-packages:/scratch/xd2/sg8812/g-interp:${PYTHONPATH:-}

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
