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

# Step 4: Interpolate the three VTU outputs onto a regular lon/lat/depth
# grid and save as NetCDF.  These .nc files are what you load for plotting
# and further analysis.
#
# Set WORK_DIR to the directory that holds your converted*.vtu files.

WORK_DIR=/Volumes/Grey/firedrake_simulations/HT/Z22/0Ma

# set -euo pipefail

# module use /g/data/fp50/modules
# module load firedrake/main-20260114
# export PYTHONPATH=/scratch/xd2/sg8812/local/lib/python3.12/site-packages:/scratch/xd2/sg8812/g-interp:${PYTHONPATH:-}

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

interp "${WORK_DIR}/converted.vtu"                  "${WORK_DIR}/converted.nc"
interp "${WORK_DIR}/converted_srts_filtered.vtu"    "${WORK_DIR}/converted_srts_filtered.nc"
interp "${WORK_DIR}/converted_tofi_filtered.vtu"    "${WORK_DIR}/converted_tofi_filtered.nc"

echo "[$(date)] All done."
