#!/bin/bash
#PBS -N vs_convert
#PBS -P xd2
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -l ncpus=1
#PBS -l mem=128GB
#PBS -l storage=scratch/xd2+gdata/fp50
#PBS -l wd
#PBS -j oe

# Step 4: Convert filtered seismic velocities (Vs, Vp) to deviation from layer mean (dlnVs, dlnVp).
#
# Set WORK_DIR to the directory containing your input netCDF files.

WORK_DIR=/Volumes/Grey/firedrake_simulations/HT/Z22/0Ma

# set -euo pipefail

# module use /g/data/fp50/modules
# module load firedrake/main-20260114
# export PYTHONPATH=/scratch/xd2/sg8812/local/lib/python3.12/site-packages:${PYTHONPATH:-}

SCRIPTS_DIR="$(dirname "$(realpath "$0")")"

echo "[$(date)] Starting Vs/Vp → dlnVs/dlnVp conversion"
python3 "${SCRIPTS_DIR}/convert_to_dlnV.py" "${WORK_DIR}/converted.nc" "${WORK_DIR}/converted_dlnV.nc"
python3 "${SCRIPTS_DIR}/convert_to_dlnV.py" "${WORK_DIR}/converted_srts_filtered.nc" "${WORK_DIR}/converted_srts_filtered_dlnV.nc"
python3 "${SCRIPTS_DIR}/convert_to_dlnV.py" "${WORK_DIR}/converted_llnl_filtered.nc" "${WORK_DIR}/converted_llnl_filtered_dlnV.nc"
echo "[$(date)] Done."
