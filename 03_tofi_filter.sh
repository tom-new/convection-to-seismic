#!/bin/bash
#PBS -N tofi_filter
#PBS -P xd2
#PBS -q normal
#PBS -l walltime=02:00:00
#PBS -l ncpus=1
#PBS -l mem=64GB
#PBS -l storage=scratch/xd2+gdata/fp50
#PBS -l wd
#PBS -j oe

# Step 3: Apply the LLNL-G3D-JPS resolution matrix to Vs and Vp.
#
# INPUT_VTU  - output of step 1 (converted.vtu)
# OUTPUT_VTU - where to write the LLNL-filtered result

INPUT_VTU=/Volumes/Grey/firedrake_simulations/HT/Z22/0Ma/converted.vtu
OUTPUT_VTU=/Volumes/Grey/firedrake_simulations/HT/Z22/0Ma/converted_tofi_filtered.vtu

# set -euo pipefail

# module use /g/data/fp50/modules
# module load firedrake/main-20260114
# export PYTHONPATH=/scratch/xd2/sg8812/local/lib/python3.12/site-packages:${PYTHONPATH:-}

SCRIPTS_DIR="$(dirname "$(realpath "$0")")"

echo "[$(date)] Starting LLNL tomographic filtering"
python3 "${SCRIPTS_DIR}/tofi_filter.py" "${INPUT_VTU}" "${OUTPUT_VTU}"
echo "[$(date)] Done."
