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

INPUT_VTU=/scratch/xd2/USERNAME/converted.vtu
OUTPUT_VTU=/scratch/xd2/USERNAME/converted_tofi_filtered.vtu

set -euo pipefail

module use /g/data/fp50/modules
module load firedrake/main-20260417
# Prepend a local g-drift checkout if you need the SLB_24 pyroliteCFMASNaCr
# dataset (not yet in the installed gdrift). Edit the path or drop the entry.
export PYTHONPATH=/scratch/xd2/USERNAME/g-drift:/scratch/xd2/USERNAME/local/lib/python3.11/site-packages:${PYTHONPATH:-}

SCRIPTS_DIR="$(dirname "$(realpath "$0")")"

echo "[$(date)] Starting LLNL tomographic filtering"
python3 "${SCRIPTS_DIR}/tofi_filter.py" "${INPUT_VTU}" "${OUTPUT_VTU}"
echo "[$(date)] Done."
