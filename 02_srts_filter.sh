#!/bin/bash
#PBS -N srts_filter
#PBS -P xd2
#PBS -q normal
#PBS -l walltime=02:00:00
#PBS -l ncpus=1
#PBS -l mem=64GB
#PBS -l storage=scratch/xd2+gdata/fp50
#PBS -l wd
#PBS -j oe

# Step 2: Apply S12RTS / S20RTS / S40RTS tomographic filters to Vs.
#
# INPUT_VTU must point to the output of step 1 (converted.vtu).
# The script writes <stem>_srts_filtered.vtu alongside the input automatically.

INPUT_VTU=/scratch/xd2/USERNAME/converted.vtu

set -euo pipefail

module use /g/data/fp50/modules
module load firedrake/main-20260417
# Prepend a local g-drift checkout if you need the SLB_24 pyroliteCFMASNaCr
# dataset (not yet in the installed gdrift). Edit the path or drop the entry.
export PYTHONPATH=/scratch/xd2/USERNAME/g-drift:/scratch/xd2/USERNAME/local/lib/python3.11/site-packages:${PYTHONPATH:-}

SCRIPTS_DIR="$(dirname "$(realpath "$0")")"

echo "[$(date)] Starting S-RTS tomographic filtering"
python3 "${SCRIPTS_DIR}/srts_filter.py" "${INPUT_VTU}"
echo "[$(date)] Done."
