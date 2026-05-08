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

# Step 1: Convert mantle convection output (PVTU) to seismic velocities (Vs, Vp).
#
# Edit the two paths below before submitting:
#   INPUT_PVTU  - path to the output_*.pvtu file from your simulation
#   OUTPUT_VTU  - where to write the converted result (a single .vtu file)

INPUT_PVTU=/scratch/xd2/USERNAME/output/output_0.pvtu
OUTPUT_VTU=/scratch/xd2/USERNAME/converted.vtu

set -euo pipefail

module use /g/data/fp50/modules
module load firedrake/main-20260417
# Prepend a local g-drift checkout if you need the SLB_24 pyroliteCFMASNaCr
# dataset (not yet in the installed gdrift). Edit the path or drop the entry.
export PYTHONPATH=/scratch/xd2/USERNAME/g-drift:/scratch/xd2/USERNAME/local/lib/python3.11/site-packages:${PYTHONPATH:-}

SCRIPTS_DIR="$(dirname "$(realpath "$0")")"

echo "[$(date)] Starting temperature → Vs/Vp conversion"
python3 "${SCRIPTS_DIR}/convert_to_vs.py" "${INPUT_PVTU}" "${OUTPUT_VTU}"
echo "[$(date)] Done."
