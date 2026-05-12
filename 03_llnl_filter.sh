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

# Step 3: LLNL-G3D-JPS resolution matrix.
#
# Pass via qsub -v:
#   NAME  run identifier (required)

set -euo pipefail

: "${NAME:?must pass -v NAME=<run-id>}"

WORK=/scratch/xd2/sg8812/kat-conversion
INPUT_VTU="${WORK}/${NAME}_converted.vtu"
OUTPUT_VTU="${WORK}/${NAME}_converted_llnl_filtered.vtu"

module use /g/data/fp50/modules
module load firedrake/main-20260417
export PYTHONPATH=/scratch/xd2/sg8812/g-drift:/scratch/xd2/sg8812/local/lib/python3.11/site-packages:${PYTHONPATH:-}

echo "[$(date)] LLNL filtering ${INPUT_VTU} -> ${OUTPUT_VTU}"
python3 "${WORK}/llnl_filter.py" "${INPUT_VTU}" "${OUTPUT_VTU}"
echo "[$(date)] Done."
