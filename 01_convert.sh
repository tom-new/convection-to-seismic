#!/usr/bin/env bash
set -euo pipefail

WORK="${1:?Usage: $0 <work_dir>}"
INPUT_PVTU="${WORK}/output/output_0.pvtu"
OUTPUT_VTU="${WORK}/converted.vtu"

echo "[$(date)] Converting ${INPUT_PVTU} -> ${OUTPUT_VTU}"
python3 "convert_to_v.py" "${INPUT_PVTU}" "${OUTPUT_VTU}"
echo "[$(date)] Done."
