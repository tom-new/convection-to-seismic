#!/usr/bin/env bash
set -euo pipefail

WORK="${1:?Usage: $0 <work_dir>}"
INPUT_VTU="${WORK}/converted.vtu"
OUTPUT_VTU="${WORK}/converted_llnl_filtered.vtu"

echo "[$(date)] LLNL filtering ${INPUT_VTU} -> ${OUTPUT_VTU}"
python3 "llnl_filter.py" "${INPUT_VTU}" "${OUTPUT_VTU}"
echo "[$(date)] Done."
