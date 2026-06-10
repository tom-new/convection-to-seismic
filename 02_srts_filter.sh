#!/usr/bin/env bash
set -euo pipefail

WORK="${1:?Usage: $0 <work_dir>}"
INPUT_VTU="${WORK}/converted.vtu"

echo "[$(date)] S-RTS filtering ${INPUT_VTU}"
python3 "srts_filter.py" "${INPUT_VTU}"
echo "[$(date)] Done."
