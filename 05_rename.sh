#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="${1:?Usage: $0 <input_dir> <output_dir>}"
OUTPUT_DIR="${2:?Usage: $0 <input_dir> <output_dir>}"

# Derive the model tag from the input directory name, e.g. "Cratons_10Ma" -> "Cratons_10Ma"
TAG="$(basename "${INPUT_DIR}")"

rename_file() {
    local input="$1"
    local output="$2"
    echo "Renaming ${input} -> ${output} ..."
    python3 rename.py "${input}" "${output}"
}

mkdir -p "${OUTPUT_DIR}"

rename_file "${INPUT_DIR}/converted.nc"                 "${OUTPUT_DIR}/${TAG}.nc"
rename_file "${INPUT_DIR}/converted_srts_filtered.nc"   "${OUTPUT_DIR}/${TAG}_S40RTS_ToFi.nc"
rename_file "${INPUT_DIR}/converted_llnl_filtered.nc"   "${OUTPUT_DIR}/${TAG}_LLNL_ToFi.nc"

echo "All done."
