INPUT_DIR=/Volumes/Grey/firedrake_simulations/HT/Z22/0Ma
OUTPUT_DIR=~/OneDrive/phd/firedrake-models/HT_Z22_sia

rename () {
    local input="$1"
    local output="$2"
    echo "Renaming ${input} → ${output} ..."
    python3 rename.py "${input}" "${output}"
}

rename "${INPUT_DIR}/converted.nc"                  "${OUTPUT_DIR}/HT_Z22.nc"
rename "${INPUT_DIR}/converted_srts_filtered.nc"    "${OUTPUT_DIR}/HT_Z22_S40RTS_ToFi.nc"
rename "${INPUT_DIR}/converted_tofi_filtered.nc"    "${OUTPUT_DIR}/HT_Z22_LLNL_ToFi.nc"

echo "All done."