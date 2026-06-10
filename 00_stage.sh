#!/bin/bash
#PBS -N stage
#PBS -P xd2
#PBS -q copyq
#PBS -l walltime=04:00:00
#PBS -l ncpus=1
#PBS -l mem=8GB
#PBS -l storage=scratch/xd2+gdata/xd2
#PBS -l wd
#PBS -j oe

# Stage a Firedrake/G-ADOPT pvtu tarball into the work directory.
#
# Pass via qsub -v:
#   NAME       run identifier; used as the directory stem  (required)
#   INPUT_TAR  absolute path to a .tar.gz of the simulations output  (required)
#
# Example:
#   qsub -v NAME=C39_3e22_MuT,INPUT_TAR=/g/data/.../foo.tar.gz 00_stage.sh

set -euo pipefail

: "${NAME:?must pass -v NAME=<run-id>}"
: "${INPUT_TAR:?must pass -v INPUT_TAR=<path/to/tarball>}"

WORK=/scratch/xd2/sg8812/kat-conversion
UNTAR_DIR="${WORK}/${NAME}_output"

cd "${WORK}"

echo "[$(date)] Untarring ${INPUT_TAR} into ${UNTAR_DIR} ..."
mkdir -p "${UNTAR_DIR}"
tar xzf "${INPUT_TAR}" -C "${UNTAR_DIR}"

echo "[$(date)] Listing pvd and pvtu:"
find "${UNTAR_DIR}" -maxdepth 3 \( -name "*.pvd" -o -name "*.pvtu" \)

PVTU="${UNTAR_DIR}/output/output_0.pvtu"
if [[ -f "${PVTU}" ]]; then
    echo "[$(date)] pvtu confirmed at: ${PVTU}"
else
    echo "[$(date)] ERROR: expected pvtu not found at ${PVTU}" >&2
    exit 2
fi

echo "[$(date)] Field names in pvtu:"
grep -o "Name=\"[^\"]*\"" "${PVTU}" | sort -u

echo "[$(date)] Done."
