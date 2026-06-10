#!/usr/bin/env bash
# pipeline.sh — extract and post-process all Cratons_XMa.tar.gz tarballs
#
# Usage:
#   ./pipeline.sh <base_dir> <output_dir> [scripts_dir]
#                 [--min-age <Ma>] [--max-age <Ma>]
#
#   base_dir    — directory containing Cratons_XMa.tar.gz files
#                 and where XMa/ subdirectories will be created
#   output_dir  — destination for the final renamed .nc files
#                 (passed to 05_rename.sh as its second argument)
#   scripts_dir — directory containing the 01–05 scripts and Python helpers
#                 (defaults to the directory of this script)
#   --min-age   — only process tarballs with age >= this value (Ma)
#   --max-age   — only process tarballs with age <= this value (Ma)
#
# Examples:
#   ./pipeline.sh /Volumes/Grey/firedrake_simulations/Cratons \
#                 ~/OneDrive/phd/firedrake-models/Cratons
#
#   ./pipeline.sh /Volumes/Grey/firedrake_simulations/Cratons \
#                 ~/OneDrive/phd/firedrake-models/Cratons \
#                 --min-age 10 --max-age 50

set -euo pipefail

# ── arguments ────────────────────────────────────────────────────────────────
BASE_DIR="${1:?Usage: $0 <base_dir> <output_dir> [scripts_dir] [--min-age N] [--max-age N]}"
OUTPUT_DIR="${2:?Usage: $0 <base_dir> <output_dir> [scripts_dir] [--min-age N] [--max-age N]}"
shift 2

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MIN_AGE=""
MAX_AGE=""

# ── parse remaining args ──────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --min-age)
            MIN_AGE="${2:?--min-age requires a value}"
            if ! [[ "${MIN_AGE}" =~ ^[0-9]+$ ]]; then
                echo "ERROR: --min-age must be a non-negative integer." >&2; exit 1
            fi
            shift 2 ;;
        --max-age)
            MAX_AGE="${2:?--max-age requires a value}"
            if ! [[ "${MAX_AGE}" =~ ^[0-9]+$ ]]; then
                echo "ERROR: --max-age must be a non-negative integer." >&2; exit 1
            fi
            shift 2 ;;
        *)
            # Treat first unrecognised positional arg as scripts_dir
            if [[ -z "${SCRIPTS_DIR_OVERRIDE:-}" ]]; then
                SCRIPTS_DIR="$1"
                SCRIPTS_DIR_OVERRIDE=1
            else
                echo "ERROR: Unexpected argument '$1'." >&2; exit 1
            fi
            shift ;;
    esac
done

if [[ -n "${MIN_AGE}" && -n "${MAX_AGE}" && "${MIN_AGE}" -gt "${MAX_AGE}" ]]; then
    echo "ERROR: --min-age (${MIN_AGE}) must be <= --max-age (${MAX_AGE})." >&2
    exit 1
fi

# ── sanity checks ─────────────────────────────────────────────────────────────
for script in 01_convert.sh 02_srts_filter.sh 03_llnl_filter.sh \
              04_interpolate.sh 05_rename.sh; do
    if [[ ! -f "${SCRIPTS_DIR}/${script}" ]]; then
        echo "ERROR: ${SCRIPTS_DIR}/${script} not found." >&2
        exit 1
    fi
done

shopt -s nullglob
tarballs=("${BASE_DIR}"/Cratons_*Ma.tar.gz)
shopt -u nullglob

if [[ ${#tarballs[@]} -eq 0 ]]; then
    echo "No tarballs matching ${BASE_DIR}/Cratons_*Ma.tar.gz found." >&2
    exit 1
fi

# ── per-tarball pipeline ──────────────────────────────────────────────────────
skipped=0
processed=0

for tarball in "${tarballs[@]}"; do
    filename="$(basename "${tarball}")"
    age="${filename#Cratons_}"
    age="${age%.tar.gz}"
    age_int="${age%Ma}"

    if ! [[ "${age_int}" =~ ^[0-9]+$ ]]; then
        echo "WARNING: Could not parse age from '${filename}', skipping." >&2
        (( skipped++ )) || true
        continue
    fi

    # ── age range filter ──────────────────────────────────────────────────────
    if [[ -n "${MIN_AGE}" && "${age_int}" -lt "${MIN_AGE}" ]]; then
        echo "Skipping ${filename} (age ${age_int} Ma < min ${MIN_AGE} Ma)"
        (( skipped++ )) || true
        continue
    fi
    if [[ -n "${MAX_AGE}" && "${age_int}" -gt "${MAX_AGE}" ]]; then
        echo "Skipping ${filename} (age ${age_int} Ma > max ${MAX_AGE} Ma)"
        (( skipped++ )) || true
        continue
    fi

    WORK="${BASE_DIR}/${age}"

    echo ""
    echo "════════════════════════════════════════════════════════"
    echo "  Processing ${filename}  ->  ${WORK}"
    echo "════════════════════════════════════════════════════════"

    # ── step 0: extract ───────────────────────────────────────────────────────
    echo "[$(date)] Extracting ${tarball} -> ${WORK}"
    mkdir -p "${WORK}"
    gtar -xzf "${tarball}" -C "${WORK}"

    # ── steps 1–5 ─────────────────────────────────────────────────────────────
    bash "${SCRIPTS_DIR}/01_convert.sh"     "${WORK}"
    bash "${SCRIPTS_DIR}/02_srts_filter.sh" "${WORK}"
    bash "${SCRIPTS_DIR}/03_llnl_filter.sh" "${WORK}"
    bash "${SCRIPTS_DIR}/04_interpolate.sh" "${WORK}"
    bash "${SCRIPTS_DIR}/05_rename.sh"      "${WORK}" "${OUTPUT_DIR}"

    echo "[$(date)] Finished ${filename}"
    (( processed++ )) || true
done

echo ""
echo "[$(date)] Done — ${processed} processed, ${skipped} skipped."