#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

# Fixed, recorded defaults for the paper surface. Override only by explicitly
# setting these task-specific environment variables before invoking the script.
IV_PATHS="${IV_PATHS:-1000000}"
IV_TAU="${IV_TAU:-320}"
IV_THREADS="${IV_THREADS:-8}"
IV_SEED="${IV_SEED:-20260811}"
IV_OUTDIR="${IV_OUTDIR:-${SCRIPT_DIR}/output/production}"
IV_PYTHON="${IV_PYTHON:-python3}"
IV_COMPILER="${IV_COMPILER:-clang++}"

mkdir -p "${IV_OUTDIR}"
MPLCONFIGDIR="${MPLCONFIGDIR:-${IV_OUTDIR}/.mplconfig}"
mkdir -p "${MPLCONFIGDIR}"

"${IV_COMPILER}" -O3 -std=c++17 -pthread \
  "${SCRIPT_DIR}/generate_iv_surface.cpp" -o "${IV_OUTDIR}/generate_iv_surface"

"${IV_OUTDIR}/generate_iv_surface" \
  --paths "${IV_PATHS}" \
  --tau "${IV_TAU}" \
  --threads "${IV_THREADS}" \
  --seed "${IV_SEED}" \
  --models both \
  --output "${IV_OUTDIR}/surface_exact.csv"

MPLCONFIGDIR="${MPLCONFIGDIR}" "${IV_PYTHON}" \
  "${SCRIPT_DIR}/plot_iv_surface.py" \
  --input "${IV_OUTDIR}/surface_exact.csv" \
  --outdir "${IV_OUTDIR}"
