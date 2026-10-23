#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT_DIR}"

KINE_FILE="input/o2sim_Kine.root"
HITS_FILE="input/o2sim_HitsFOC.root"
OUTPUT_FILE="output/Results.root"
GEO_FILE="../FOCAL/simulation/geometryFiles/geometry_Sheets.txt"

make -C "${ROOT_DIR}"

mkdir -p output
rm -rf output/*

if [[ ! -f O2FoCalGenQC.so ]]; then
  echo "ERROR: O2FoCalGenQC.so was not built in ${ROOT_DIR}" >&2
  exit 1
fi
if [[ ! -f "${KINE_FILE}" ]]; then
  echo "ERROR: kine file not found: ${KINE_FILE}" >&2
  exit 1
fi
if [[ ! -f "${HITS_FILE}" ]]; then
  echo "ERROR: hits file not found: ${HITS_FILE}" >&2
  exit 1
fi
if [[ ! -f "${GEO_FILE}" ]]; then
  echo "ERROR: geometry file not found: ${GEO_FILE}" >&2
  exit 1
fi

# runGenQC(kine, hits, output, geometry)
root -l -b -q "macro/runGenQC.C(\"${KINE_FILE}\",\"${HITS_FILE}\",\"${OUTPUT_FILE}\",\"${GEO_FILE}\")"
