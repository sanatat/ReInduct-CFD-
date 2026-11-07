#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
mkdir -p v8/results
python3 -m pytest -q v4/tests v5/tests v6/tests v7/tests v8/tests | tee v8/results/all_pytests.log
python3 -m v8.evidence | tee v8/results/evidence_console.json
python3 -m v7.validation.interfacial_law > v8/results/nickel_interface_table_console.json
# Compile a real, independent C++ kinetics test, but DO NOT confuse it with the OpenFOAM app.
g++ -std=c++17 -O2 -Wall -Wextra v8/cpp_kinetics_check.cpp -o v8/results/independent_cpp_check
v8/results/independent_cpp_check | tee v8/results/independent_cpp_check.log
if command -v wmake >/dev/null && command -v blockMesh >/dev/null; then
  printf 'OpenFOAM commands present: runtime must be separately validated with bash v8/openfoam/Allrun\n' > v8/results/openfoam_status.txt
else
  printf 'OpenFOAM NOT VERIFIED: wmake/blockMesh unavailable in this environment.\n' > v8/results/openfoam_status.txt
fi
printf 'v8 numerical tests and C++ stand-alone check completed; OpenFOAM status is in v8/results/openfoam_status.txt\n'
