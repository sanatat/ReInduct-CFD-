#!/usr/bin/env bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$(dirname "$0")/.."
mkdir -p v7/results/cpp
OPENBLAS_NUM_THREADS=1 python3 -m pytest -q v4/tests v5/tests v6/tests v7/tests | tee v7/results/all_pytests.log
python3 -m v7.validation.nio_holdout > v7/results/kinetics_console.json
python3 -m v7.validation.interfacial_law > v7/results/kinetic_law_console.json
python3 -m v7.validation.closed_reactor > v7/results/conservation_console.json
python3 -m v7.validation.battery_phase_envelope > v7/results/battery_phase_console.json
python3 -m v7.mapping --out v7/results/em_mapping > v7/results/map_console.json
python3 -m v7.validation.mapping_convergence > v7/results/mapping_refinement_console.json
k=$(python3 - <<'PY'
from v7.validation.interfacial_law import rate_table_from_all_data
T=rate_table_from_all_data()
print(T[T[:,0]==800,1][0])
PY
)
g++ -std=c++17 -O2 -Wall -Wextra -pedantic v7/cpp_check.cpp -o v7/results/cpp/check_kinetics
v7/results/cpp/check_kinetics "$k" > v7/results/cpp/check_kinetics.json
printf 'Python tests, C++ unit smoke check, mapping and experimental holdouts: COMPLETED\n'
if ! command -v wmake >/dev/null; then
 echo 'OpenFOAM runtime: NOT VERIFIED (wmake unavailable)' | tee v7/results/openfoam_runtime_status.txt
else
 bash v7/openfoam/Allrun
fi
