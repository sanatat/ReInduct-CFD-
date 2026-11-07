#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
mkdir -p v9/results
python3 -m pytest -q v4/tests v5/tests v6/tests v7/tests v8/tests v9/tests | tee v9/results/full_pytest.log
python3 -m v9.plasma_energy > v9/results/plasma_console.json
python3 -m v9.lco_process > v9/results/lco_console.json
python3 -m v9.validation_protocol > v9/results/kinetics_console.json
python3 -m v9.openfoam.preflight | tee v9/results/openfoam_preflight.log
# C++ scientific programming independently compiled and run (NOT OpenFOAM C++).
g++ -std=c++17 -O2 -Wall -Wextra v8/cpp_kinetics_check.cpp -o v9/results/cpp_kinetics_check
v9/results/cpp_kinetics_check | tee v9/results/cpp.log
python3 -m v9.evidence
