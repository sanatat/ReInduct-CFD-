#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
mkdir -p v10/results
export OPENBLAS_NUM_THREADS=1
python3 -m pytest -q v4/tests v5/tests v6/tests v7/tests v8/tests v9/tests v10/tests | tee v10/results/all_tests.log
python3 -m v10.plasma_transport_1d | tee v10/results/plasma1d_summary.log
python3 -m v10.mhd_hartmann_heat | tee v10/results/hartmann_summary.log
python3 -m v10.lco_published_rates > v10/results/published_rates_summary.log
python3 -m v10.external_phases > v10/results/independent_phase_summary.log
python3 -m v10.plot_v10 | tee v10/results/plot_manifest.log
# Confirm real standalone C++ compiled and runs, but NOT OpenFOAM C++.
g++ -std=c++17 -O2 -Wall -Wextra -pedantic v8/cpp_kinetics_check.cpp -o v10/results/cpp_kinetics_check
v10/results/cpp_kinetics_check > v10/results/cpp_kinetics_log.txt
# Attempt a real OpenFOAM run; keep fail-closed status distinct from numerical tests.
set +e
bash v10/openfoam/Allrun > v10/results/openfoam_attempt.log 2>&1
foam_rc=$?
set -e
printf '%s\n' "$foam_rc" > v10/results/openfoam_exit_code.txt
python3 -m v10.evidence
if [[ "$foam_rc" -eq 0 ]]; then
  echo 'OpenFOAM execution succeeded and log checks passed. Review detailed outputs.'
else
  echo "OpenFOAM status: BLOCKED/FAILED (exit=$foam_rc); no claim of OpenFOAM verification."
fi
