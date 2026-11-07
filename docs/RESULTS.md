# Execution record — 3 November 2026

The C++ code in this repository was compiled locally with `g++ -std=c++17 -O2 -Wall -Wextra -pedantic`. Runs complete successfully. This page reports actual **synthetic model outputs**, not experimental facts or validated industrial predictions. `OpenFOAM` was unavailable in the execution environment and no OpenFOAM CFD case has been executed.

## Baseline

Command: `./build/reactor1d --out results/base`

| Quantity | Computed value |
|---|---:|
| 1D FV cells | 60 |
| Time interval | 120 seconds |
| Time step | 0.01 seconds |
| Mean solid conversion | 0.6537382907 |
| Mean solid temperature | 884.5730256 K |
| Maximum solid temperature | 907.1472802 K |
| Absolute integrated H2 balance residual | 6.87095 × 10^-12 mol/m² |
| Analytical skin depth for assumptions in source | 0.00159155 m |

See `results/base/temperature.png`, `conversion.png`, `conversion_history.png`.

## Grid-refinement check

Fixed dt = 0.005 seconds, T = 120 seconds.

| Cells | Mean X | Maximum solid T (K) |
|---:|---:|---:|
| 30 | 0.6520451274 | 907.1282173 |
| 60 | 0.6534747968 | 907.1452388 |
| 120 | 0.6541942299 | 907.1489913 |

Difference in mean conversion between 60 and 120 cells: 0.0007194 absolute (~0.11% of the finest-grid value). This comparison checks grid sensitivity, but is **not** a mathematically certified error bound.

## Time-step sensitivity

120 cells fixed, T = 120 seconds.

| dt (s) | Mean X | Maximum solid T (K) |
|---:|---:|---:|
| 0.01 | 0.6544578 | 907.1512 |
| 0.005 | 0.6541942 | 907.1490 |
| 0.0025 | 0.6540624 | 907.1479 |

The trends are consistent with first-order explicit time stepping. Difference in mean conversion between dt=0.005 and dt=0.0025 is ~0.000132 absolute.

## Illustrative sensitivity to flow and heat source

Simulation results from `python3 scripts/parameter_study.py`:

| Gas superficial velocity (m/s) | 300 kW/m³ | 600 kW/m³ | 900 kW/m³ |
|---:|---:|---:|---:|
| 0.008 | 51.50% | 54.48% | 57.32% |
| 0.016 | 61.36% | 65.37% | 69.28% |
| 0.032 | 66.09% | 70.36% | 74.44% |

Numbers are **outputs of a deliberately illustrative simplified physical model**, with unfitted kinetic and heating coefficients; do not generalize these trends to industrial pyrolysis or black-mass recovery. Visualization in `results/parameter_study/conversion_sensitivity.png`.

## Passed programmatic tests

- Species concentrations nonnegative.
- Solid conversion stays in [0,1].
- Inlet H2 accounting closes against outlet, void inventory and reduced oxide within 1e-9 relative tolerance.
- `A=0` implies no reaction and no H2O produced.
- `A=0, heat=0` and uniform temperature maintain initial thermal equilibrium.
- Increasing prescribed heat input increases peak temperature and conversion under the fixed illustrative model.

## Remaining blockers before claiming expertise in TU-BRIDGE requirements

1. No genuine CFD velocity-pressure calculation yet.
2. No genuine multi-region or multiphase OpenFOAM case yet.
3. No electromagnetic field solve, Lorentz coupling or MHD model; skin depth is analytic reference only.
4. No induction/plasma heating calculated from fields or radical plasma chemistry.
5. No independent experimental validation or actual black-mass composition.

Do not describe these as completed in a CV.
