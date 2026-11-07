# ReInduct-CFD v7 | Evidence and acceptance report
Date: 2026-11-05

## Scope and verdict

| Requested stage | Outcome | Evidence | What remains |
|---|---|---|---|
| 1 - Compile/run OpenFOAM | Candidate solver/geometry tool/annular wedge/CI ready; **BLOCKED locally** | `openfoam_runtime_status.txt`, `openfoam/Allrun`, workflow | Run `wmake`, `blockMesh`, `checkMesh`, solver; address any real logs |
| 2 - Two-way flow/chemistry and units | Conservative molar source equations + iterative C++ port written, 0-D closed ledger verified | `closed_reactor.py`, `zeroD_conservation.json`, 58 Python tests | C++ field solver not compiled; test field-level integrated balances on true mesh |
| 3 - Induction to CFD | Conservative mapping to per-cell `qEM` + source convergence on synthetic wedge | `em_mapping/power_mapping.json`, `mapping_convergence.json` | Rerun using real OpenFOAM C,V and verify solver power and flow/temperature |
| 4 - Experimental kinetics | Table-III transcription (19 points), temperature and H2 fraction holdouts, empirical interface law and literature-based LiCoO2 product-phase scenarios | `kinetics_holdout.json`, `interfacial_law.json`, `battery_phase/phase_scenarios.json` | Independent experimental *reactor* data and calibrated surface area |

## Key local results

- Python (v4+v5+v6+v7): **58 passed**. Log: `v7/results/all_pytests.log`.
- Independent C++17 compile/run with g++: completed, hydrogen inventory residual O(1e-14), `cpp/check_kinetics.json`.
- Zero-dimensional coupled species+heat mass error: O(1e-15 mol/m3 bulk), energy error O(1e-7 J/m3) and no negative inventory.
- Source mapping on 12x48 sector: reference 30.086 W, raw interpolated 30.085 W, power-conserved 30.086 W (this reference is dependent on EM grid choice).
- Quadrature mapping convergence at 36x144: local relative volume L2 = 0.0004275; power equality enforced by normalization.
- Entire 800 deg C held out from fitting: global Arrhenius MAE 0.6344 um/s, PCHIP MAE 0.5761 um/s (4 observations). Prediction visibly overestimates at 800 C.
- H2=0.26 pressure-fraction held out: Arrhenius MAE 0.2362 um/s, PCHIP MAE 0.0784 um/s. Training and test from SAME source publication.
- LiCoO2/C high and low phase stoichiometries atom balanced; no experimental phase yield predictions.

## Important risks

1. The **OpenFOAM C++ application is not compiled or run** in this runtime. API compatibility and wedge `checkMesh` cannot be guaranteed, even with syntax-reviewed source.
2. The original v6 reaction NiO and its packed-bed parameters are an illustrative chemistry surrogate, not measured black mass. Empirical dense NiO Ni-layer rate ≠ bulk packed-bed intrinsic rate. The area density parameter is assumed.
3. The underlying EM coil distribution, conductivity and frequency are synthetic. Conserving the EM integral does not validate real induction heating.
4. Transport model assumes a stationary bed, prescribed porosity, simplified constant-density incompressibility and no moving particles, true plasma chemistry or 3D reactor.
5. No independent-laboratory reactor measurement has been used. No claim of external reactor validation or industrial-scale digital twin.

## Primary sources

- Hidayat et al. 2009 (DOI https://doi.org/10.1007/s11663-008-9212-0) Table III, reduction of *dense NiO*, original experimental rate units um/s.
- Recovery of Cobalt and Lithium by Carbothermic Reduction of LiCoO2 Cathode Material: A Kinetic Study (2023), https://pmc.ncbi.nlm.nih.gov/articles/PMC9821372/ for phase observations and two-stage mechanism.
- OpenFOAM v2606 release: https://www.openfoam.com/current-release
- OpenFOAM solverInfo: https://doc.openfoam.com/2606/tools/post-processing/function-objects/utilities/solverInfo/

## Acceptance gates for TRUE completion of stages 1-4

- A: OpenFOAM CI build has successful compiler link output and `End` plus finite residuals in solver log.
- B: `checkMesh` Mesh OK and source `qEM` count matches real mesh cells.
- C: Conservation residuals for gas+solid and energy are checked from **OpenFOAM** output, not just the 0D solver.
- D: At least one independent reactor experiment (temperature, gas composition or conversion) with uncertainty is available for true external validation.

Do not conflate test count, algebraic balance and physical predictive validity.
