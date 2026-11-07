# ReInduct-CFD v4 — executed numerical evidence and scientific limitations

**Date:** 5 November 2026. **Scope:** research-training numerical simulation, *not industrially validated* reactor CFD. Physical constants and reference reactor geometry remain illustrative.

## Executed v4 coupled case

- 2D MAC pressure-projection incompressible Navier–Stokes, with localized T_s/X-dependent Brinkman drag.
- Reactive finite-volume H2/H2O transport and NiO solid conversion; local thermal nonequilibrium two temperatures.
- Time-dependent temperature/conversion -> drag -> momentum -> advection -> reaction feedback, a limited **two-way** coupling (not density/viscosity feedback).
- 1D magnetoquasistatic (MQS) heating re-evaluated from solid temperature every 100 time steps.
- Domain 0.12×0.02 m; 24×10 cells; dt 0.0020 s; duration 3.00 s.
- Mean reacted NiO fraction: **0.005909 (0.591%)**.
- Maximum solid temperature: **872.04 K**.
- Mean axial velocity: **0.013023 m/s**.
- Maximum discrete velocity divergence: **6.47e-16 s⁻¹**.
- H2+reacted-solid integral balance relative residual: **1.81e-16**.
- Total modeled sensible energy balance relative residual: **1.99e-14**.

> Both small residuals are tests of the implemented finite-volume bookkeeping, **not** measurements of physical accuracy. Gas chemistry, material properties, EM conductivity and kinetics are NOT calibrated for a real battery black mass.

## Grid, timestep, and ablations (all executed)

| Case | Mesh | dt [s] | Mean conversion | Max Ts [K] | Mean U [m/s] |
|---|---|---:|---:|---:|---:|
| grid_20x8 | 20×8 | 0.0020 | 0.5885% | 867.15 | 0.0129975 |
| grid_24x10 | 24×10 | 0.0020 | 0.5909% | 872.04 | 0.0130228 |
| grid_32x12 | 32×12 | 0.0020 | 0.5857% | 875.48 | 0.0127534 |
| time_dt001 | 24×10 | 0.0010 | 0.5903% | 872.05 | 0.0130223 |
| fixed_drag | 24×10 | 0.0020 | 0.5910% | 872.04 | 0.0130318 |
| no_reaction | 24×10 | 0.0020 | 0.0000% | 872.04 | 0.0130212 |
| no_induction | 24×10 | 0.0020 | 0.5786% | 850.01 | 0.0130340 |

**Interpretation:** The coarse-to-fine difference in maximum temperature is still several kelvin. This is **not a mesh-independent prediction**. Decreasing dt from 0.002 to 0.001 s causes small output changes at fixed mesh. Switching off temperature-dependent drag causes a measurable (but small) velocity/conversion difference: the coupling is present but modest under these particular conditions.

## Published kinetic benchmark — separate from reactor model

Data [µm/s] manually transcribed from *Hidayat et al. (2009),* Table III, 100% H2, DOI [10.1007/s11663-008-9212-0](https://doi.org/10.1007/s11663-008-9212-0).

An Arrhenius log-linear regression fitted **only at 500 °C and 675 °C** gives a derived apparent Eₐ = **65.85 kJ/mol**. Holdout 600 °C has absolute error **0.108 µm/s**, whereas the high-temperature (800/900/1000 °C) holdout MAE is **17.65 µm/s**. The extrapolation **fails** because the measurements are nonmonotonic; do not present this regression as validated high-T reactor kinetics. Ni-layer-growth data in µm/s cannot be inserted directly into volumetric reactor rates without morphology, specific area and a suitable rate model.

This is an independent **literature comparison**, not experimental validation of the coupled CFD field solution.

## Actual battery cathode reaction — conservative chemistry scope

The published overall carbothermal LCO pathway (RSC review, DOI [10.1039/D3CS00898C](https://doi.org/10.1039/D3CS00898C)) is:

`4 LiCoO2 + 3 C -> 2 Li2CO3 + 4 Co + CO2`

For the illustrative batch of 10 mol LiCoO2 and 7.5 mol C at assumed 60% accessible conversion:

- cobalt metal: **353.60 g**
- lithium carbonate: **221.66 g**
- residual mass balance: **2.27e-13 g**
- atomic Li/Co/C/O balances: zero by construction

**This is a stoichiometric yield calculation**, *not* a predicted or measured recovery, kinetic model, or full NMC recycling process.

## OpenFOAM execution status and audit

- An OpenCFD v2606 candidate C++ solver and 2D case are in `v4/openfoam_v2606/`.
- `bash v4/openfoam_v2606/Allrun` was attempted here; stopped with `OPENFOAM_NOT_INSTALLED: Missing blockMesh`.
- Thus **NO v4 OpenFOAM reactive run, C++ wmake compile, checkMesh result or runtime log is claimed**.
- The separate portable C++17 NiO empirical rate law **did compile and pass cross-language unit tests**.
- C++ OpenFOAM candidate is uncompiled and may need corrections in its exact API, source-term units and boundary configuration.

## What is still lacking for TU-BRIDGE

1. Install compatible OpenFOAM, compile custom solver and run checkMesh, timestep/mesh refinement, source balance and physical inlet/outlet reactor geometry.
2. Couple actual coil-geometry 2D/3D Maxwell field to Joule heating, not 1D slab heating.
3. Implement complete Eulerian gas-solid or CFD-DEM particle motion if required; current solid continuum is fixed.
4. Implement experimentally grounded black-mass chemistry (e.g., LCO/NMC/LFP) and calibrated temperature-dependent thermophysical properties, not only stoichiometry.
5. Model discharge plasma transport/ionization, electrons, and radicals only after reliable cross-validation; the separate v3 0D radical and Hartmann tests are not plasma-reactor validation.
6. Compare reaction, outlet composition, temperature, pressure drop with independent measured reactor datasets; current empirical benchmark is an independent NiO-growth dataset.

## Reproducibility

```bash
python3 -m pip install -r requirements.txt
make v4test
make v4run
make v4chem
make v4benchmark
make v4studies
```

Baseline v3 tools and tests preserved; see `README_V3_FA.md` and `docs/V3_VALIDATION.md`.
