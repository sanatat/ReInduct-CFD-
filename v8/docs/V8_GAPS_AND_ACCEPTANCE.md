# V8 verification status and research gaps (2026-11-05)

## What was executed locally

- `bash v8/run_v8.sh`: **76/76 Python tests**, independent C++17 program compiled/run, evidence plots generated.
- 2D axisymmetric electromagnetic/thermal finite-volume prototype with two-way temperature/conductivity updates, grey-wall radiation sink, gas heat-transfer closure.
- Kinetic cross-scale Ni/NiO consistency tests; ideal-gas state/molar fraction checks.
- LiCoO2 two-stage approximate benchmark on four SAME-STUDY 240-minute endpoint estimates, two training and two held-out temperatures. This is not external validation.
- Offline OpenFOAM-like temperature scalar field read and mesh-aware heating feedback mapped with power conservation.

## Gaps and exit criteria

| Gap | Status | Sufficient evidence to mark closed |
|---|---|---|
| Compiled OpenFOAM v2606 solver | **BLOCKED**: wmake/blockMesh unavailable | Real CI `wmake` compilation/link, successful `blockMesh`, `checkMesh`, complete solver logs |
| PISO and reactive CFD field coupling | **SOURCE CANDIDATE ONLY**: provisional flux projection and SI fix | Solver pressure-velocity/species/energy residuals, mesh/time-step refinement, divergence / integrated species-heat ledgers |
| NiO cross-scale rate law | **DIMENSIONAL BUG FIXED**, product-layer porosity/interface area unknown | Measured porosity, reaction area density, independent NiO packed-bed rate curves |
| H2 fraction gas state | **FIXED Tg**; compressible multi-component flow **NOT** done | Replace fixed density with coupled continuity, EOS and temperature-dependent transport/energy in tested OpenFOAM code |
| Radiation/high-temperature | **NUMERICALLY TESTED** separate EM–thermal FV and lumped OpenFOAM sink | Real wall boundary radiation/validated emissivity and temperature maps |
| EM↔CFD | **OFFLINE BIDIRECTIONAL FILE INTERFACE TESTED**, OpenFOAM coupling **NOT RUN** | Run both solvers iteratively on actual CFD mesh and verify field/energy convergence |
| Real LiCoO2/C reaction | **MECHANISM-INFORMED PROTOTYPE**, only same-study rough endpoints | Experimental, digitized time series and measured material properties; independent reactor holdouts |
| Full plasma/MHD | **NOT IMPLEMENTED** | Electron-energy/electron collision rates, charge/energy conservation, plasma benchmark |
| Moving particles / 3D CFD-DEM | **NOT IMPLEMENTED** | Particle momentum, collisions, gas/particle momentum exchange and validation as needed for target reactor |
| Independent reactor experimental validation | **NOT AVAILABLE** | Documented independent reactor measurements with uncertainty, MAE/RMSE/coverage |

## Selected quantitative evidence

See `v8/results/EVIDENCE_MANIFEST.json` for actual values (the values are refreshed on execution). Only locally executed fields and numbers are included. No OpenFOAM runtime verification is claimed.

Scientific references: Hidayat et al. 2009 DOI 10.1007/s11663-008-9212-0; 2023 LiCoO2/C kinetic paper https://pmc.ncbi.nlm.nih.gov/articles/PMC9821372/.
