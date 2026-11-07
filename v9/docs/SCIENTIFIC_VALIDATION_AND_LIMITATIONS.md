# ReInduct-CFD v9: scientific validation, scope and remaining blockers

**Date:** 2026-11-06. **Target:** TU-BRIDGE DP.MUL.3. This file is a record of executed evidence and excludes unverified claims.

## 1. What is numerically verified

- Re-executed all Python unit/regression tests (v4–v9), with separately compiled and executed standalone C++17 kinetics check.
- Introduced a 0D energy-accounted effective H₂ ↔ 2H chemistry model. Atom inventory `2[H2]+[H]` and energy `Ee+Eg+(D/2)[H]` balance the imposed electron power. This model is deliberately **not** a plasma discharge model; it omits electron-energy distribution, ion chemistry, electric potential, plasma transport and gas-phase velocity.
- Introduced kinetic LiCoO₂/C batch **scenarios** (liquid-free chemical recycling) with constrained carbon availability and temperature-dependent observed product phase classes: Li₂CO₃+Co at 700–900°C; Li₂O+Co at 1000–1100°C. The high-temperature CO gas speciation is assumed solely to close a theoretical balance. The model prohibits interpolating through 900–1000°C, where product phases require experimental resolution.
- Added 4-fold leave-one-temperature-out diagnostics on **four approximate 240-min endpoints from the same study**, making clear these do not constitute external validation. Each fold uses different held-out temperature.
- Added fail-closed independent reactor data ingestion requiring source DOI, experiment/ reactor IDs, same chemistry, feed carbon/LCO ratio, observed conversion and uncertainty. **No independently measured reactor time-series has been supplied.** Empty template intentionally fails, rather than filling synthetic observations.

## 2. OpenFOAM status

- The candidate v7 `reactiveFixedBedFoamV7` C++ solver and v9 `Allrun` remain included for OpenFOAM.com v2606.
- `v9/openfoam/Allrun` was **attempted locally**; it exits code 2 with `ERROR: wmake unavailable; runtime NOT VERIFIED`. Logs are retained.
- CI `.github/workflows/v9-verification.yml` includes a separate OpenFOAM setup/build/mesh/solve job using `gerlero/setup-openfoam@v1` with version 2606 and log upload on failure. This workflow has **not been run on the user's GitHub repository in this turn**; its job must pass before any OpenFOAM experience claim.
- OpenFOAM incompressible gas still uses *fixed* density, approximate species properties, a PISO-like solver, and simplified radiation source. Full compressible transport, gas solid particle kinematics, and plasma coupling remain **not verified**. No convergence log is currently available.

## 3. Study/provenance limitations

- LCO carbon-reduction kinetics article: Nuraeni et al., Metallurgical and Materials Transactions B 54 (2023), 602–620, DOI `10.1007/s11663-022-02712-1`. The study reported distinct apparent activation energies around 121 and 95 kJ/mol for diffusion/nucleation stages. v8/v9 stage prefactors are assumed/fitted only to approximate summary endpoints, not a complete measured time series.
- Independent chemistry context is available for **H₂ reduction of spent battery cathodes**, paper at https://pmc.ncbi.nlm.nih.gov/articles/PMC9548595/ (apparent Ea ~84.86 kJ/mol over 500–700°C), but that system is **not** directly comparable numerically to LiCoO₂ + C and was excluded from quantitative validation.
- NiO dense-plate hydrogen reduction is a **surrogate** and not a predictive battery black-mass rate law; interfacial active area and product porosity are not known.
- No reactor data are fabricated or relabelled as experimental observations.

## 4. What would close the remaining gaps

| Gap | Exit evidence |
|---|---|
| OpenFOAM execution | Verified CI or local `log.wmake_solver`, `log.blockMesh`, `log.checkMesh`, `log.solver`, completed run, divergence and budget diagnostics |
| OpenFOAM coupling | Stable mesh/time-step studies; detailed gas/solid species, enthalpy and flux balances on OpenFOAM outputs |
| LiCoO₂ reactor kinetics | Independent raw measured conversion-vs-time series with material composition, temperature and carbon ratio; fit/calibration/holdout and uncertainty |
| EM–CFD production loop | Actual CFD cell mesh mapping and thermal remapping on saved solver fields over multiple converged iterations; monitor total deposited power |
| Plasma physics | Electron energy equation, charge neutrality/Poisson, charged species transport, verified plasma benchmark and coupling to EM & reactor; current 0D chemistry does not meet this |
| Moving solids | Particle motion/collisions and interphase momentum transfer if reactor design actually requires it; a fixed bed does not necessarily require DEM |
| Reactor external validation | Real measured reactor temperature and product conversion from an independent apparatus; report MAE, RMSE, uncertainty and holdout splits |

**Safe CV description:** “Developed and tested Python/C++ reduced-order battery recycling reactor and EM–thermal models, including mass/energy conservation checks, kinetic scenario modelling, and uncertainty-aware same-study validation. Prepared (not yet run) an OpenFOAM v2606 fixed-bed reactor case.”
