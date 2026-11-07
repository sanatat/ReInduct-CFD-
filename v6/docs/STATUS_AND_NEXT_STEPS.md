# v6 research integrity / technical status (2026-11-05)

| Required competency | Demonstrated in locally run model | Remaining limitations |
|---|---|---|
| C++ programming | Compiled kinetics kernel + numerical equivalence Python; legacy C++17 reactor | No *compiled* OpenFOAM custom extension |
| Reactive transport / CFD | FV axisymmetric diffusion/advection/reaction + v4 periodic-channel incompressible MAC CFD | v6 axisymmetric reactor uses pressure-driven plug/Ergun closure; full compressible reactive NS not delivered |
| Gas-solid multiphase | Distinct T_g and T_s, reaction and solid conversion, transport and Ergun closure | No resolved phase motion, granular stress, collisions, Eulerian-Eulerian or CFD-DEM model |
| Electromagnetic induction | 2D axisymmetric azimuthal current, Joule heating, A_theta, magnetic field and thermally dependent sigma | Coil is Gaussian current band; not CAD winding, measured magnetic properties or nonlinear 3D induction |
| Magnetohydrodynamics/plasma | Legacy Hartmann verification and illustrative motional force from v5 | No self-consistent plasma ionization/charge, electron energy, radiation or plasma MHD coupling |
| Thermodynamics/kinetics | NiO/H2 equilibrium estimates 500-1000K; group-heldout kinetic comparisons | NiO layer growth data not converted/calibrated to volumetric bed reaction; T ranges limited |
| Battery black mass | LiCoO2/C stoichiometric flowsheet and carbon/separation scenario | No measured black-mass process data, phase equilibria or lithium/cobalt real yields |
| Verification / reproducibility | Manufactured solution, spatial convergence, discrete balance and 43 passing tests | Not experimental validation or independent multi-lab replication |
| OpenFOAM | Case static preflight; local/Docker runner; unexecuted GitHub CI | Need first successful wmake, checkMesh and actual solver logs, likely code debugging |

## Prioritized next experiments

1. Publish repository, run GitHub Action `v6-validation.yml`, fix any OpenFOAM v2606 compilation/API/runtime failure based on real logs. **Until then, claim only static preflight.**
2. Convert axisymmetric reaction case into an axisymmetric wedge OpenFOAM mesh (or a validated axisymmetric solver) and verify pressure/continuity/energy, not simply re-label Cartesian heat mapping.
3. Replace NiO surrogate with a published chemically validated LIB cathode/black-mass reduction and actual property/kinetic dataset; validate temperature and species outputs against experimental curves.
4. Use a coupled granular or Eulerian-Eulerian gas-solid momentum model, then test grid/particle independence.
5. Add plasma only after solving a charge/electron-energy/conductivity model with physically justified ionization data and dimensional consistency. Proxies should not be called discharge simulations.

### Caveat on figures and numerical residuals

Mass/energy residuals are internal solver accounting checks. The absolute errors are expected to be tiny for the finite-volume scheme; they do not measure real-world accuracy. Peak temperature and conversion vary with assumed electrical conductivity and reaction law. They must not be presented as experimental predictions in an academic application.
