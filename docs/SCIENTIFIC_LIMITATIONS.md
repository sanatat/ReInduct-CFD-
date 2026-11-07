# Scientific scope and unresolved gaps

The active model is a 1D porous-bed gas/solid reduced-order numerical surrogate with a separate axisymmetric induction-heating solver. It should not be called a full multiphase CFD, OpenFOAM, plasma/MHD, or validated industrial battery-recycling model.

| Topic | Present capability | Evidence still needed |
|---|---|---|
| OpenFOAM | Historical code and cases only | Port active v13 physics, compile with `wmake`, check mesh, run, demonstrate convergence and conservation |
| Reactive CFD | Reduced-order molar transport and Darcy–Ergun pressure closure | Full PIMPLE/Navier–Stokes with compressible mixed-gas transport |
| Multiphase & multiregion | Gas/solid nonequilibrium temperatures in a porous-bed surrogate | Independent phases/meshes or validated solid mobility and interphase coupling |
| Plasma/MHD | Illustrative H₂/H dissociation and separate Hartmann/plasma prototypes | Calibrated electron/ion kinetic chemistry, physical boundary conditions and a coupled reactor model |
| Electromagnetics | Axisymmetric induction model with refreshed thermal feedback | Full simultaneous EM–CFD coupling and grid/time-step study of mapped sources |
| Solid thermochemistry | Temperature-dependent *gas sensible* Shomate properties | Solid phase heat capacities and formation enthalpies; complete reaction ΔH(T) |
| Battery-recycling kinetics | Literature-inspired effective LiCoO₂/C rate proxies | Independently calibrated kinetic law for realistic black-mass composition and surface area |
| Reactor validation | Code-level tests, benchmarks and preliminary sensitivity | Independent matched experiment with temperature, conversion, outlet gas and uncertainty |

## Key references

- NIST CO: https://webbook.nist.gov/cgi/cbook.cgi?ID=C630080&Mask=1
- NIST CO₂: https://webbook.nist.gov/cgi/cbook.cgi?Name=carbon+dioxide&cTG=on
- NIST N₂: https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Plot=on&Type=JANAFG
- NIST H₂: https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=C0F
- Nuraeni et al. (2023): https://pmc.ncbi.nlm.nih.gov/articles/PMC9821372/

## Safe description for applications

> Developed and numerically tested a reduced-order reactive gas–solid battery-material recycling model with pressure–enthalpy iteration, temperature-dependent sensible-gas properties, induction-heating feedback and Python/C++17 regression tests. OpenFOAM operation and experimental reactor validation are not yet demonstrated.
