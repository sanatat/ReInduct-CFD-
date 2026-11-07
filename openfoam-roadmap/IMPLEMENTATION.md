# OpenFOAM development tasks — NOT implemented in this repository

Use a single, precisely recorded distribution and version for each OpenFOAM case. These are milestones, not accomplished features.

## Stage 1 — verified flow and passive scalar
- Create a simple 2D duct with a thin, stationary porous resistance zone; configure velocity inlet and pressure outlet.
- Use `blockMesh`, `checkMesh`, `foamRun` with a selected flow solver module on Foundation versions >=11; adapt for actual locally installed version.
- Prescribe passive scalar at inlet and examine advection. Report average residence time, inlet/outlet mass-flow closure, 3-grid comparison.
- Do not turn on chemistry, heating or multiphase physics yet.

## Stage 2 — reactive species
- Add H2 and H2O transport and an explicitly stoichiometric NiO reaction sink/source.
- Reuse transparent simple rates from C++ pilot; verify positivity and outlet species flux and X budget.
- Compare 1D centerline/bulk results with C++ model at identical temperature, porosity and flow assumptions.
- Document the precise solver class, data format, and boundary dictionaries for the installed OpenFOAM release.

## Stage 3 — separate gas and solid regions
- Build gas and solid meshes sharing a consistent interface and use conjugate heat transfer (CHT) with `foamMultiRun` / appropriate modules on Foundation OpenFOAM.
- Benchmark steady pure conduction across two layers using analytic thermal resistance, BEFORE reactive terms.
- If physically modeling gas permeation *inside* a solid bed, a volume-averaged two-temperature porous model is more appropriate than a solid block obstructing gas flow; do not confuse the two.
- Add interphase heat transfer and reaction source once thermal verification passes.

## Stage 4 — electromagnetics
- Define induction frequency, current, conductor geometry, field boundary conditions and temperature-dependent conductivity.
- Solve low-frequency Maxwell (Eddy-current) equations using a compatible library/solver; mapping q_Joule = 0.5 Re(J·E*) or equivalent to energy.
- Verify skin depth and integral power in a known slab/cylinder benchmark. Then map to OpenFOAM with energy conservation (total EM losses ≈ total volumetric heat input).
- Use Q prescribed as an exploratory alternative, labeled as such; **never call that a solved EM field**.

## Stage 5 — MHD/plasma only after EM verification
- Add J×B to momentum equation only with justified conductivity, field and plasma regime.
- Add temperature-dependent ionization and radical kinetics using defensible literature data; radicals are not equivalent to molecular hydrogen H2.
- Report controls with/without MHD/plasma, chemical mass balance, species positivity, power balance.

## Stage 6 — uncertainty and true validation
- Fit kinetics and physical parameters to published experiments; avoid leaking fitted values into independent validation data.
- Report model sensitivity to mesh, time step, inlet temperature, gas speed, and heating power.
- Evaluate actual black-mass multicomponent recovery only after parameterization and experimental comparison.

## Suggested outputs for each successful step
1. Case dictionaries and custom C++ source with reproducible shell instructions.
2. Screenshot of geometry, velocity vectors, T field, conversion contours, but only if actually computed.
3. `postProcessing/` and CSV tables for mass/energy budgets and grid study.
4. `RESULTS.md`: run date, CPU/RAM, distribution version, exact command, parameters, results and limitations.

## Linux Mint and resource notes
Linux Mint 22.x is Ubuntu 24.04-based, but binary installation compatibility must be checked against repository and runtime dependencies. A supported Ubuntu container or directly supported binary is usually simpler than compiling a large CFD package on <8 GB RAM. OpenFOAM Foundation and OpenCFD forks differ in naming and syntax; use their matching official docs. Keep early meshes around 5k–20k cells, single process before trying parallel execution. Do not assume that a complex plasma run is feasible on this hardware.
