# ReInduct-CFD Pro v10 — Research Evidence Summary

**Project type:** Independently developed numerical research portfolio for high-temperature battery-recycling process modeling. **Scope:** verification-scale solvers and literature comparisons; **not** an experimentally validated industrial digital twin.

## Computational contributions

- **Spatial plasma-transport prototype:** Implemented a one-dimensional, dimensionless periodic finite-volume transport system for H₂, H radicals, H₂⁺ ions, electrons, electron energy, and gas energy. Charge-dependent drift uses an electrostatic field computed from the discrete Poisson equation. Dissociation, ionization and recombination are mass- and charge-consistent. The energy audit includes stored electrostatic field energy and distinguishes work by the imposed electric field from internal redistribution. Compared with the former zero-dimensional chemistry model, this adds actual spatial transport and charge separation.
- **MHD–thermal benchmark:** Solved dimensional plane Hartmann flow using a finite-difference boundary-value formulation. The Lorentz braking force affects the velocity field; viscous and Joule dissipation provide the source for the conductive temperature profile. Verified flow against the analytical Hartmann velocity profile, mesh refinement, and mechanical power balance.
- **Battery recycling kinetics:** Transcribed all 10 stage-specific rate constants and reported regression R² values from Tables II and III of a 2023 LiCoO₂/C reduction study. Estimated activation energies of 121.17 and 95.28 kJ/mol, close to the paper's 121 and 95 kJ/mol. Evaluated leave-one-temperature-out performance for **reported rate coefficients from the same paper**, not independent reactor measurements.
- **Cross-study applicability:** Compared model-predicted product classes with independent published **product-recovery** findings in an air-atmosphere LiCoO₂/C study (2022). Carbonate appearance is consistent; cobalt oxidation state and operating atmosphere differ. This qualitative mismatch is reported as a model applicability limit, not hidden or converted into a numeric conversion error.
- **OpenFOAM v2606:** Prepared a non-destructive build/mesh/solve pipeline that requires actual compilation, computed cell geometry, conservative EM heat-source mapping, meshing diagnostics and parsed solver logs. The present environment has no `wmake` or `blockMesh`, so **OpenFOAM execution has not been verified**.

## Reproducibility and verification

The combined v4–v10 regression suite executes **108 Python tests** successfully in the current runtime, plus an independent C++17 kinetics check. The default 1D plasma model produces maximum inventory errors of approximately 1.3×10⁻¹⁵ (H atoms), 6.1×10⁻¹⁸ (net charge), and 2.4×10⁻¹² (total energy including field). For the 101-point Hartmann benchmark, the relative velocity L2 error is approximately 6.43×10⁻⁵. Outputs, raw arrays, input assumptions, logs and SHA-256 hashes are included under `v10/results/`.

## Key unsolved requirements

OpenFOAM source compilation and execution; physically calibrated plasma transport including electrodes, electron energy distribution and reaction cross-sections; full EM/plasma/gas–solid CFD coupling; independent LiCoO₂/C reaction-conversion/time data from a second reactor with comparable conditions. These tasks require a compatible scientific software environment and/or external measurements and **are not represented as completed work**.

## Suggested accurate CV phrasing

“Developed numerically verified gas–solid reactor and EM/thermal research prototypes in Python and C++, including a conservative 1D charged-species transport/Poisson model, a benchmarked Hartmann MHD–thermal solver, and published-data-based LiCoO₂/C kinetic analysis. Documented model transferability limits and built an OpenFOAM v2606 runtime-validation workflow (execution pending).”

## Primary references

1. Nuraeni et al., *Recovery of Cobalt and Lithium by Carbothermic Reduction of LiCoO₂ Cathode Material: A Kinetic Study*, *Metallurgical and Materials Transactions B* (2023). DOI: [10.1007/s11663-022-02712-1](https://doi.org/10.1007/s11663-022-02712-1).
2. She et al., *Product control and a study of the structural change process during the recycling of lithium-ion batteries based on the carbothermic reduction method*, *Journal of Chemical Research* (2022). DOI: [10.1177/17475198211066533](https://doi.org/10.1177/17475198211066533).
3. [OpenFOAM.com v2606](https://www.openfoam.com/current-release).
