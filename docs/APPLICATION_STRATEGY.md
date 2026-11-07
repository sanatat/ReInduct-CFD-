# Application strategy for TU-BRIDGE DP.MUL.3

**Key caution:** A reduced-order C++ transport solver is a meaningful first step but does not meet the advertised hands-on reactive CFD/OpenFOAM experience by itself. This project cannot erase that gap merely by being listed on a CV.

## Four working sessions, prioritized for 8 Nov 2026 submission

**Session 1 (prototype, achieved):** Run and check 1-D C++ reactive-bed model; preserve grid and time-step study. Capture results, limitations and reproducible scripts. Verify inputs in README and save repo to GitHub.

**Session 2 (OpenFOAM basics, pending):** Install a supported distribution, run a standard flow tutorial unmodified and then with one modified boundary condition, record exact version and settings. Begin a simple 2D duct case with 5k–20k cells. Document `blockMesh`, `checkMesh`, convergence, velocity fields and outlet flux closure.

**Session 3 (reactive scalar, pending):** Introduce a simple advected scalar or a documented reactive-flow tutorial. Show reaction on/off controls, species balance, and conservation checks. Match first-order kinetics to reduced-order reference only when units align.

**Session 4 (thermal + literature, pending):** Use a supported gas/solid heat-transfer tutorial and verify analytic two-layer thermal resistance. Compare relevant published NiO kinetics, explicitly state what is and is not reusable for black mass/plasma.

**Submission rule:** Clearly show exactly which stages actually ran. If OpenFOAM stage cannot run before the deadline, present the C++ verification project and describe subsequent OpenFOAM work as a proposed research direction, not prior experience.

## Honest CV bullet for the present deliverable

> Developed and tested a C++17 reduced-order finite-volume benchmark of reactive H2–NiO transport in a porous heated bed, featuring two-temperature gas/solid coupling, Arrhenius-type reaction kinetics, hydrogen conservation checks, and mesh/time-step sensitivity studies.

## Stronger CV wording only after independently running OpenFOAM stages

> Implemented and benchmarked an OpenFOAM reactive-flow case with documented species conservation and mesh convergence, informed by a separately verified C++ reactive porous-bed model.

Do not use the stronger bullet until the OpenFOAM case and verification actually exist.

## Relevance to requested skill gaps

| Skill | Present | Further proof needed |
|---|---|---|
| Reactive CFD | Reaction and transport in 1-D C++ *surrogate* | Real momentum/energy/species CFD solver |
| Multiphase / multiregion | Two continuum temperatures in same 1D bed | Coupled regions or explicit two-phase equations |
| High temperature | Temperature-dependent Arrhenius reaction and thermal model | Radiative/high-temp reactor validation |
| Induction/plasma | Prescribed localized heat and textbook skin depth | Maxwell field, Joule source, plasma physics |
| Electromagnetics / MHD | Scalar skin-depth equation | E, B, J fields; JxB forces if physically relevant |
| C++ | Working compiled finite-volume code | Integration into OpenFOAM C++ objects |
| OpenFOAM | Learning roadmap | Actual installed, executed, documented case |
| Thermodynamics / kinetics | Arrhenius kinetics + illustrative reaction enthalpy | Published parameter fits, gas/species properties, reaction network |
