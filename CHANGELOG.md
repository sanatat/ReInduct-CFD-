# Changelog

## v13.1 — GitHub publication edition (2026-11-07)

**Publication and reproducibility changes only: the numerical physics is still v13.**

- Retained the functioning v13 research core and historically required numerical modules.
- Removed cached bytecode, redundant historical results, old OpenFOAM CI claims and superseded workflows from the publication snapshot.
- Added a single `make verify` command and fail-closed publication audit.
- Added one GitHub Actions workflow that explicitly tests **numerical code, not OpenFOAM**.
- Added source-provenance documentation, limitations, a safe citation template, contribution instructions and a security/safety disclaimer.
- Reproduced 153 passing Python tests and a successful C++17 thermochemistry cross-check in the preparation environment.
- Noted OpenFOAM runtime, plasma–MHD and independent reactor experimental validation as outstanding.

## v13.0 — Numerical research core

Pressure–enthalpy Picard iteration, NIST Shomate sensible-gas properties, EM/thermal refresh feedback, C++17 regression and parameter screening. Prior histories are retained under `docs/history/`.