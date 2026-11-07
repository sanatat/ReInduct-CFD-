# v7 verified numerical and candidate OpenFOAM port

See `../README_V7_FA.md` and `docs/VERIFICATION_V7.md`.

`bash v7/run_v7.sh` runs Python verification, independent C++ smoke check, source mapping, holdout metrics and phase-envelope checks. `bash v7/openfoam/Allrun` is intentionally fail-closed unless OpenFOAM.com v2606 actually exists and runs. No source or output in this archive claims the OpenFOAM solver was successfully compiled locally.

Physics retained: v6 axisymmetric EM and fixed-bed surrogate. New: table-based surface growth -> bulk rate (explicit area assumption), consistent molar species source units for an OpenFOAM candidate, conservative cell heat-source map with geometric export, 800C held-out NiO kinetics, battery qualitative phase scenarios.
