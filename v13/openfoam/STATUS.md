# OpenFOAM verification boundary

The v13 numerical core is a Python FV/Ergun reduced-order model with C++17 thermochemical cross-check. It is **not an OpenFOAM executable**. The legacy v11 `reactiveLCOFoamV11` candidate distributed for historical continuity is not physically equivalent to v13 and is **not compiled or runtime-verified**. The `v13/openfoam/Allrun` script fails closed instead of reporting success.

The historical case may be tested independently using `bash v12/openfoam/Allrun` on an OpenFOAM.com v2606 environment; success would be a smoke test of the historical candidate **only**. A true v13 OpenFOAM implementation requires porting temperature-dependent species enthalpies, the pressure/energy Picard update and EM refresh feedback, plus performing solver/mesh convergence and physical validation.
