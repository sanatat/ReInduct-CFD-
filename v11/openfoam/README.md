# OpenFOAM.com v2606 candidate — NOT COMPILED

`reactiveLCOFoamV11.C` is an *unverified* custom solver intended to replace the older NiO-specific V7 candidate. It includes a three-reaction LCO/C scheme, dimensioned gas/solid inventories, a PIMPLE outer loop with pressure-correction calls, and multi-field outputs. The model is deliberately **incompressible, fixed-bed and constant-property**; this means it is not a complete gas mass/momentum-conserving reactive industrial solver. No plasma terms are present.

Only `bash v11/openfoam/Allrun` in a **genuine OpenFOAM.com v2606** environment can establish compilation/runtime compatibility. The script fails if any dependency is missing, writes its logs to `v11/results/openfoam_runtime`, uses a sandboxed case copy, exports mesh-derived cell centroids, maps qEM with bulk power conservation, and runs the candidate. It does **not** claim runtime verification in this ZIP.

The standalone v11 axial FV model in Python is the module with actual gas/solid element and energy conservation tests and EM temperature feedback. The OpenFOAM solver is a prospective integration target; independently verified results **must not** be attributed to it.

* Kinetic stage constants: Nuraeni et al., DOI `10.1007/s11663-022-02712-1`, Tables II and III; effective rates are NOT independently calibrated local kinetics.
* PIMPLE API reference: https://api.openfoam.com/2606/classFoam_1_1pimpleControl.html
* The v11 OpenFOAM branch still needs compilation, divergence assessments, true state-dependent properties, and experimental verification.
