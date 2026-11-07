# v7 OpenFOAM 2606 candidate (NOT COMPILED in the supplied environment)

`bash v7/openfoam/Allrun` does all steps **on a machine with OpenFOAM.com v2606**:
1. `wmake` C++ solver and mesh geometry exporter.
2. `blockMesh`/`checkMesh` for annular axisymmetric 4-degree wedge, 12 radial by 48 axial cells.
3. `exportCellGeometry` writes actual per-cell (x,y,z,V) in original mesh order.
4. `v7.mapping` solves v6 EM model, conservatively maps heating into OpenFOAM cells, writes `0/qEM` with SI units.
5. Executes the C++ reactive fixed-bed solver, captures residuals and checks that the process really ends.

The `Allrun` script **fails** if any dependency is missing or a step fails. A passing Python test **cannot** replace this.

### Physical assumptions
* Constant-density incompressible superficial velocity and Ergun drag; no particle translation or plasma.
* Gas species H2 and H2O as **mol per m3 of pore gas** (not mass fractions); porosity 0.45.
* Interstitial inventory epsilon*c [mol per m3 bulk]; NiO solid capacity 450 mol per m3 bulk.
* Reaction R [mol per m3 bulk per s], H2 source -R/epsilon, H2O +R/epsilon, X rate +R/capacity.
* Joule heating qEM [W/m3 bulk], temperature K; reaction heat -deltaH*R. Thermophysical constants are illustrative.
* Numerical feedback loops are Picard-type within timesteps. This is NOT a proven fully converged transient or validated digital twin.

### Known risk
The solver is a **candidate**, written against the v2606 API but uncompiled here. Mesh wedge face orientation, OpenFOAM operator overloads and field I/O may require corrections in GitHub CI. Never mark this runtime-complete without actual logs.


## v8 corrections / limitations
- Corrected **Ni product-layer molar density** in the interface-to-volume source.
- H2 mole fraction now computed from **gas temperature Tg**, not solid Ts.
- Added lumped grey-wall radiation loss qRad [W/m3 bulk].
- Tentative pressure-velocity flux projection revised and log of continuity divergence added.
- **OpenFOAM has NOT been compiled/executed here**; these edits are source candidates.
- Gas density remains constant for this incompressible PISO candidate; variable-density high-temperature
  flow still requires a consistent compressible solver. The v8 Python mixture properties are standalone.
- The thermal EM feedback test is run separately (see v8/thermal_em.py); it is not yet an
  in-memory bidirectional EM/OpenFOAM co-simulation.
