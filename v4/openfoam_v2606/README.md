# OpenFOAM v2606 custom solver: candidate port (UNCOMPILED)

`reactiveFixedBedFoam.C` is source-code for the **OpenCFD/OpenFOAM.com v2606**
API family. It is an *experimental port* of a subset of the separate Python
multiphysics model. The environment used to create v4 has **no OpenFOAM installation**.
It was therefore **not compiled or executed**, and syntax/API adjustments may be needed.
It must not be called an OpenFOAM CFD result until run and checked.

After installing the matching version and `source /path/to/OpenFOAM.../etc/bashrc`:

```bash
bash v4/openfoam_v2606/Allrun
```

Expected successful stages: `wmake`, `blockMesh`, `checkMesh`, `reactiveFixedBedFoam`.
The logs appear in `v4/openfoam_v2606/log.wmake` and the case folder.
Fix runtime and physics issues based on actual logs; mesh, Courant,
reactive source bounding, and conservation must be tested before publication.

## Included physics

- Incompressible laminar PISO for velocity / kinematic pressure
- Simplified localized gas-solid momentum sink depending on Ts and X
- NiO + H2 -> Ni + H2O H2/H2O and solid conversion X
- Simplified two-temperature LTNE exchange
- Imported 2-D Python induction **qEM** scalar heat source (one-way/offline coupling only)

## Known limitations

- Not a granular Euler–Euler model, no MHD, no plasma.
- Uses normalized gas H2 and H2O variables, not rigorously computed thermodynamics.
- Hard-coded illustrative properties.
- Needs compatibility checking with exact OpenCFD v2606 headers and libraries.
- Source terms need quantitative conservation audit against transported fluxes.
- Constant-density flow is not suitable for industrial H2 temperature changes.

Official documentation: https://doc.openfoam.com/2606/quickstart/ and
https://doc.openfoam.com/2606/tools/processing/numerics/fvoptions/

## V5 induction heat mapping

Run `python3 -m v5.openfoam.export_foam_heat` to regenerate `channelReactive2D/0/qEM` for the 30×12×1 mesh. qEM carries units [W/m³] and is divided by `rhoCpS` in the solid temperature equation. This is a one-way static heat source exported before the run, **not** iterative EM/CFD feedback inside OpenFOAM. The C++ source has not been compiled in the offline environment.
