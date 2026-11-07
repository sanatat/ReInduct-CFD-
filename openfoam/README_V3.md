# OpenFOAM adapter (not executed on this system)

The package contains two case skeletons:

* `cavity2D` — laminar lid-driven cavity with `icoFoam`.
* `cavityScalar2D` — same cavity plus **passive** H2 transport via a `scalarTransport` function object (OpenCFD / openfoam.com syntax). It is **not** a reacting gas or multiregion case.

OpenCFD docs: https://doc.openfoam.com/2606/tools/post-processing/function-objects/solvers/scalarTransport/

Run after sourcing a compatible OpenFOAM environment:

```bash
bash openfoam/run_cavity_scalar.sh
```

`icoFoam` and `scalarTransport` may be absent or configured differently in OpenFOAM Foundation releases. Do not interpret a successful `icoFoam` run as validated reactive CFD: the external 2D Python solver is the only executable reacting + EM/T coupling in this package.

**The remaining OpenFOAM milestone:** port the model's reaction stoichiometry and gas/solid energy equations to an appropriate solver/fvModel implementation; validate against a simpler benchmark, then compare against this package. This code does not yet provide that port.
