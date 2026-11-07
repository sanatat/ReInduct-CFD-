# OpenFOAM case status

`openfoam/cavity2D` is a full legacy/ESI-style **non-reacting incompressible 2-D lid-driven cavity** case. It demonstrates OpenFOAM dictionaries, mesh, boundary conditions and PISO momentum/pressure solution. It is NOT a model of the reactor, and it has NOT been executed in this container because OpenFOAM is not installed. It is supplied as a first reproducible OpenFOAM exercise and may require adaptation to your distribution.

Run after sourcing an installed OpenFOAM environment:

```bash
make openfoam
```

For Foundation distributions that lack `icoFoam`, use the installed modern tutorial and its supported `foamRun -solver ...` procedure. **Do not claim a successful OpenFOAM run without a log.**

For reactive CFD, prefer starting from a supported reacting solver tutorial for the exact installed OpenFOAM release, rather than assuming configuration syntax is stable across releases. Develop source terms and region coupling incrementally. Track solver version and model fidelity carefully.
