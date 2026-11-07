# OpenFOAM.com v2606 actual-runtime verification

Run `bash v10/openfoam/Allrun` **after** sourcing OpenFOAM.com v2606.
Only `wmake`, `blockMesh`, `checkMesh`, geometry export, power mapping, the
real `reactiveFixedBedFoamV7` executable, and parsed solver logs can produce
`v10/results/openfoam_runtime/VERIFIED.txt`.

- This script copies the OpenFOAM case before running to protect the source.
- `OPENFOAM_NOT_VERIFIED` on missing executables is **not** a successful run.
- It uses the pre-existing fixed-bed *NiO/H2* surrogate, not a validated LiCoO2 recycling solver.
- A completed solver log passes a numerical sanity screen, not independent
  physical or experimental validation. Real OpenFOAM residual/budget and
  grid/time-step verification should be added before scientific claims.

### Run with OpenFOAM v2606 active on Linux Mint

```bash
source /opt/openfoam2606/etc/bashrc  # change to your actual installed path
bash v10/openfoam/Allrun
```

### Optional community container (requires Docker + Internet)

A Docker image tagged v2606 is published by microfluidica/openfoam; check the
image environment provides `wmake` and run this exact script:

```bash
docker run --rm -v "$PWD:/work" -w /work microfluidica/openfoam:2606 \
  bash -lc 'bash v10/openfoam/Allrun'
```

The corresponding `.github/workflows/v10-verification.yml` likewise runs
OpenFOAM in CI and uploads evidence on both success and failure. **No CI job has
been executed in the user's account by creating these files.**
