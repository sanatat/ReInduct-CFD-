#!/usr/bin/env bash
# Passive scalar OpenCFD branch tutorial scaffold. UNEXECUTED in this package.
set -euo pipefail
cd "$(dirname "$0")/cavityScalar2D"
for bin in blockMesh icoFoam; do
  command -v "$bin" >/dev/null || { echo "Missing $bin; source a compatible OpenFOAM installation" >&2; exit 2; }
done
blockMesh > log.blockMesh 2>&1
checkMesh > log.checkMesh 2>&1 || { cat log.checkMesh; exit 3; }
icoFoam > log.icoFoam 2>&1
if ! grep -q 'scalarTransport' log.icoFoam; then echo 'WARNING: confirm that scalarTransport executed (see log.icoFoam)' >&2; fi
