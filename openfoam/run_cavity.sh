#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/cavity2D"
command -v blockMesh >/dev/null || { echo 'ERROR: OpenFOAM environment not sourced; blockMesh unavailable' >&2; exit 2; }
command -v icoFoam >/dev/null || { echo 'ERROR: icoFoam not available in this distribution; try an installed legacy/ESI release or adapt to foamRun.' >&2; exit 2; }
blockMesh > log.blockMesh 2>&1
if command -v checkMesh >/dev/null; then checkMesh > log.checkMesh 2>&1; fi
icoFoam > log.icoFoam 2>&1
printf 'Case finished; see %s/log.icoFoam\n' "$PWD"
