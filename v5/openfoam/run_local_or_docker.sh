#!/usr/bin/env bash
# Launches an actual OpenCFD v2606 install or reproducible community Docker image.
# No fake success: exits nonzero when compilation, meshing or solver fail.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
if command -v blockMesh >/dev/null 2>&1 && command -v wmake >/dev/null 2>&1; then
  echo 'OpenFOAM native tools found; building and running the v2606 solver.'
  bash v4/openfoam_v2606/Allrun
elif command -v docker >/dev/null 2>&1; then
  echo 'OpenFOAM native tools not found; attempting v2606 Docker build.'
  echo 'Container image: microfluidica/openfoam:2606 (community, not Foundation official).'
  docker run --rm --entrypoint bash -v "$ROOT:/work" -w /work microfluidica/openfoam:2606 -lc '
    set -euo pipefail
    if ! command -v wmake >/dev/null 2>&1; then
      for f in /opt/openfoam*/etc/bashrc /opt/OpenFOAM/*/etc/bashrc /usr/lib/openfoam/*/etc/bashrc; do
        if [[ -f "$f" ]]; then source "$f"; break; fi
      done
    fi
    command -v blockMesh
    command -v wmake
    bash v4/openfoam_v2606/Allrun
  '
else
  echo 'OPENFOAM_BLOCKED: neither a sourced native OpenFOAM v2606 with wmake/blockMesh nor Docker is installed.' >&2
  echo 'Install Docker on Linux Mint or source OpenFOAM v2606 and rerun: make v5openfoam' >&2
  exit 2
fi
