"""Honest OpenFOAM availability report, NOT a compile/runtime replacement."""
from __future__ import annotations
import json,os,shutil,subprocess
from pathlib import Path

def report(out='v9/results/openfoam_environment.json'):
    required=['wmake','blockMesh','checkMesh']
    found={exe:shutil.which(exe) for exe in required}
    enabled=all(found.values())
    r={'WM_PROJECT_VERSION':os.environ.get('WM_PROJECT_VERSION'),
       'executables':found,'prerequisites_present':enabled, 'custom_solver_requires_local_build':True,
       'solver_compiled_and_executed':False,
       'next_step':'bash v9/openfoam/Allrun (only where OpenFOAM v2606 is installed)',
       'warning':'preflight does not certify OpenFOAM API compatibility or runtime; success requires preserved logs from full Allrun'}
    p=Path(out);p.parent.mkdir(exist_ok=True,parents=True);p.write_text(json.dumps(r,indent=2)+'\n')
    return r
if __name__=='__main__':print(json.dumps(report(),indent=2))
