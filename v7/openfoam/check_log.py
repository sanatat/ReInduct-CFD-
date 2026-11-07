"""Fail-closed assessment of real OpenFOAM build/mesh/solver logs."""
from __future__ import annotations
from pathlib import Path
import json,re,sys

def check(logdir,case):
    logdir=Path(logdir);case=Path(case)
    files={k:(logdir/f'log.{k}') for k in ('wmake_solver','wmake_geometry','blockMesh','checkMesh','geometry','mapping','solver')}
    exist={k:p.exists() and p.stat().st_size>0 for k,p in files.items()}
    lines=files['solver'].read_text() if exist['solver'] else ''
    mesh=files['checkMesh'].read_text() if exist['checkMesh'] else ''
    check={'logs_exist':all(exist.values()),'solver_finished':bool(re.search(r'(?m)^End\s*$',lines)),
           'mesh_success':bool(re.search(r'Mesh OK\.',mesh)),
           'timestep_count':len(re.findall(r'Time\s*=\s*\S+',lines)),
           'coupling_iterations_reported':bool(re.search(r'coupling\s+\d+',lines)),
           'solver_no_fatal':not bool(re.search(r'FOAM FATAL|Floating point exception|Segmentation fault',lines,re.I)),
           'qEM_field_written':(case/'0/qEM').is_file()}
    check['runtime_verified']=all([check['logs_exist'],check['solver_finished'],check['mesh_success'],check['timestep_count']>0,check['coupling_iterations_reported'],check['solver_no_fatal'],check['qEM_field_written']])
    (logdir/'runtime_evidence.json').write_text(json.dumps(check,indent=2)+'\n')
    if not check['runtime_verified']:raise SystemExit('OpenFOAM runtime evidence INSUFFICIENT')
    return check
if __name__=='__main__':print(json.dumps(check(sys.argv[1],sys.argv[2]),indent=2))
