"""Fail-closed OpenFOAM build+runtime evidence inspector. Never inferred from CI config."""
from pathlib import Path
import sys,re,json

def check(path):
    p=Path(path)
    required=('wmake_solver','wmake_geometry','blockMesh','checkMesh','geometry','mapping','solver')
    missing=[s for s in required if not (p/f'log.{s}').is_file() or (p/f'log.{s}').stat().st_size==0]
    if missing:raise RuntimeError(f'Missing OpenFOAM runtime logs: {missing}')
    logs={s:(p/f'log.{s}').read_text(errors='replace') for s in required}
    if not re.search(r'Mesh OK\.',logs['checkMesh']):raise RuntimeError('checkMesh did not report Mesh OK')
    if not re.search(r'(?m)^End\s*$',logs['solver']):raise RuntimeError('solver did not complete')
    times=re.findall(r'Time\s*=\s*([0-9.]+)',logs['solver'])
    if len(times)<2:raise RuntimeError('fewer than 2 solver timesteps')
    if re.search(r'FOAM FATAL|floating point exception|segmentation fault',logs['solver'],re.I):
        raise RuntimeError('solver fatal failure')
    rows=re.findall(r'V11 meanLCO=([\d.eE+-]+) maxTs=([\d.eE+-]+) maxDivPhi=([\d.eE+-]+) EMpower=([\d.eE+-]+)',logs['solver'])
    if len(rows)<2:raise RuntimeError('missing reactor diagnostics')
    import math
    vals=[[float(q) for q in row] for row in rows]
    if any(not all(math.isfinite(x) for x in v) for v in vals):raise RuntimeError('nonfinite output')
    if any(v[0]<0 or v[1]<200 or v[2]>1.0 or v[3]<=0 for v in vals):
        raise RuntimeError('failed physical sanity checks (including loose divergence threshold)')
    if not (p/'work/case/0/qEM').exists():raise RuntimeError('qEM not mapped')
    out={'openfoam_runtime_verified':True,'timesteps':len(times),'reactor_output_records':len(rows),
         'max_divergence_1_s':max(x[2] for x in vals),'EM_power_W':vals[-1][3],
         'note':'Numerical runtime evidence ONLY. Not reactor experimental validation.'}
    (p/'runtime_evidence.json').write_text(json.dumps(out,indent=2)+'\n')
    return out

if __name__=='__main__':
    try:print(json.dumps(check(sys.argv[1]),indent=2))
    except Exception as ex:
        print(f'NOT_VERIFIED: {ex}',file=sys.stderr);sys.exit(2)
