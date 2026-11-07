"""Fail-closed numerical diagnostic of *actual* OpenFOAM solver log.

These are numerical sanity checks; they do NOT certify reactor physics.
"""
from pathlib import Path
import json,re,sys,math

def check(logdir):
    logdir=Path(logdir)
    f=logdir/'log.solver'
    if not f.exists():raise ValueError('No solver log generated')
    s=f.read_text()
    if not re.search(r'(?m)^End\s*$',s):raise ValueError('Solver did not report End')
    fields=('meanX','minH2','maxTs','totalEMPower','totalRadiativeLoss','maxDivergencePhi')
    rex=r'(?<![A-Za-z])([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?)'
    rows=[]
    for line in s.splitlines():
        if 'meanX=' not in line:continue
        vals={}
        for field in fields:
            match=re.search(r'\b'+field+r'=\s*'+rex,line)
            if match is None:raise ValueError('Missing '+field+' in runtime diagnostic')
            value=float(match.group(1))
            if not math.isfinite(value):raise ValueError('Nonfinite '+field)
            vals[field]=value
        if not 0<=vals['meanX']<=1:raise ValueError('Mean conversion outside [0,1]')
        if vals['minH2'] < -1e-7 or vals['maxTs']<=0 or vals['totalEMPower']<0:
            raise ValueError('Invalid concentrations, temperature or power')
        rows.append(vals)
    if len(rows)<2:raise ValueError('Need >=2 reported reactor timesteps')
    data={'status':'NUMERICAL_SANITY_ONLY','time_steps_checked':len(rows),
          'first':rows[0],'last':rows[-1],
          'limits':'No independent reactor experiment, no formal grid/time-step study, no plasma coupling confirmed'}
    (logdir/'numerical_sanity.json').write_text(json.dumps(data,indent=2)+'\n')
    return data
if __name__=='__main__':
    print(json.dumps(check(sys.argv[1]),indent=2))
