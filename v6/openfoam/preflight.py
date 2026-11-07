"""Offline structure and dataset preflight; does NOT compile or execute OpenFOAM."""
from pathlib import Path
import re,json,shutil

ROOT=Path(__file__).resolve().parents[2]
CASE=ROOT/'v4/openfoam_v2606/channelReactive2D'
SOLVER=ROOT/'v4/openfoam_v2606/reactiveFixedBedFoam'

def inspect():
    fields=['U','p','H2','H2O','Tg','Ts','X','qEM']
    existing={name:(CASE/'0'/name).is_file() for name in fields}
    bmd=(CASE/'system/blockMeshDict').read_text()
    mesh=re.search(r'hex\s*\([^)]*\)\s*\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)',bmd)
    dims=[int(k) for k in mesh.groups()] if mesh else None
    qtext=(CASE/'0/qEM').read_text()
    nums=re.search(r'internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\);',qtext,re.S)
    count=int(nums.group(1)) if nums else None
    vals=[float(k) for k in nums.group(2).split()] if nums else []
    qty_match=bool(nums and len(vals)==count and dims and count==dims[0]*dims[1]*dims[2])
    report={'status':'static files present; runtime not attempted',
            'environment_wmake':bool(shutil.which('wmake')),
            'environment_blockMesh':bool(shutil.which('blockMesh')),
            'environment_checkMesh':bool(shutil.which('checkMesh')),
            'expected_fields':existing,'mesh_shape_xyz':dims,'qEM_list_count':count,
            'qEM_cell_count_verified':qty_match,'qEM_nonnegative':bool(vals and min(vals)>=0),
            'qEM_dimensions_W_m3':bool(re.search(r'dimensions\s*\[\s*1\s+-1\s+-3\s+0\s+0\s+0\s+0\s*\]',qtext)),
            'solver_source_present':(SOLVER/'reactiveFixedBedFoam.C').is_file(),
            'solver_make_present':(SOLVER/'Make/options').is_file() and (SOLVER/'Make/files').is_file(),
            'runScript_present':(ROOT/'v4/openfoam_v2606/Allrun').is_file(),
            'openfoam_compiled':False,'openfoam_solver_run':False}
    report['static_preflight_passed']=bool(all(existing.values()) and qty_match and
                report['qEM_nonnegative'] and report['qEM_dimensions_W_m3'] and report['solver_source_present'] and
                report['solver_make_present'] and report['runScript_present'])
    return report

if __name__=='__main__':
    report=inspect()
    dest=ROOT/'v6/results/openfoam_static_audit.json'
    dest.parent.mkdir(exist_ok=True,parents=True)
    dest.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if not report['static_preflight_passed']:raise SystemExit(1)
