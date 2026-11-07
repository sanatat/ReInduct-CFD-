"""Fail-closed static checks of v11 OpenFOAM candidate (NOT compilation)."""
from pathlib import Path
import json,re

def lint(root=None):
    if root is None:
        root=Path(__file__).resolve().parents[2]
    root=Path(root)
    case=root/'v11/openfoam/case'
    solver=(root/'v11/openfoam/solver/reactiveLCOFoamV11.C').read_text()
    schemes=(case/'system/fvSchemes').read_text()
    solution=(case/'system/fvSolution').read_text()
    errors=[]
    # Only fields occurring inside fv scheme's div(phi,F) actually need a scheme.
    fields=re.findall(r'fvm::div\(phi\s*,\s*(\w+)\)',solver)
    for fld in set(fields):
        if not re.search(r'div\(phi,\s*'+re.escape(fld)+r'\)\s+Gauss\s+\w+\s*;',schemes):
            errors.append('missing scheme div(phi,'+fld+')')
    for name in ['CO','CO2','N2','Tg','Ts','LCO','Carbon','Li2O','CoO','CoMetal','Li2CO3','qEM','p','U']:
        path=case/'0'/name
        if name=='qEM' and not path.exists() and (case/'0/qEM.template').exists():
            path=case/'0/qEM.template'
        if not path.is_file():errors.append('missing initial field '+name)
        elif not re.search(r'object\s+'+name+r'\s*;',path.read_text()):errors.append('wrong field object '+name)
    for name in ['nCorrectors','nOuterCorrectors']:
        if not re.search(r'\b'+name+r'\s+[1-9]',solution):errors.append('missing PIMPLE '+name)
    if 'readScalar(kin.lookup("Ea1_J_mol"))' not in solver: errors.append('missing LCO kinetics')
    result={'static_case_syntax_check_passed':not errors,
       'checked_divergence_fields':sorted(set(fields)),
       'errors':errors,'openfoam_build_verified':False,
       'limitations':['This is a text-level preflight, not an OpenFOAM parser.',
         'The legacy v11 solver is INCOMPRESSIBLE; not a verified compressible gas-flow solver.',
         'The v12 Python model is separate from this C++ solver.']}
    return result
if __name__=='__main__':
    r=lint(); print(json.dumps(r,indent=2)); raise SystemExit(0 if not r['errors'] else 2)
