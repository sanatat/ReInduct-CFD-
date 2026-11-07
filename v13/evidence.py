"""Fail-closed provenance: numerical verification != physical validation."""
import hashlib,json,platform,shutil,sys
from pathlib import Path
import numpy,scipy

ROOT=Path(__file__).resolve().parents[1]

def write(out='v13/results/EVIDENCE_MANIFEST.json'):
    sources=[ROOT/'v13'/s for s in ('reactor.py','compressible_transport.py','thermo_pressure.py','gas_properties.py','em_feedback.py','study.py','cpp/gas_shomate.cpp')]
    values={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    case=ROOT/'v13/results/final_5s/evidence.json'
    evidence=json.loads(case.read_text()) if case.exists() else None
    toolchain={k:bool(shutil.which(k)) for k in ('wmake','blockMesh','checkMesh','docker','g++')}
    from v12.openfoam.case_lint import lint as legacy_lint
    static_check=legacy_lint()
    manifest={
      'version':'13.0', 'scope':'Numerically verified reduced-order research surrogate',
      'python':sys.version.split()[0],'platform':platform.platform(),
      'numpy':numpy.__version__,'scipy':scipy.__version__,
      'source_sha256':values,
      'numerical_evidence':{k:evidence[k] for k in (
       'mean_LCO_conversion','max_element_balance_mol_m2',
       'relative_energy_ledger_error','max_end_step_EOS_residual_Pa',
       'max_pressure_energy_Picard_iterations',
       'max_EM_feedback_relative_residual','max_EM_feedback_iterations_per_refresh')}
       if evidence else 'run model first',
      'toolchain':toolchain,
      'legacy_openfoam_static_case_lint_passed':bool(static_check['static_case_syntax_check_passed']),
      'openfoam_v13_runtime_verified':False,
      'openfoam_v11_legacy_runtime_verified':False,
      'full_plasma_mhd_verified':False,
      'experimental_reactor_validated':False,
      'experimental_black_mass_recovery_validated':False,
      'reason':'No validated OpenFOAM v13 solver, full plasma-MHD or independent matched reactor data.'}
    dest=ROOT/out;dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest
if __name__=='__main__':
    m=write();print(json.dumps({'toolchain':m['toolchain'],'sha256_count':len(m['source_sha256']),
      'openfoam_v13_runtime_verified':m['openfoam_v13_runtime_verified']},indent=2))
