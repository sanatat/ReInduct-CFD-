"""Create bounded, evidence-backed release manifest. No false successful statuses."""
from pathlib import Path
import json,re,hashlib,datetime

def main():
    root=Path('v10/results')
    runlog=(root/'all_tests.log').read_text()
    m=re.search(r'(\d+) passed',runlog)
    tests=int(m.group(1)) if m else None
    foamexit=int((root/'openfoam_exit_code.txt').read_text().strip())
    foamverified=foamexit==0 and (root/'openfoam_runtime/VERIFIED.txt').exists()
    from v10.lco_published_rates import SOURCE
    reports={'plasma':json.loads((root/'plasma1d.json').read_text())['report'],
             'hartmann':json.loads((root/'hartmann.json').read_text()),
             'rate':json.loads((root/'rate_validation.json').read_text()),
             'external':json.loads((root/'external_phase_audit.json').read_text())}
    data={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'reexecuted_python_tests_passed':tests,
          'standalone_cpp_compiled_and_ran':(root/'cpp_kinetics_log.txt').exists(),
          'openfoam_v2606_runtime_verified':foamverified,
          'openfoam_exit_code':foamexit,
          'external_reactor_conversion_validated':False,
          'physically_calibrated_plasma':False,
          'published_rate_source':SOURCE,
          'rate_models_are_from_one_published_experiment':True,
          'independent_phase_study_is_qualitative_only':True,
          'primary_evidence':{
            'plasma_H_atom_inventory_error':reports['plasma']['absolute_H_atom_inventory_error'],
            'plasma_integrated_charge_error':reports['plasma']['absolute_integrated_charge_error'],
            'plasma_total_energy_including_field_error':reports['plasma']['absolute_total_energy_including_Poisson_field_error'],
            'Hartmann_relative_L2_flow_error':reports['hartmann']['relative_L2_flow_error'],
            'Hartmann_mechanical_budget_error':reports['hartmann']['mechanical_budget_relative_error'],
            'early_stage_kinetic_activation_J_mol':reports['rate']['early_fit']['E_a_J_per_mol'],
            'late_stage_kinetic_activation_J_mol':reports['rate']['late_fit']['E_a_J_per_mol'],
            'independent_study_cobalt_phase_mismatch':reports['external']['metallic_cobalt_vs_cobalt_oxide_mismatch'],
          },
          'file_sha256':{}}
    for file in sorted(root.rglob('*')):
        if file.is_file() and 'runtime/work' not in str(file) and file.name not in ('evidence_manifest.json',):
            data['file_sha256'][str(file.relative_to(root))]=hashlib.sha256(file.read_bytes()).hexdigest()
    (root/'evidence_manifest.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:v for k,v in data.items() if k!='file_sha256'},indent=2))
    return data

if __name__=='__main__':main()
