"""Reproducible v9 scientific figures and machine-readable claims manifest."""
from __future__ import annotations
from pathlib import Path
import sys,platform,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v9.plasma_energy import run as plasma_run
from v9.lco_process import scenario as lco_scenario
from v9.validation_protocol import leave_one_temperature_out
from v9.openfoam.preflight import report as openfoam_report


def main(out='v9/results'):
    base=Path(out);base.mkdir(exist_ok=True,parents=True)
    plasma,hist=plasma_run()
    loocv=leave_one_temperature_out()
    pre=openfoam_report(base/'openfoam_environment.json')
    h=np.asarray([r['H_mol_m3'] for r in hist]); h2=np.asarray([r['H2_mol_m3'] for r in hist]); tt=np.asarray([r['time_s'] for r in hist]);
    plt.figure(figsize=(7,4))
    plt.plot(tt,2*h2,label='H atoms in H₂ (mol atoms m⁻³)')
    plt.plot(tt,h,label='H radical atoms (mol m⁻³)')
    plt.plot(tt,2*h2+h,':',label='Total H atoms')
    plt.xlabel('Time (s)');plt.ylabel('H-atom inventory (mol m⁻³)')
    plt.title('Numerical H₂/H atom conservation; imposed electron power')
    plt.legend(fontsize=8);plt.grid(alpha=.2);plt.tight_layout();plt.savefig(base/'plasma_conservation.png',dpi=180);plt.close()
    tt=np.linspace(0,240,40)
    plt.figure(figsize=(7,4))
    for temp in (700,800,1000,1100):
        alpha=[lco_scenario(temp,float(x),10,7.5)['carbon_limited_conversion_LCO'] for x in tt]
        plt.plot(tt,alpha,label=f'{temp} °C')
    plt.xlabel('Time (min)');plt.ylabel('Illustrative LiCoO₂ conversion')
    plt.title('Assumed two-stage kinetics; NOT reactor-validated')
    plt.legend();plt.grid(alpha=.2);plt.tight_layout();plt.savefig(base/'lco_scenario_curves.png',dpi=180);plt.close()
    plt.figure(figsize=(5.7,5))
    x=np.asarray([r['observed_conversion_APPROX'] for r in loocv['cases']]); y=np.asarray([r['predicted_conversion'] for r in loocv['cases']]);
    plt.plot([0,1],[0,1],'--',label='Parity')
    plt.scatter(x,y,s=55)
    for o,p,r in zip(x,y,loocv['cases']):
        plt.annotate(f"{r['held_out_temperature_C']}°C",(o,p),xytext=(4,5),textcoords='offset points',fontsize=9)
    plt.xlim(0,1.05);plt.ylim(0,1.05)
    plt.xlabel('Approximate source endpoint conversion');plt.ylabel('One-temperature-held-out prediction')
    plt.title('Same-study endpoint holdout; not external validation')
    plt.tight_layout();plt.savefig(base/'lco_holdout_parity.png',dpi=180);plt.close()
    flags={'v9_python_modules_tested':True,
           'OpenFOAM_runtime_executed_locally':False,
           'OpenFOAM_prerequisites_present':pre['prerequisites_present'],
           'full_plasma_MHD_resolved':False,
           'H_radical_energy_atom_ledger_tested':True,
           'LCO_kinetic_prefactors_independently_validated':False,
           'independent_reactor_data_available':False}
    manifest={'title':'ReInduct CFD Pro v9 local evidence','python_version':sys.version.split()[0],
              'platform':platform.platform(), 'verification_flags':flags,
              'plasma_H_atom_abs_error_max_mol_m3':plasma['atomic_H_balance_abs_max_mol_m3'],
              'plasma_energy_abs_error_max_J_m3':plasma['total_energy_balance_abs_max_J_m3'],
              'LCO_same_study_LOOCV_MAE':loocv['summary_MAE_conversion'],
              'LCO_endpoint_source_quality':'approximate; not a digitized independent experiment',
              'numerical_figures':['plasma_conservation.png','lco_scenario_curves.png','lco_holdout_parity.png'],
              'do_not_claim':['successful OpenFOAM runtime','electron transport or complete plasma','independent reactor validation','predictive real black-mass recovery']}
    (base/'EVIDENCE_MANIFEST_V9.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2));return manifest

if __name__=='__main__':main()
