"""Run verification studies, held-out published rate comparisons and plots."""
from pathlib import Path
from dataclasses import replace
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v11.integrated_reactor import Config, simulate
from v10.lco_published_rates import build as published_rate_check

def execute(out='v11/results'):
    dest=Path(out);dest.mkdir(parents=True,exist_ok=True)
    c=Config()
    baseline=simulate(c,dest/'base')
    noem=simulate(replace(c,enable_em=False),dest/'no_em')
    nochem=simulate(replace(c,enable_chemistry=False),dest/'no_chem')
    # Explicit Euler stability: dt decreases with grid size, with fixed physical EM grid.
    mesh=[]
    for nz,dt in [(16,.03),(24,.03),(40,.03)]:
        run=simulate(replace(c,nz=nz,dt_s=dt,duration_s=60.,em_grid_nr=16,em_grid_nz=24),dest/f'grid_{nz}')
        mesh.append({'nz':nz,'dt_s':dt,'conversion':run['mean_conversion_LCO'],
                     'max_Ts_K':run['max_final_Ts_K'],
                     'max_element_error':run['max_abs_element_residual'],
                     'relative_energy_error':run['relative_energy_residual']})
    # Fixed-grid time-step refinement for time accuracy separate from mesh refinement.
    timesteps=[]
    for dt in (.16,.08,.04):
        run=simulate(replace(c,nz=20,dt_s=dt,duration_s=48.,em_grid_nr=16,em_grid_nz=24),dest/f'dt_{str(dt).replace(".", "p")}')
        timesteps.append({'dt_s':dt,'conversion':run['mean_conversion_LCO'],'max_Ts_K':run['max_final_Ts_K']})
    rates=published_rate_check(dest/'published_rates')
    checks={'baseline':{'conversion':baseline['mean_conversion_LCO'],'max_Ts_K':baseline['max_final_Ts_K'],
       'mean_Ts_K':baseline['mean_final_Ts_K'],'element_error':baseline['max_abs_element_residual'],
       'relative_energy_error':baseline['relative_energy_residual'],
       'max_ideal_gas_pressure_deviation':baseline['max_inferred_ideal_gas_pressure_deviation_fraction'],
       'em_initial_W':baseline['em_updates'][0]['annulus_EM_power_W'],
       'em_last_W':baseline['em_updates'][-1]['annulus_EM_power_W']},
       'no_EM':{'conversion':noem['mean_conversion_LCO'],'max_Ts_K':noem['max_final_Ts_K']},
       'no_chem':{'conversion':nochem['mean_conversion_LCO']},
       'grid_study':mesh,'time_step_study':timesteps,
       'literature_same_study_early_holdout_logMAE':rates['early_holdouts']['logMAE'],
       'literature_same_study_late_holdout_logMAE':rates['late_holdouts']['logMAE'],
       'independent_experimental_reactor_validation':False,
       'openfoam_run':False,
       'kinetic_transfer_warning':'Published stage model k in min^-1 are not local gas-solid mechanistic constants.'}
    (dest/'v11_summary.json').write_text(json.dumps(checks,indent=2)+'\n')
    fig,axs=plt.subplots(2,1,figsize=(8.2,7.8),sharex=True)
    for name,label in [('base','EM + chemistry'),('no_em','Without EM')]:
        r=json.loads((dest/name/'evidence.json').read_text())
        s=r['snapshots']
        axs[0].plot([x['time_s'] for x in s],[x['mean_Tsolid_K'] for x in s],marker='.',label=label)
        axs[1].plot([x['time_s'] for x in s],[100*x['mean_LCO_conversion'] for x in s],marker='.',label=label)
    axs[0].set_ylabel('Mean solid T (K)');axs[0].legend();axs[0].grid(alpha=.25)
    axs[1].set_ylabel('LiCoO2 conversion (%)');axs[1].set_xlabel('Time (s)');axs[1].grid(alpha=.25)
    fig.suptitle('V11 conservative illustrative model: coupled vs no induction')
    fig.tight_layout();fig.savefig(dest/'integrated_comparison.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.2,4))
    ax.plot([r['nz'] for r in mesh],[r['conversion']*100 for r in mesh],'-o')
    ax.set_xlabel('Axial FV cells');ax.set_ylabel('LiCoO2 conversion (%)')
    ax.set_title('Spatial sensitivity (fixed dt=0.03s; EM grid fixed)')
    ax.grid(alpha=.25);fig.tight_layout();fig.savefig(dest/'grid_study.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.2,4))
    ax.plot([r['dt_s'] for r in timesteps],[r['conversion']*100 for r in timesteps],'-o')
    ax.set_xlabel('Time step (s)');ax.set_ylabel('LiCoO2 conversion (%)')
    ax.set_title('Time-step sensitivity at fixed grid (20 cells)')
    ax.grid(alpha=.25);fig.tight_layout();fig.savefig(dest/'time_step_study.png',dpi=180);plt.close(fig)
    return checks

if __name__=='__main__':
    s=execute();print(json.dumps(s['baseline'],indent=2))
