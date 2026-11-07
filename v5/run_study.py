"""Generate all v5 evidence from executed local code. No OpenFOAM claims."""
from __future__ import annotations
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from v5.em2d import convergence
from v5.kinetics_validation import assess
from v5.thermochemistry import report as thermo_report
from v4.battery_chemistry import batch
from v4.coupled_reactor import run as reactor

def run_all():
    out=Path('v5/results');out.mkdir(parents=True,exist_ok=True)
    em,mms,mesh=convergence(out/'em2d')
    kinetics=assess(out/'kinetics')
    thermo=thermo_report(out/'thermo')
    noem=reactor(nx=24,ny=10,dt=.002,end=1.2,em=False,em_backend='2d',thermo_backend='tabulated',out=out/'reactor_no_em')
    with_em=reactor(nx=24,ny=10,dt=.002,end=1.2,em=True,em_backend='2d',thermo_backend='tabulated',out=out/'reactor_em2d')
    with_mhd=reactor(nx=24,ny=10,dt=.002,end=1.2,em=True,em_backend='2d',thermo_backend='tabulated',mhd_sigma_g=5.,out=out/'reactor_mhd_proxy')
    q=np.load(out/'reactor_em2d'/'fields.npz')
    fig,axs=plt.subplots(1,3,figsize=(12,3.5))
    for a,title,arr in zip(axs,['Solid T [K]','Conversion X','Inductive heat [W/m³]'],
                            [q['Ts'],q['X'],np.broadcast_to(q['qem'],q['Ts'].shape)]):
        im=a.imshow(arr.T,origin='lower',extent=(0,.12,0,.02),aspect='auto',cmap='viridis')
        a.set_title(title);a.set_xlabel('x [m]');fig.colorbar(im,ax=a)
    axs[0].set_ylabel('y [m]')
    fig.tight_layout();fig.savefig(out/'reactor_em2d'/'coupled_fields.png',dpi=155);plt.close(fig)
    combined={
      'EM_2D':em,'manufactured_solution':mms,'EM_coil_mesh_convergence':mesh,
      'kinetics_26pct_H2_holdout_MAE_um_s':{name:value['heldout_MAE_um_s'] for name,value in kinetics.items() if isinstance(value,dict) and 'heldout_MAE_um_s' in value},
      'thermochemistry_Kp_at_800_K':thermo['values'][1]['Kp'],
      'reactor_with_em':with_em,'reactor_without_em':noem,'reactor_with_illustrative_mhd_proxy':with_mhd,
      'mhd_mean_velocity_relative_change':(with_mhd['mean_velocity_m_s']/with_em['mean_velocity_m_s']-1.),
      'mean_conversion_absolute_difference_with_em':with_em['mean_conversion']-noem['mean_conversion'],
      'max_temperature_K_difference_with_em':with_em['max_solid_temperature_K']-noem['max_solid_temperature_K'],
      'lco_atomic_balance_checked':batch()['atomic_balance'],
      'openfoam_2606_compiled':False,'openfoam_2606_run':False,
      'reactor_experimental_validation':False,
      'full_plasma_mhd_coupling':False}
    (out/'v5_summary.json').write_text(json.dumps(combined,indent=2)+'\n')
    print('FINAL_V5_SUMMARY',json.dumps({
       'qavg_W_m3':with_em['last_joule_heat_source_mean_W_m3'],
       'max_solid_T_K':with_em['max_solid_temperature_K'],
       'conversion':with_em['mean_conversion'],
       'delta_T_with_em':combined['max_temperature_K_difference_with_em'],
       'h_balance':with_em['h2_balance_relative'],
       'energy_balance':with_em['energy_balance_relative']},indent=2))
    return combined
if __name__=='__main__':run_all()
