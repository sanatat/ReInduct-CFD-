"""Reproducible scientific evidence for the v6 axisymmetric EM+reactor prototypes."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v6.axisymmetric_em import Grid,manufactured,run as em_run
from v6.axisymmetric_bed import run as bed_run
from v6.kinetics_cv import crossvalidate
from v6.battery_flowsheet import study as battery_study

def main(out='v6/results'):
    d=Path(out);d.mkdir(exist_ok=True,parents=True)
    mms=[];power=[]
    for n in [12,18,26,40]:
        g=Grid(nr=n,nz=round(n*1.6))
        m=manufactured(g)
        mms.append({'radial_cells':n,'axial_cells':g.nz,'relative_L2':m['relative_L2'],
                    'algebraic_residual':m['algebraic_relative_residual']})
    for n in [16,24,36,48]:
        g=Grid(nr=n,nz=round(n*1.6))
        s,f=em_run(g)
        power.append({'radial_cells':n,'axial_cells':g.nz,
                      'power_W':s['joule_power_W'],'max_B_T':s['peak_B_T'],
                      'residual':s['algebraic_relative_residual']})
    s,f=em_run(Grid(nr=48,nz=77),write=True,out=d/'axisymmetric')
    with (d/'axisymmetric'/'verification.json').open('w') as fp:
        json.dump({'manufactured':mms,'power_refinement':power,
                   'power_relative_change_last_two':abs(power[-1]['power_W']/power[-2]['power_W']-1),
                   'claim':'discretization verification; no experimental EM validation'},fp,indent=2)
    fig,axs=plt.subplots(1,3,figsize=(14,4))
    img=axs[0].pcolormesh(f['z']*1000,f['r']*1000,f['q'].T/1e6,shading='auto')
    fig.colorbar(img,ax=axs[0],label='MW/m³');axs[0].set_title('Axisymmetric induced Joule heat')
    axs[0].set_ylabel('Radius [mm]');axs[0].set_xlabel('Axial z [mm]')
    axs[1].loglog([x['radial_cells'] for x in mms],[x['relative_L2'] for x in mms],'o-')
    axs[1].set_title('Manufactured PDE verification');axs[1].set_xlabel('N radial');axs[1].set_ylabel('relative L2')
    axs[2].plot([x['radial_cells'] for x in power],[x['power_W'] for x in power],'o-')
    axs[2].set_title('Integrated Joule power refinement');axs[2].set_xlabel('N radial');axs[2].set_ylabel('W')
    for a in axs[1:]:a.grid(alpha=.2)
    fig.tight_layout();fig.savefig(d/'axisymmetric'/'em_verification.png',dpi=150);plt.close(fig)
    comparisons={}
    for tag, opts in [('induction_reaction',dict(induction=True,reaction=True)),
                      ('reaction_no_induction',dict(induction=False,reaction=True)),
                      ('nonreacting_no_induction',dict(induction=False,reaction=False))]:
        comparisons[tag]=bed_run(nz=24,nr=10,end_s=2.0,dt_s=.004,
                                 em_nr=36,em_nz=58,out=d/tag,**opts)
    for n in [12,18,24]:
        sr=bed_run(nz=n,nr=max(5,n//2),end_s=.8,dt_s=.004,
                   em_nr=20,em_nz=32,induction=True,out=d/f'grid_{n}',write=True)
        comparisons[f'grid_{n}']=sr
    (d/'reactor_cases.json').write_text(json.dumps(comparisons,indent=2)+'\n')
    base=np.load(d/'induction_reaction'/'fields.npz')
    fig,axs=plt.subplots(1,3,figsize=(13.5,3.9))
    for ax,k,lab in zip(axs,['Ts','X','qEM'],['Solid temperature [K]','NiO model conversion [-]','Joule heat [MW/m³]']):
        v=base[k] if k!='qEM' else base[k]/1e6
        mesh=ax.pcolormesh(base['z']*1000,base['r']*1000,v.T,shading='auto')
        fig.colorbar(mesh,ax=ax);ax.set_title(lab);ax.set_xlabel('z [mm]');ax.set_ylabel('r [mm]')
    fig.tight_layout();fig.savefig(d/'induction_reaction'/'reactor_fields.png',dpi=150);plt.close(fig)
    cv=crossvalidate(out=d/'kinetics')
    battery_study(out=d/'battery')
    report={'em':s,'reactor':comparisons,'kinetics_CV_MAE_um_s':{
        k:v['pooled_MAE_um_s'] for k,v in cv.items() if isinstance(v,dict)},
        'em_verification_mms_last':mms[-1],'power_refinement_last_two':power[-2:],
        'disclaimers':['Source parameters approximate/illustrative; no experimental reactor validation',
                       'NiO/H2 is not battery black mass thermochemistry',
                       'Axisymmetric bed flow: Ergun closure, not full CFD Navier-Stokes',
                       'OpenFOAM code in project still pending first external compile/run',
                       'Kinetics CV uses a single study and a different observable (Ni growth rate).']}
    (d/'v6_evidence.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'MMS_error_last':mms[-1]['relative_L2'],
     'power_relative_change_last_two':abs(power[-1]['power_W']/power[-2]['power_W']-1),
     'reactor_Tmax_K':comparisons['induction_reaction']['max_solid_temperature_K'],
     'reactor_conversion':comparisons['induction_reaction']['mean_NiO_conversion'],
     'kinetics_CV_MAE_um_s':report['kinetics_CV_MAE_um_s']},indent=2))
    return report
if __name__=='__main__':main()
