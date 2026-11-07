"""Build actual figures and an evidence manifest from locally computed v8 results."""
from __future__ import annotations
import json, shutil, subprocess
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v8.thermal_em import solve as em_thermal
from v8.lco_holdout import main as holdout
from v8.physics import ideal_gas_mixture,molar_density_ratio

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'v8/results'

def generate():
    OUT.mkdir(parents=True,exist_ok=True)
    er,iterations=em_thermal(write=True,out=OUT/'em_thermal')
    ho=holdout(OUT/'lco_holdout.json')
    gas=ideal_gas_mixture(.6,.2,850,101325)
    packed=np.load(OUT/'em_thermal/fields.npz')
    for name,title,unit in [('Ts','EM–thermal converged solid temperature','K'),('qEM','Induction Joule heating source','W/m³')]:
        fig,ax=plt.subplots(figsize=(7.5,3.5))
        im=ax.imshow(packed[name],origin='lower',aspect='auto',extent=[packed['r'][0]*1000,packed['r'][-1]*1000,packed['z'][0]*1000,packed['z'][-1]*1000]);fig.colorbar(im,ax=ax,label=unit)
        ax.set_xlabel('Radius (mm)');ax.set_ylabel('Axial position (mm)');ax.set_title(title)
        fig.tight_layout();fig.savefig(OUT/f'{name}_v8.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.8,3.6))
    ax.semilogy([x['iteration'] for x in iterations],[x['max_temperature_change_K'] for x in iterations],marker='o')
    ax.set_xlabel('EM↔thermal fixed-point iteration');ax.set_ylabel('Maximum temperature change (K)');ax.grid(True,alpha=.25)
    fig.tight_layout();fig.savefig(OUT/'em_thermal_convergence.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.7,3.8))
    for name,marker in [('train','o'),('held_out','s')]:
        rows=ho[name]['rows'];ax.scatter([x['observed_alpha_240min_approx'] for x in rows],[x['predicted_alpha_240min'] for x in rows],marker=marker,label=name)
        for x in rows:ax.annotate(str(x['T_C'])+'°C',(x['observed_alpha_240min_approx'],x['predicted_alpha_240min']),xytext=(3,5),textcoords='offset points',fontsize=8)
    ax.plot([0,1],[0,1],linestyle='--',label='perfect agreement');ax.set_xlim(0,1.03);ax.set_ylim(0,1.06)
    ax.set_xlabel('Article-reported approximate conversion at 240 min');ax.set_ylabel('Illustrative model conversion');ax.legend()
    fig.tight_layout();fig.savefig(OUT/'lco_temperature_holdout.png',dpi=170);plt.close(fig)
    report={'version':'8.0','date':'2026-11-05','model_corrections':[
         'Ni product layer molar density replaces NiO density in volume reaction mapping',
         'H2 ideal gas mole fraction uses gas temperature Tg',
         '2D axisymmetric EM / thermal under-relaxed fixed point with FV conduction and grey-wall radiation',
         'Offline geometry-aware OpenFOAM qEM feedback map from ASCII temperature field',
         'Approximate LiCoO2 two-stage kinetics with temperature holdout',
         'PISO-like pressure candidate with projected phi and continuity log; NOT COMPILED'],
      'measurable_results':{'em_thermal':er,'lco_held_out':ho['held_out'],
         'ni_vs_nio_molar_density_ratio':molar_density_ratio(),'gas_mixture_850K':gas},
      'research_claims':{'locally_numerically_verified':True,'openfoam_compiled':bool(shutil.which('wmake')) and bool(shutil.which('blockMesh')),
           'openfoam_runtime_verified':False,'independent_reactor_experiment_validated':False,
           'fully_coupled_plasma_mhd':False,'moving_particle_cfd_dem':False},
      'warning':'This model contains illustrative/uncalibrated assumptions; test pass != process validation.'}
    (OUT/'EVIDENCE_MANIFEST.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
if __name__=='__main__':
    r=generate();print(json.dumps({'em_thermal':r['measurable_results']['em_thermal'],'lco_held_out':r['measurable_results']['lco_held_out'],'research_claims':r['research_claims']},indent=2))
