"""Compact repeatable v4 resolution/feedback studies; avoid false convergence claims."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v4.coupled_reactor import run

def main():
    base=Path('v4/results');base.mkdir(exist_ok=True,parents=True)
    experiments={
        'grid_20x8':dict(nx=20,ny=8,dt=.002,end=3.),
        'grid_24x10':dict(nx=24,ny=10,dt=.002,end=3.),
        'grid_32x12':dict(nx=32,ny=12,dt=.002,end=3.),
        'time_dt001':dict(nx=24,ny=10,dt=.001,end=3.),
        'fixed_drag':dict(nx=24,ny=10,dt=.002,end=3.,back_coupling=False),
        'no_reaction':dict(nx=24,ny=10,dt=.002,end=3.,reaction=False),
        'no_induction':dict(nx=24,ny=10,dt=.002,end=3.,em=False),
    }
    rows=[]
    for tag,params in experiments.items():
        summary=run(**params,out=str(base/'studies'/tag))
        rows.append(dict(case=tag,nx=summary['nx'],ny=summary['ny'],dt=summary['dt_s'],
                        mean_conversion=summary['mean_conversion'],max_Ts_K=summary['max_solid_temperature_K'],
                        mean_U=summary['mean_velocity_m_s'],div=summary['max_divergence_1_s'],
                        mass_rel=summary['h2_balance_relative'],energy_rel=summary['energy_balance_relative']))
    (base/'studies'/'overview.json').write_text(json.dumps(rows,indent=2)+'\n')
    fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    a=next(r for r in rows if r['case']=='grid_20x8')
    for ax,col,title in zip(axes.ravel(),['mean_conversion','max_Ts_K','mean_U','energy_rel'],
                           ['Mean NiO conversion (fraction)','Max solid temperature [K]','Mean axial U [m/s]','Energy balance relative residual']):
        names=[x['case'] for x in rows]
        values=[x[col] for x in rows]
        ax.bar(range(len(names)),values)
        ax.set_xticks(range(len(names)),names,rotation=65,ha='right',fontsize=7)
        ax.set_title(title,fontsize=10)
        ax.grid(axis='y',alpha=.2)
    fig.suptitle('V4 sensitivity screening: synthetic model, not experimental validation',fontsize=11)
    fig.savefig(base/'studies'/'sensitivity.png',dpi=170);plt.close(fig)
    a=np.load(base/'coupled'/'fields.npz')
    fig,axes=plt.subplots(1,3,figsize=(12,3.8),constrained_layout=True)
    for ax,field,title in zip(axes,['Ts','X','Ux'],['Solid T [K]','NiO conversion','Axial velocity [m/s]']):
        h=ax.imshow(a[field].T,origin='lower',aspect='auto',extent=[0,.12,0,.02])
        ax.set_title(title);ax.set_xlabel('x [m]');ax.set_ylabel('y [m]')
        fig.colorbar(h,ax=ax)
    fig.suptitle('Computed fields (v4 coupled finite-volume/projection model)',fontsize=11)
    fig.savefig(base/'coupled'/'reactor_fields.png',dpi=170);plt.close(fig)
    print('Wrote',str(base/'studies'/'overview.json'))

if __name__=='__main__':main()
