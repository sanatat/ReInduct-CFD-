#!/usr/bin/env python3
"""Reproduce v3 numerical verification and sensitivity tables/figures."""
from pathlib import Path
import sys,csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from cfd.mac_channel import simulate as channel
from reactor2d.reactor import simulate as reactor


def make_plot(x,y,xlab,ylab,title,path):
    fig,ax=plt.subplots(figsize=(7.2,4.4))
    ax.plot(x,y,'o-',lw=2)
    ax.set(xlabel=xlab,ylabel=ylab,title=title)
    ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(path,dpi=160);plt.close(fig)

def heatmap(a,x,y,xlab,ylab,title,path):
    fig,ax=plt.subplots(figsize=(7.4,4.2))
    im=ax.pcolormesh(x,y,a.T,shading='auto')
    fig.colorbar(im,ax=ax)
    ax.set(xlabel=xlab,ylabel=ylab,title=title)
    fig.tight_layout();fig.savefig(path,dpi=160);plt.close(fig)

def main():
    out=ROOT/'results/v3/evidence';out.mkdir(parents=True,exist_ok=True)
    print('\n=== 2-D MAC NS + Brinkman gas-flow benchmark ===',flush=True)
    flow=channel(nx=48,ny=16,dt=.003,end=7,porous_drag=.7,out=out/'ns_brinkman')
    print('\n=== 2-D MAC Hartmann MHD benchmark ===',flush=True)
    mhd=channel(nx=48,ny=16,dt=.003,end=5,B_T=.1,out=out/'hartmann')
    print('\n=== Coupled MAC velocity + EM/T + redox ===',flush=True)
    coupled=reactor(nx=48,ny=16,dt=.003,end=12,velocity_file=out/'ns_brinkman/velocity.npz',out=out/'reacting')
    print('\n=== One-way vs two-way sigma(T) sensitivity ===',flush=True)
    fixed=reactor(nx=48,ny=16,dt=.003,end=12,velocity_file=out/'ns_brinkman/velocity.npz',feedback=False,out=out/'fixed_sigma')
    print('\n=== 2-D reaction/heat mesh-refinement study ===',flush=True)
    grid=[]
    for n,m in [(24,8),(48,16),(72,24)]:
        s=reactor(nx=n,ny=m,dt=.0025,end=6,velocity_model='plug',out=out/f'grid_{n}x{m}')
        grid.append((n,m,s['mean_conversion'],s['max_solid_T_K'],s['hydrogen_balance_rel'],s['energy_balance_rel']))
    with (out/'mesh_refinement.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['nx','ny','mean_conversion','max_solid_T_K','h_balance_relative','energy_balance_relative']);w.writerows(grid)
    make_plot([n*m for n,m,*_ in grid],[i[2] for i in grid],
              'Number of cells','Mean NiO conversion [-]','2-D FVM mesh refinement',out/'mesh_refinement.png')
    print('\n=== Time-step sensitivity (constant geometry) ===',flush=True)
    time=[]
    for dt in [.004,.002,.001]:
        s=reactor(nx=36,ny=12,dt=dt,end=4,velocity_model='plug',out=out/f'dt_{dt}')
        time.append((dt,s['mean_conversion'],s['max_solid_T_K']))
    with (out/'timestep_refinement.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['dt_s','mean_conversion','max_solid_T_K']);w.writerows(time)
    print('\n=== Induction field / conversion sensitivity ===',flush=True)
    fields=[]
    for h0 in [0,3000,6000,10000,14000]:
        s=reactor(nx=36,ny=12,dt=.003,end=5,h0=h0,no_em=h0==0,out=out/f'H0_{h0}')
        fields.append((h0,s['mean_conversion'],s['max_solid_T_K'],s['cumulative_EM_heat_J_per_m']))
    with (out/'H0_sensitivity.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['H0_A_m','mean_conversion','max_solid_T_K','cumulative_EM_heat_J_per_m']);w.writerows(fields)
    make_plot([a for a,*_ in fields],[a[1] for a in fields],
              'Tangential H-field boundary amplitude [A/m]', 'Maximum solid temperature [K]',
              'Induction amplitude sensitivity (illustrative)',out/'H0_temperature.png')
    f=np.load(out/'reacting/fields.npz')
    heatmap(f['NiO_conversion'],f['x_m'],f['y_m'],'Axial x [m]','Cross-channel y [m]',
            'Solid conversion X after 12 s',out/'conversion_2d.png')
    heatmap(f['Tsolid_K'],f['x_m'],f['y_m'],'Axial x [m]','Cross-channel y [m]',
            'Solid temperature [K] after 12 s',out/'solid_temperature_2d.png')
    flowdata=np.load(out/'ns_brinkman/velocity.npz')
    heatmap(flowdata['u_xfaces_m_s'][:-1,:],f['x_m'],f['y_m'],
            'Periodic channel x [m]','Cross-channel y [m]',
            '2-D Navier-Stokes x-velocity [m/s]',out/'flow_velocity_2d.png')
    with (out/'reacting/history.csv').open() as fp:
        hist=np.genfromtxt(fp,delimiter=',',names=True)
    make_plot(hist['time_s'],hist['conversion_mean'],'Simulated time [s]',
              'Mean NiO conversion [-]','Reaction conversion history',out/'conversion_history.png')
    make_plot(f['y_m'],f['qem_W_m3'],'Cross-channel y [m]',
              'Volumetric induced Joule heat [W/m³]','Induction energy source profile',out/'induction_source.png')
    report={'flow':flow,'mhd':mhd,'reactor':coupled,'fixed_sigma':fixed,
            'grid_refinement':grid,'time_refinement':time,'h0_sensitivity':fields,
            'scope':'all illustrative model-based numerical tests; not experimental validation'}
    (out/'study.json').write_text(json.dumps(report,indent=2)+'\n')
    print('\nSTUDY COMPLETE:',out)
if __name__=='__main__':main()
