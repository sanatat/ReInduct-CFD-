"""Plots computed from v10 solver outputs; no illustrative or fabricated data."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v10.plasma_transport_1d import Settings,initialize
from v10.lco_published_rates import EARLY,LATE,arrhenius_fit

def create(out='v10/results'):
    p=Path(out);p.mkdir(parents=True,exist_ok=True)
    c=Settings();ini=initialize(c)
    a=np.load(p/'plasma1d.npz')
    fig,ax=plt.subplots(figsize=(9,4.8))
    ax.plot(ini['x'],ini['h'],label='Initial H radicals')
    ax.plot(a['x'],a['h'],label='Final H radicals')
    ax.set(xlabel='Periodic axial coordinate (dimensionless)',ylabel='H radical concentration (normalized)',title='1D reaction–drift–diffusion–Poisson: H radical transport')
    ax.legend();ax.grid(alpha=.25);fig.tight_layout();fig.savefig(p/'plasma_radical_profile.png',dpi=160);plt.close(fig)

    fig,ax=plt.subplots(figsize=(9,4.8))
    ax.plot(a['x'],a['E'],label='Electric field')
    ax.plot(a['x'],a['potential'],label='Electrostatic potential')
    ax.set(xlabel='Periodic axial coordinate (dimensionless)',ylabel='Normalized field / potential',title='Self-consistent Poisson solution (non-dimensional verification)')
    ax.legend();ax.grid(alpha=.25);fig.tight_layout();fig.savefig(p/'plasma_field.png',dpi=160);plt.close(fig)

    d=np.load(p/'hartmann.npz')
    fig,ax=plt.subplots(figsize=(9,4.8))
    ax.plot(d['y']*1e3,d['exact'],label='Analytical Hartmann velocity')
    ax.plot(d['y']*1e3,d['u'],'--',label='Finite-difference velocity')
    ax.set(xlabel='Channel transverse position (mm)',ylabel='Axial velocity (m/s)',title='Hartmann flow: analytical vs numerical, B = 3 T')
    ax.legend();ax.grid(alpha=.25);fig.tight_layout();fig.savefig(p/'hartmann_comparison.png',dpi=160);plt.close(fig)

    fig,ax=plt.subplots(figsize=(9,4.8))
    for data,label in ((EARLY,'Early stage GB (published)'),(LATE,'Late stage A3 (published)')):
        T=np.array([r[0]+273.15 for r in data]);k=np.array([r[2] for r in data]);fit=arrhenius_fit(data)
        tx=np.linspace(min(T),max(T),100)
        slope=-fit['E_a_J_per_mol']/8.31446261815324
        intercept=np.mean(np.log(k)-slope/T)
        ax.scatter(1000/T,np.log(k),label=label)
        ax.plot(1000/tx,intercept+slope/tx,label=f'Arrhenius {label.split()[0]} fit')
    ax.set(xlabel='1000/T (K⁻¹)',ylabel='ln(k / min⁻¹)',title='Published LiCoO2+C stage rate constants (same-study regression)')
    ax.legend();ax.grid(alpha=.25);fig.tight_layout();fig.savefig(p/'lco_published_rate_arrhenius.png',dpi=160);plt.close(fig)
    return [str(x) for x in p.glob('*.png') if x.name in {'plasma_radical_profile.png','plasma_field.png','hartmann_comparison.png','lco_published_rate_arrhenius.png'}]

if __name__=='__main__':
 print('\n'.join(create()))
