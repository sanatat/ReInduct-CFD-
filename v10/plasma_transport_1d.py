"""PERIODIC 1-D *dimensionless* drift-diffusion-Poisson plasma verification.

H2, H, H2+, electrons, electron/neutral energy are conservative FV variables.
The Poisson operator uses a documented *dimensionless* coupling coefficient,
NOT physical vacuum permittivity or a validated reactor discharge model.
Averaged external electric field injects work J.E into electron energy.
Hydrogen, charge and total (thermal + chemical) energy budgets are tested.
The periodic boundary is chosen for code verification, not reactor realism.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import json
import numpy as np

DISSOC = 3.0 # energy units per H2 molecule equivalent
ION = 6.0    # energy units per H2+ ion equivalent

@dataclass(frozen=True)
class Settings:
    n: int = 80
    length: float = 1.0
    dt: float = 0.0004
    duration: float = 0.8
    e_external: float = 0.4
    poisson_coupling: float = 2.0
    mu_e: float = 0.2
    mu_i: float = 0.045
    d_e: float = 0.0002
    d_i: float = 0.0001
    d_neutral: float = 0.0001
    flow_speed: float = 0.04
    dissociation_rate: float = 0.1
    ionization_rate: float = 0.06
    dissoc_recomb_rate: float = 0.5
    radical_recomb_rate: float = 0.01
    electron_to_gas_exchange: float = 0.02
    external_heating: float = 0.03


def electric_field(ni: np.ndarray, ne: np.ndarray, cfg: Settings, return_faces=False):
    """Discrete-periodic Poisson, compatible with FV face-flux divergence.

    A non-zero net periodic charge is impossible under periodic boundary
    conditions; the initial condition and reactions keep mean charge at zero.
    """
    dx=cfg.length/cfg.n
    rho=cfg.poisson_coupling*(ni-ne)
    if abs(float(rho.mean()))>1e-10:
        raise ValueError('Periodic Poisson equation requires zero net charge')
    freq=np.fft.fftfreq(cfg.n, d=dx)
    eig=(2*np.sin(np.pi*freq*dx)/dx)**2
    phat=np.zeros(cfg.n,dtype=complex)
    hat=np.fft.fft(rho-rho.mean())
    phat[1:]=hat[1:]/eig[1:]
    phi=np.fft.ifft(phat).real
    eface=cfg.e_external-(np.roll(phi,-1)-phi)/dx
    ecell=.5*(eface+np.roll(eface,1))
    if return_faces:return ecell,phi,eface
    return ecell,phi


def transport(c, v, diff, dx, return_flux=False):
    """Periodic upwind + central diffusion FV: conserved integral exactly."""
    vf=.5*(v+np.roll(v,-1))
    right=np.roll(c,-1)
    flux=vf*np.where(vf>=0,c,right)-diff*(right-c)/dx
    rhs=-(flux-np.roll(flux,1))/dx
    return (rhs,flux) if return_flux else rhs


def initialize(cfg=Settings()):
    x=(np.arange(cfg.n)+.5)*cfg.length/cfg.n
    # ne and ni begin charge-neutral pointwise, with a small density variation.
    ne=0.009*(1+0.14*np.sin(2*np.pi*x/cfg.length))
    ni=ne.copy()
    return dict(x=x,h2=1+0.04*np.cos(2*np.pi*x/cfg.length),
                h=0.025*(1+0.09*np.sin(4*np.pi*x/cfg.length)),
                ni=ni,ne=ne,ee=4.0+0.3*np.cos(2*np.pi*x/cfg.length),
                eg=1.2+0.1*np.sin(2*np.pi*x/cfg.length))


def integrals(s, cfg):
    dx=cfg.length/cfg.n
    atoms=np.sum(2*s['h2']+s['h']+2*s['ni'])*dx
    charge=np.sum(s['ni']-s['ne'])*dx
    _,_,eface=electric_field(s['ni'],s['ne'],cfg,return_faces=True)
    electrostatic=np.sum((eface-cfg.e_external)**2)*dx/(2*cfg.poisson_coupling)
    energy=np.sum(s['ee']+s['eg']+.5*DISSOC*s['h']+ION*s['ni'])*dx+electrostatic
    return np.array([atoms,charge,energy])


def run(cfg=Settings(), save_path=None):
    if cfg.n<12 or cfg.dt<=0 or cfg.duration<=0 or cfg.length<=0:raise ValueError('Invalid grid/time')
    nstep=round(cfg.duration/cfg.dt)
    if abs(nstep*cfg.dt-cfg.duration)>1e-12:raise ValueError('duration must be multiple of dt')
    s=initialize(cfg)
    dx=cfg.length/cfg.n
    initial=integrals(s,cfg)
    added_energy=0.
    max_atom=0.;max_charge=0.;max_energy=0.;largest_field=0.;max_cfl=0.
    snapshots=[]
    for n in range(nstep):
        e,phi,eface=electric_field(s['ni'],s['ne'],cfg,return_faces=True)
        vi=cfg.flow_speed+cfg.mu_i*e
        ve=cfg.flow_speed-cfg.mu_e*e
        cfl=cfg.dt*(max(np.max(np.abs(ve)),np.max(np.abs(vi)),cfg.flow_speed)/dx +
                    2*max(cfg.d_e,cfg.d_i,cfg.d_neutral)/dx**2)
        max_cfl=max(max_cfl,float(cfl))
        if cfl>=.8:raise ValueError('Explicit FV stability margin exceeded')
        largest_field=max(largest_field,float(np.max(np.abs(e))))
        h2,h,ni,ne,ee,eg=(s[k] for k in ('h2','h','ni','ne','ee','eg'))
        # Derive current at precisely the same FV faces used in charge transport.
        ni_rhs,fi=transport(ni,vi,cfg.d_i,dx,return_flux=True)
        ne_rhs,fe=transport(ne,ve,cfg.d_e,dx,return_flux=True)
        face_charge_current=fi-fe
        face_work=face_charge_current*eface
        # J.E work includes change in stored electrostatic field energy. Only
        # uniform applied-field work is external; electrostatic field energy
        # is included in the total-energy diagnostic.
        joule_per_cell=.5*(face_work+np.roll(face_work,1))
        applied_work=cfg.e_external*np.mean(face_charge_current)
        # Reaction extents per volume per (dimensionless) time; all rates nonnegative.
        rd=cfg.dissociation_rate*ne*h2
        ri=cfg.ionization_rate*ne*h2
        rr=cfg.dissoc_recomb_rate*ne*ni
        rh=cfg.radical_recomb_rate*h*h
        # Correct reaction stoichiometry:
        # H2 + e -> 2H + e; H2+e -> H2+ + 2e;
        # H2+ + e -> 2H; 2H -> H2.
        # Electron power P_ext + j.E accounts external work (not predictive discharge physics).
        external=cfg.external_heating+joule_per_cell
        xfer=cfg.electron_to_gas_exchange*(ee-eg)
        increments={
            'h2':transport(h2,np.full(cfg.n,cfg.flow_speed),cfg.d_neutral,dx)-rd-ri+rh,
            'h':transport(h,np.full(cfg.n,cfg.flow_speed),cfg.d_neutral,dx)+2*rd+2*rr-2*rh,
            'ni':ni_rhs+ri-rr,
            'ne':ne_rhs+ri-rr,
            'ee':transport(ee,ve,cfg.d_e,dx)+external-DISSOC*rd-ION*ri-xfer,
            'eg':transport(eg,np.full(cfg.n,cfg.flow_speed),cfg.d_neutral,dx)
                 +(ION-DISSOC)*rr+DISSOC*rh+xfer,
        }
        for k in increments:
            s[k] += cfg.dt*increments[k]
            if np.min(s[k])<-1e-9:raise FloatingPointError(f'negative {k} at iteration {n}')
        added_energy+=cfg.dt*(cfg.external_heating*cfg.length+applied_work*cfg.length)
        now=integrals(s,cfg)
        max_atom=max(max_atom,abs(now[0]-initial[0]))
        max_charge=max(max_charge,abs(now[1]-initial[1]))
        max_energy=max(max_energy,abs(now[2]-initial[2]-added_energy))
        if (n+1)%(max(1,nstep//10))==0:
            snapshots.append({'t':(n+1)*cfg.dt,'mean_H':float(s['h'].mean()),
                              'mean_ions':float(s['ni'].mean()),
                              'mean_electrons':float(s['ne'].mean()),
                              'max_electric_field':float(np.max(np.abs(e))),
                              'total_energy':float(now[2])})
    e,phi=electric_field(s['ni'],s['ne'],cfg)
    report={'model':'dimensionless conservative 1D drift-diffusion-Poisson, NOT validated plasma',
            'settings':asdict(cfg),'steps':nstep,'max_CFL_plus_diffusion':max_cfl,
            'absolute_H_atom_inventory_error':max_atom,
            'absolute_integrated_charge_error':max_charge,
            'absolute_total_energy_including_Poisson_field_error':max_energy,
            'integrated_external_heater_and_applied_field_work':added_energy,
            'initial':initial.tolist(),'final':integrals(s,cfg).tolist(),
            'max_field':largest_field,
            'limitation':'periodic dimensionless transport, assumed rates and transport coefficients; no cathode sheath, EEDF or benchmark plasma experiment'}
    if save_path is not None:
        from pathlib import Path
        p=Path(save_path);p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(json.dumps({'report':report,'snapshots':snapshots},indent=2)+'\n')
        np.savez(p.with_suffix('.npz'),**s,E=e,potential=phi)
    return report, s, snapshots

if __name__=='__main__':
    r,_,_=run(save_path='v10/results/plasma1d.json')
    print(json.dumps({k:v for k,v in r.items() if k!='settings'},indent=2))
