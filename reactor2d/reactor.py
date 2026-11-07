#!/usr/bin/env python3
"""2-D conservative FVM gas/solid reactive-transport with temperature-feedback MQS.

SIMPLIFICATIONS: prescribed (plug or fully developed Poiseuille) gas velocity;
no numerical Navier-Stokes solve inside this reactor module. An independent
2-D incompressible NS benchmark lives in cfd/. The interpenetrating gas and
solid continua share cells (local thermal non-equilibrium), not separate meshes.
Linear slab MQS uses axial-mean solid T, not coil geometry nor plasma MHD.
All property/kinetic parameters are illustrative, NOT verified battery data.
"""
from __future__ import annotations
import argparse
import csv
import json
import sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from reactor2d.thermo import reaction_rate, equilibrium_constant
from em.variable_induction import solve_variable, conductivity_from_temperature

SB=5.670374419e-8

def x_flux(C, ufaces, dx, diffusivity, epsilon, inlet):
    """Conservative species flux at every face, including inlet/outlet."""
    nx,ny=C.shape
    u=np.broadcast_to(ufaces,(nx+1,ny))
    flux=np.zeros((nx+1,ny))
    flux[0,:]=u[0,:]*np.where(u[0,:]>=0,inlet,C[0,:])-epsilon*diffusivity*2*(C[0,:]-inlet)/dx
    flux[1:nx,:]=u[1:nx,:]*np.where(u[1:nx,:]>=0,C[:-1,:],C[1:,:])-epsilon*diffusivity*(C[1:,:]-C[:-1,:])/dx
    flux[nx,:]=u[nx,:]*np.where(u[nx,:]>=0,C[-1,:],inlet)
    return flux

def x_heat_flux(T,ufaces,dx,k,cp,inlet):
    nx,ny=T.shape
    u=np.broadcast_to(ufaces,(nx+1,ny))
    flux=np.zeros((nx+1,ny))
    flux[0,:]=u[0,:]*cp*np.where(u[0,:]>=0,inlet,T[0,:])-2*k*(T[0,:]-inlet)/dx
    flux[1:nx,:]=u[1:nx,:]*cp*np.where(u[1:nx,:]>=0,T[:-1,:],T[1:,:])-k*(T[1:,:]-T[:-1,:])/dx
    flux[nx,:]=u[nx,:]*cp*np.where(u[nx,:]>=0,T[-1,:],inlet)
    return flux

def y_flux(C,dy,k,vfaces=None,advective_scale=1.):
    nx,ny=C.shape
    flux=np.zeros((nx,ny+1))
    flux[:,1:ny]=-k*(C[:,1:]-C[:,:-1])/dy
    if vfaces is not None:
        v=vfaces[:,1:ny]
        flux[:,1:ny]+=advective_scale*v*np.where(v>=0,C[:,:-1],C[:,1:])
    return flux

def divergence(fx,fy,dx,dy):
    return (fx[1:,:]-fx[:-1,:])/dx+(fy[:,1:]-fy[:,:-1])/dy

def ergun_pressure_drop(U, *, L=0.12, eps=0.45, dp=0.002,
                        mu=4e-5, rho=0.20):
    """Postprocessed diagnostic, NOT coupled to momentum/pressure field."""
    if min(U,L,eps,dp,mu,rho)<=0 or eps>=1: raise ValueError('invalid Ergun input')
    linear=150*mu*(1-eps)**2*U/(eps**3*dp**2)
    quadratic=1.75*rho*(1-eps)*U**2/(eps**3*dp)
    return L*(linear+quadratic)

def simulate(*, nx=48,ny=16,dt=0.003,end=12., Umean=0.016,
             h0=10000., no_reaction=False, no_em=False, feedback=True,
             velocity_model='plug', out='results/v3/reactor', em_steps=100,
             emissivity=0.20,velocity_file=None,
             kinetic_Ea=85000.,kinetic_A=1500.,kinetic_dH=-2200.,kinetic_dS=50.):
    if (nx<6 or ny<5 or min(dt,end,Umean)<=0 or h0<0 or em_steps<1
        or velocity_model not in ('plug','poiseuille') or not 0<=emissivity<=1
        or kinetic_Ea<0 or kinetic_A<0):
        raise ValueError('invalid inputs')
    L=0.12; height=0.02; dx=L/nx;dy=height/ny;dV=dx*dy
    eps=0.45; N0=150.;cp_g=450.; cp_s=1.188e6
    kg=0.008;ks=1.0;D=1e-5;hvol=600.;deltaH=kinetic_dH
    tin=850.; tamb=825.;cH2_in=15.;cH2O_in=0.
    y=(np.arange(ny)+0.5)*dy
    x=(np.arange(nx)+0.5)*dx
    if velocity_model=='plug':
        u=np.full(ny,Umean)
    else:
        eta=y/height;u=6*Umean*eta*(1-eta)
        u*=Umean/np.mean(u) # discrete mean agrees with user input
    u_faces=np.broadcast_to(u[None,:],(nx+1,ny)).copy()
    v_faces=np.zeros((nx,ny+1))
    if velocity_file is not None:
        with np.load(velocity_file) as f:
            ux=f['u_xfaces_m_s'];vy=f['v_yfaces_m_s']
            if ux.shape != (nx+1,ny) or vy.shape != (nx,ny+1):
                raise ValueError('flow mesh differs from reactive scalar mesh')
            u_faces=ux.copy();v_faces=vy.copy()
            if np.max(np.abs(v_faces[:,[0,-1]]))>1e-10:
                raise ValueError('flow must be impermeable at y walls')
        u=u_faces[:-1,:].mean(axis=0)
    advCFL=float(np.max(np.abs(u_faces))*dt/(eps*dx) + np.max(np.abs(v_faces))*dt/(eps*dy))
    diffCFL=float(dt*(D*(2/dx**2+2/dy**2)))
    heatCFL=float(dt*(kg*(2/dx**2+2/dy**2))/(eps*cp_g))
    exchCFL=float(dt*hvol/(eps*cp_g))
    if advCFL+diffCFL/eps > 0.85 or heatCFL+exchCFL>0.85:
        raise ValueError(f'unstable dt: adv={advCFL:.3g}; species={diffCFL/eps:.3g}; gas heat={heatCFL+exchCFL:.3g}')
    h2=np.zeros((nx,ny));h2o=np.zeros((nx,ny))
    Tg=np.full((nx,ny),tin);Ts=np.full((nx,ny),tin); X=np.zeros((nx,ny))
    initial_energy=float(np.sum(eps*cp_g*Tg+cp_s*Ts)*dV)
    accumulated_h2=0.;accumulated_energy=0.; cum_qem=0.; cum_rxn=0.;cum_rad=0.
    nr=int(np.ceil(end/dt));dt=end/nr
    times=[]; conversions=[];temps=[];massres=[];energyres=[];qinput=[]
    em=solve_variable(np.full(301,1e6),length=height,h0=max(h0,1e-12))
    qprofile=np.interp(y,em['y'],em['q']) if not no_em else np.zeros(ny)
    em_refreshes=0
    if no_em: qprofile[:]=0
    for step in range(nr):
        if feedback and not no_em and step%em_steps==0:
            ts_mean_y=np.mean(Ts,axis=0)
            temperature_em=np.interp(em['y'],y,ts_mean_y)
            sigma=conductivity_from_temperature(temperature_em)
            em=solve_variable(sigma,length=height,h0=h0)
            qprofile=np.interp(y,em['y'],em['q'])
            em_refreshes+=1
        fh2x=x_flux(h2,u_faces,dx,D,eps,cH2_in)
        fh2y=y_flux(h2,dy,eps*D,v_faces)
        fh2ox=x_flux(h2o,u_faces,dx,D,eps,cH2O_in)
        fh2oy=y_flux(h2o,dy,eps*D,v_faces)
        # Operator-split conservative species transport.
        h2transport=h2-dt*divergence(fh2x,fh2y,dx,dy)/eps
        h2otransport=h2o-dt*divergence(fh2ox,fh2oy,dx,dy)/eps
        if min(float(h2transport.min()),float(h2otransport.min()))<-1e-9:
            raise FloatingPointError('negative transported concentration; reduce dt')
        if no_reaction:
            rate=np.zeros_like(h2)
        else:
            rate=reaction_rate(Ts,h2transport,h2otransport,X,n0=N0,
                               preexp=kinetic_A,Ea=kinetic_Ea,
                               dh_ref=kinetic_dH,ds_ref=kinetic_dS)
            rate=np.minimum(rate,np.maximum(0,eps*h2transport/dt))
            rate=np.minimum(rate,np.maximum(0,N0*(1-X)/dt))
        h2=h2transport-dt*rate/eps
        h2o=h2otransport+dt*rate/eps
        X=X+dt*rate/N0
        # Full, face-based conservative gas sensible energy.
        fTgx=x_heat_flux(Tg,u_faces,dx,kg,cp_g,tin)
        fTgy=y_flux(Tg,dy,kg,v_faces,cp_g)
        inter=hvol*(Ts-Tg)    # equal/opposite interphase transfer
        newTg=Tg+dt*(-divergence(fTgx,fTgy,dx,dy)+inter)/(eps*cp_g)
        # Conduction in the continuous solid, boundaries insulated except radiation.
        fTsx=np.zeros((nx+1,ny))
        fTsx[1:nx,:]=-ks*(Ts[1:,:]-Ts[:-1,:])/dx
        fTsy=y_flux(Ts,dy,ks)
        radiative=np.zeros_like(Ts)
        radiative[:,0]+=emissivity*SB*(Ts[:,0]**4-tamb**4)/dy
        radiative[:,-1]+=emissivity*SB*(Ts[:,-1]**4-tamb**4)/dy
        Qem=np.broadcast_to(qprofile,(nx,ny))
        Qrxn=-deltaH*rate
        newTs=Ts+dt*(-divergence(fTsx,fTsy,dx,dy)+Qem+Qrxn-inter-radiative)/cp_s
        if min(float(newTg.min()),float(newTs.min()))<=0 or not np.all(np.isfinite(newTs)):
            raise FloatingPointError('unphysical temperatures; reduce dt')
        Tg,Ts=newTg,newTs
        # Integrated boundary flux accounting, per unit depth.
        accumulated_h2+=dt*dy*np.sum(fh2x[0,:]-fh2x[-1,:])
        boundary_g=dy*np.sum(fTgx[0,:]-fTgx[-1,:])
        accumulated_energy+=dt*(boundary_g+dV*np.sum(Qem+Qrxn-radiative))
        cum_qem+=dt*dV*np.sum(Qem)
        cum_rxn+=dt*dV*np.sum(Qrxn)
        cum_rad+=dt*dV*np.sum(radiative)
        if step%max(1,nr//100)==0 or step==nr-1:
            inventory=float(dV*np.sum(eps*h2+N0*X))
            energy_now=float(np.sum(eps*cp_g*Tg+cp_s*Ts)*dV)
            mres=accumulated_h2-inventory
            eres=accumulated_energy-(energy_now-initial_energy)
            times.append((step+1)*dt)
            conversions.append(float(X.mean()))
            temps.append(float(Ts.max()))
            massres.append(mres)
            energyres.append(eres)
            qinput.append(float(np.sum(Qem)*dV))
    final_energy=float(np.sum(eps*cp_g*Tg+cp_s*Ts)*dV)
    summary={
        'model':'2-D conservative transport, two-temperature LTNE, slab MQS induction',
        'not_implemented':'reactor NS momentum coupling, granular multiphase tracking, full coil FEM, MHD feedback, plasma, validated NiO kinetics',
        'nx':nx,'ny':ny,'dt_s':dt,'duration_s':end,'velocity_model':velocity_model if velocity_file is None else 'computed_MAC_velocity',
        'velocity_field_source':str(velocity_file) if velocity_file else 'prescribed',
        'max_transverse_velocity_m_s':float(np.max(np.abs(v_faces))),
        'average_superficial_velocity_m_s':float(u.mean()),'em_temperature_feedback':feedback,
        'em_field_recomputations':em_refreshes, 'reaction_enabled':not no_reaction,
        'em_enabled':not no_em, 'emissivity':emissivity,
        'mean_conversion':float(X.mean()),'max_solid_T_K':float(Ts.max()),
        'mean_solid_T_K':float(Ts.mean()),'outlet_H2_mol_m3':float(h2[-1,:].mean()),
        'hydrogen_balance_abs_mol_per_m':abs(float(massres[-1])),
        'hydrogen_balance_rel':abs(float(massres[-1]))/max(abs(accumulated_h2),1e-12),
        'energy_balance_abs_J_per_m':abs(float(energyres[-1])),
        'energy_balance_rel':abs(float(energyres[-1]))/max(abs(accumulated_energy),1e-12),
        'thermal_energy_delta_J_per_m':final_energy-initial_energy,
        'net_boundary_plus_sources_energy_J_per_m':accumulated_energy,
        'cumulative_EM_heat_J_per_m':cum_qem,
        'cumulative_reaction_heat_J_per_m':cum_rxn,
        'cumulative_radiation_loss_J_per_m':cum_rad,
        'Ergun_pressure_drop_Pa_diagnostic_only':ergun_pressure_drop(float(u.mean())),
        'K_equilibrium_at_mean_solid_T':float(equilibrium_constant(Ts.mean(),kinetic_dH,kinetic_dS)),
        'kinetic_A_s-1':kinetic_A,'kinetic_Ea_J_mol':kinetic_Ea,
        'thermo_dH_J_mol':kinetic_dH,'thermo_dS_J_mol_K':kinetic_dS,
        'CFL_adv':advCFL,'CFL_species':diffCFL/eps,'CFL_gas_energy_plus_exchange':heatCFL+exchCFL,
        'disclaimer':'illustrative parameters, not predictive validated battery recycling'
    }
    path=Path(out);path.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path/'fields.npz',x_m=x,y_m=y,H2=h2,H2O=h2o,
                        Tgas_K=Tg,Tsolid_K=Ts,NiO_conversion=X,qem_W_m3=qprofile,
                        u_superficial_m_s=u,u_xfaces_m_s=u_faces,v_yfaces_m_s=v_faces)
    (path/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (path/'history.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['time_s','conversion_mean','Tsolid_max_K','H_balance_mol_m','E_balance_J_m','EM_power_W_m'])
        w.writerows(zip(times,conversions,temps,massres,energyres,qinput))
    print(json.dumps({k:summary[k] for k in ['mean_conversion','max_solid_T_K','hydrogen_balance_rel','energy_balance_rel','Ergun_pressure_drop_Pa_diagnostic_only']},indent=2))
    return summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--nx',type=int,default=48);p.add_argument('--ny',type=int,default=16)
    p.add_argument('--dt',type=float,default=0.003);p.add_argument('--end',type=float,default=12.)
    p.add_argument('--velocity',type=float,default=.016);p.add_argument('--h0',type=float,default=10000.)
    p.add_argument('--no-reaction',action='store_true');p.add_argument('--no-em',action='store_true')
    p.add_argument('--fixed-sigma',action='store_true')
    p.add_argument('--velocity-model',choices=('plug','poiseuille'),default='plug')
    p.add_argument('--out',default='results/v3/reactor')
    p.add_argument('--flow-file',default=None,help='Optional velocity.npz from cfd/mac_channel.py')
    p.add_argument('--Ea',type=float,default=85000.,help='illustrative Arrhenius activation energy [J/mol]')
    p.add_argument('--A',type=float,default=1500.,help='illustrative Arrhenius prefactor [1/s]')
    a=p.parse_args()
    simulate(nx=a.nx,ny=a.ny,dt=a.dt,end=a.end,Umean=a.velocity,h0=a.h0,
             no_reaction=a.no_reaction,no_em=a.no_em,feedback=not a.fixed_sigma,
             velocity_model=a.velocity_model,out=a.out,velocity_file=a.flow_file,kinetic_Ea=a.Ea,kinetic_A=a.A)
if __name__=='__main__':main()
