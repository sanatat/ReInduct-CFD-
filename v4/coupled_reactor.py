#!/usr/bin/env python3
"""V4: 2-D projection-NS / species / LTNE solid-gas energy / reaction model.

Momentum <-> chemistry coupling is through a local T,X-dependent Brinkman
resistance, NOT gas expansion/compressible flow. Periodic hydrodynamic x, open
reactive scalar inflow/outflow x. This is an academic *representative* model,
not a validated induction/plasma/particle-resolving industrial reactor.

Both balances are computed from integrated FINITE-VOLUME FACE fluxes, not by
clipping a balance until it looks small. Fundamental inputs are illustrative.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from cfd.mac_channel import poisson_neumann_periodic
from reactor2d.reactor import x_flux, y_flux, x_heat_flux, divergence
from reactor2d.thermo import reaction_rate
from em.variable_induction import solve_variable, conductivity_from_temperature
from v5.em2d import run as solve_coil_2d
from v5.thermochemistry import thermo as reference_thermo
from v5.mhd2d import motional_lorentz

SIGMA_SB=5.670374419e-8

def run(nx=24,ny=10,dt=0.002,end=3.0,back_coupling=True,reaction=True,em=True,
        out='v4/results/coupled', beta_drag=0.008, em_backend='1d', thermo_backend='illustrative', mhd_sigma_g=0.):
    if nx<8 or ny<6 or min(dt,end)<=0 or beta_drag<0: raise ValueError('invalid arguments')
    if em_backend not in ('1d','2d'):raise ValueError('unsupported EM backend')
    if thermo_backend not in ('illustrative','tabulated'):raise ValueError('unsupported thermo backend')
    if mhd_sigma_g<0 or (mhd_sigma_g>0 and (not em or em_backend!='2d')):
        raise ValueError('MHD proxy requires 2D induction EM enabled')
    length=0.12; height=0.02;dx=length/nx;dy=height/ny;vol=dx*dy
    nu=5.0e-5;eps=0.45;D=1e-5; cpg=450.;cps=1.188e6; kg=0.008;ks=1.
    hgs=600.; ns0=150.; h2in=15.;Tin=850.;Tenv=825.; dH=-2200.
    # Drag scale is illustrative [s^-1], not an independently calibrated Ergun fit.
    nominal_drag=0.8
    drive=12*nu*0.016/height**2
    nsteps=int(np.ceil(end/dt));dt=end/nsteps
    visc_cfl=dt*nu*(2/dx**2+2/dy**2)
    if visc_cfl>.4: raise ValueError('viscous CFL too large')
    ux=np.zeros((nx,ny)); vy=np.zeros((nx,ny+1));p=np.zeros((nx,ny))
    H2=np.zeros((nx,ny));H2O=np.zeros((nx,ny));X=np.zeros((nx,ny))
    Tg=np.full((nx,ny),Tin);Ts=Tg.copy()
    Einit=np.sum(eps*cpg*Tg+cps*Ts)*vol
    matter_in=0.;energy_in=0.;reacted=0.; accumulated_em=0.;cum_rad=0.
    times=[];conversion=[];tmax=[];uavg=[];balance_mass=[];balance_energy=[]; div_hist=[]
    emfield=None; qprof=np.zeros(ny); qmap=np.zeros((nx,ny))
    bxmap=np.zeros((nx,ny),complex);bymap=np.zeros((nx,ny),complex)
    cumulative_mhd_mech_energy=0.
    y=(np.arange(ny)+.5)*dy
    for step in range(nsteps):
        # Real two-way feedback: T_s and conversion X -> drag field -> U,V ->
        # reactant supply and advection -> changes in T_s and X.
        zone=((np.arange(nx)+.5)*dx>=.30*length)&((np.arange(nx)+.5)*dx<=.70*length)
        zone_ux=((np.arange(nx)*dx>=.30*length)&(np.arange(nx)*dx<=.70*length))[:,None]
        if back_coupling:
            drag=nominal_drag*np.maximum(0.2,1+beta_drag*(Ts-Tin))*(1-.3*X)
        else:
            drag=np.full_like(Ts,nominal_drag)
        du=zone_ux*np.roll((drag+np.roll(drag,1,axis=0))*.5,0,axis=0)
        dv=zone[:,None]*np.concatenate((drag[:,:1],.5*(drag[:,:-1]+drag[:,1:]),drag[:,-1:]),axis=1)
        # MAC incompressible NS predictor (periodic axial velocity, wall no-slip).
        ub=np.concatenate((-ux[:,:1],ux[:,:-1]),axis=1)
        ut=np.concatenate((ux[:,1:],-ux[:,-1:]),axis=1)
        Lu=(np.roll(ux,-1,0)-2*ux+np.roll(ux,1,0))/dx**2+(ut-2*ux+ub)/dy**2
        vatu=.25*(vy+np.roll(vy,1,0))[:,:-1]+.25*(vy+np.roll(vy,1,0))[:,1:]
        dudx=np.where(ux>=0,(ux-np.roll(ux,1,0))/dx,(np.roll(ux,-1,0)-ux)/dx)
        dudy=np.where(vatu>=0,(ux-ub)/dy,(ut-ux)/dy)
        us=ux+dt*(nu*Lu-ux*dudx-vatu*dudy+drive-du*ux)
        Lvy=(np.roll(vy,-1,0)-2*vy+np.roll(vy,1,0))[:,1:-1]/dx**2+(vy[:,:-2]-2*vy[:,1:-1]+vy[:,2:])/dy**2
        uatv=.25*(ux+np.roll(ux,-1,0))[:,:-1]+.25*(ux+np.roll(ux,-1,0))[:,1:]
        dvdx=np.where(uatv>=0,(vy[:,1:-1]-np.roll(vy,1,0)[:,1:-1])/dx,(np.roll(vy,-1,0)[:,1:-1]-vy[:,1:-1])/dx)
        dvdy=np.where(vy[:,1:-1]>=0,(vy[:,1:-1]-vy[:,:-2])/dy,(vy[:,2:]-vy[:,1:-1])/dy)
        vs=vy.copy();vs[:,1:-1]+=dt*(nu*Lvy-uatv*dvdx-vy[:,1:-1]*dvdy-dv[:,1:-1]*vy[:,1:-1])
        rhs=(np.roll(us,-1,0)-us)/dx+(vs[:,1:]-vs[:,:-1])/dy
        dp=poisson_neumann_periodic(rhs/dt,dx,dy)
        ux=us-dt*(dp-np.roll(dp,1,0))/dx
        vy=vs;vy[:,1:-1]-=dt*(dp[:,1:]-dp[:,:-1])/dy
        p+=dp;p-=p.mean()
        div=(np.roll(ux,-1,0)-ux)/dx+(vy[:,1:]-vy[:,:-1])/dy
        maxdiv=float(np.max(np.abs(div)))
        uf=np.concatenate([ux,ux[:1,:]],axis=0)
        # Density is constant in momentum and molar concentrations are reported
        # per volume of gas; species fluxes use superficial velocity.
        cfl=dt*(np.max(np.abs(uf))/dx+np.max(np.abs(vy))/dy)/eps
        diff=dt*2*D*(1/dx**2+1/dy**2)
        if cfl+diff>.75:raise FloatingPointError('species CFL exceeded')
        fhx=x_flux(H2,uf,dx,D,eps,h2in);fhy=y_flux(H2,dy,eps*D,vy)
        fwx=x_flux(H2O,uf,dx,D,eps,0.);fwy=y_flux(H2O,dy,eps*D,vy)
        hc=H2-dt*divergence(fhx,fhy,dx,dy)/eps
        wc=H2O-dt*divergence(fwx,fwy,dx,dy)/eps
        if min(hc.min(),wc.min()) < -1e-9:raise FloatingPointError('negative species before reaction')
        if reaction:
            if thermo_backend=='tabulated':
                # Curated NiO+H2 equilibrium constants and reaction enthalpy,
                # not a calibration of the empirical volumetric kinetic prefactor.
                props=reference_thermo(Ts)
                K=props['Kp']
                driving=np.maximum(hc-wc/K,0.)/15.
                rate=ns0*1500.*np.exp(-85000./(8.314462618*Ts))*driving*np.maximum(1-X,0.)
            else:
                rate=reaction_rate(Ts,hc,wc,X,n0=ns0,preexp=1500.,Ea=85000.,dh_ref=-2200.,ds_ref=50.)
            rate=np.minimum(rate,np.minimum(eps*np.maximum(hc,0.)/dt,ns0*(1-X)/dt))
        else:rate=np.zeros_like(hc)
        H2=hc-dt*rate/eps;H2O=wc+dt*rate/eps;X+=dt*rate/ns0
        # Induction field periodically refreshed by solid temperature (1-D MQS slab)
        if em and (step%100==0):
            if em_backend=='2d':
                # 2-D vector potential with temperature-dependent bed sigma.
                # q in W/m3 returned on (ny,nx) grid; reactor is (nx,ny).
                _, em2fields = solve_coil_2d(nx=nx,ny=ny,temperature=Ts.T,write=False)
                qmap=em2fields['q'].T.copy()
                bxmap=em2fields['Bx'].T.copy();bymap=em2fields['By'].T.copy()
            else:
                if emfield is None: syg=np.full(151,1e6)
                else:syg=conductivity_from_temperature(np.interp(emfield['y'],y,Ts.mean(axis=0)))
                emfield=solve_variable(syg,length=height,h0=10000.)
                qprof=np.interp(y,emfield['y'],emfield['q'])
        qem=qmap if em_backend=='2d' else np.broadcast_to(qprof,(nx,ny))
        mhd_thermal=np.zeros((nx,ny))
        if mhd_sigma_g>0:
            # Weakly conducting gas demonstration: motional JxB source, then
            # reproject to enforce incompressibility after explicit forcing.
            uc=.5*(ux+np.roll(ux,-1,axis=0))
            vc=.5*(vy[:,:-1]+vy[:,1:])
            fx,fy,mech=motional_lorentz(uc,vc,bxmap,bymap,sigma_g=mhd_sigma_g,rho_g=.10)
            cumulative_mhd_mech_energy+=dt*vol*np.sum(mech)
            mhd_thermal=-mech   # motional ohmic dissipation [W/m3] returned to gas
            ux+=.5*dt*(fx+np.roll(fx,1,axis=0))
            vy[:,1:-1]+=.5*dt*(fy[:,:-1]+fy[:,1:])
            divf=(np.roll(ux,-1,0)-ux)/dx+(vy[:,1:]-vy[:,:-1])/dy
            correction=poisson_neumann_periodic(divf/dt,dx,dy)
            ux-=dt*(correction-np.roll(correction,1,0))/dx
            vy[:,1:-1]-=dt*(correction[:,1:]-correction[:,:-1])/dy
            div=(np.roll(ux,-1,0)-ux)/dx+(vy[:,1:]-vy[:,:-1])/dy
            maxdiv=float(np.max(np.abs(div)))
            uf=np.concatenate([ux,ux[:1,:]],axis=0)
            cfl=dt*(np.max(np.abs(uf))/dx+np.max(np.abs(vy))/dy)/eps
            if cfl+diff>.75:raise FloatingPointError('species CFL exceeded after MHD projection')
            # qEM = induction heat; MHD motional heat is accounted separately.

        # Gas convection must be treated with same superficial flux as species;
        # heat flux uses cpg and effective conduction coefficient kg.
        fgx=x_heat_flux(Tg,uf,dx,kg,cpg,Tin)
        fgy=y_flux(Tg,dy,kg,vy,cpg)
        exch=hgs*(Ts-Tg)
        Tg_new=Tg+dt*(-divergence(fgx,fgy,dx,dy)+exch+mhd_thermal)/(eps*cpg)
        fsx=np.zeros((nx+1,ny));fsx[1:-1]=-ks*(Ts[1:]-Ts[:-1])/dx
        fsy=y_flux(Ts,dy,ks)
        rad=np.zeros_like(Ts)
        rad[:,0]+=.2*SIGMA_SB*(Ts[:,0]**4-Tenv**4)/dy
        rad[:,-1]+=.2*SIGMA_SB*(Ts[:,-1]**4-Tenv**4)/dy
        rheat=(-reference_thermo(Ts)['delta_H_Jmol']*rate if thermo_backend=='tabulated' else -dH*rate)
        Ts_new=Ts+dt*(-divergence(fsx,fsy,dx,dy)+qem+rheat-exch-rad)/cps
        if min(Tg_new.min(),Ts_new.min())<=0 or not np.all(np.isfinite(Tg_new)):
            raise FloatingPointError('unphysical temperature')
        Tg=Tg_new;Ts=Ts_new
        # Integrated conservation: H2 inventory + Ni-produced H atoms.
        matter_in+=dt*dy*np.sum(fhx[0,:]-fhx[-1,:])
        boundary_energy=dy*np.sum(fgx[0,:]-fgx[-1,:])
        energy_in+=dt*(boundary_energy+vol*np.sum(qem+rheat-rad+mhd_thermal))
        reacted+=dt*vol*np.sum(rate)
        accumulated_em+=dt*vol*np.sum(qem)
        cum_rad+=dt*vol*np.sum(rad)
        if step%max(1,nsteps//70)==0 or step==nsteps-1:
            mres=matter_in-vol*np.sum(eps*H2+ns0*X)
            eres=energy_in-(vol*np.sum(eps*cpg*Tg+cps*Ts)-Einit)
            times.append((step+1)*dt);conversion.append(float(X.mean()));tmax.append(float(Ts.max()))
            uavg.append(float(ux.mean()));balance_mass.append(float(mres));balance_energy.append(float(eres));div_hist.append(maxdiv)
        if np.max(np.abs(ux))>3 or np.any(X>1+1e-10) or np.any(H2< -1e-9):
            raise FloatingPointError('unstable/negative conversion or flow')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out/'fields.npz',H2=H2,H2O=H2O,X=X,Tg=Tg,Ts=Ts,Ux=ux,Vy=vy,
                        qem=qem if em_backend=='2d' else qprof,x_m=(np.arange(nx)+.5)*dx,y_m=y)
    data=np.array([times,conversion,tmax,uavg,balance_mass,balance_energy,div_hist]).T
    np.savetxt(out/'history.csv',data,delimiter=',',comments='',header='time_s,mean_conversion,max_Tsolid_K,mean_U_m_s,H_balance_mol_per_m,energy_balance_J_per_m,max_divergence_per_s')
    summary={'model':'V5 EM-coupled 2D incompressible MAC + FVM reactive LTNE model',
      'em_backend':em_backend,'thermo_backend':thermo_backend, 'coupling_mode':'two-way via T,X-dependent porous resistance and advective supply' if back_coupling else 'one-way fixed Darcy drag',
      'experimentally_validated':False,'openfoam_executed':False,'niO_surrogate':True,
      'nx':nx,'ny':ny,'dt_s':dt,'end_s':end,'mean_conversion':float(X.mean()),
      'mean_velocity_m_s':float(ux.mean()),'max_solid_temperature_K':float(Ts.max()),
      'max_divergence_1_s':maxdiv,
      'last_joule_heat_source_mean_W_m3':float(np.mean(qem)),
      'mhd_gas_sigma_S_m':float(mhd_sigma_g),
      'cumulative_mhd_motional_mechanical_energy_J_per_m':float(cumulative_mhd_mech_energy),'h2_balance_abs_mol_per_m':abs(float(balance_mass[-1])),
      'h2_balance_relative':abs(float(balance_mass[-1]))/max(abs(matter_in),1e-12),
      'gas_plus_solid_energy_balance_abs_J_per_m':abs(float(balance_energy[-1])),
      'energy_balance_relative':abs(float(balance_energy[-1]))/max(abs(energy_in),1e-12),
      'max_h2_nonnegative':bool(np.min(H2)>=-1e-12),
      'max_conversion_within_one':bool(np.max(X)<=1+1e-12),
      'total_consumed_H2_mol_per_m':float(reacted),'cumulative_induction_energy_J_per_m':float(accumulated_em),
      'cumulative_radiation_loss_J_per_m':float(cum_rad),
      'limitations':'periodic-x incompressible flow; open-x species; fixed density, simplified conductivity/kinetics, no resolved particle motion, 2D long-conductor induction only (not full physical coil), no plasma discharge, no experimentally validated reactor'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    return summary

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--nx',type=int,default=24);ap.add_argument('--ny',type=int,default=10)
    ap.add_argument('--dt',type=float,default=.002);ap.add_argument('--end',type=float,default=3.)
    ap.add_argument('--no-feedback',action='store_true');ap.add_argument('--no-reaction',action='store_true');ap.add_argument('--no-em',action='store_true')
    ap.add_argument('--out',default='v4/results/coupled');ap.add_argument('--em-backend',choices=['1d','2d'],default='1d');ap.add_argument('--thermo-backend',choices=['illustrative','tabulated'],default='illustrative');ap.add_argument('--mhd-sigma-g',type=float,default=0.);args=ap.parse_args()
    run(nx=args.nx,ny=args.ny,dt=args.dt,end=args.end,back_coupling=not args.no_feedback,
        reaction=not args.no_reaction,em=not args.no_em,out=args.out,em_backend=args.em_backend,thermo_backend=args.thermo_backend,mhd_sigma_g=args.mhd_sigma_g)
