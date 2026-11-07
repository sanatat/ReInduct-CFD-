"""Axisymmetric finite-volume gas-solid fixed-bed training reactor, powered by EM v6.

Solves FV species and 2-temperature energy in cylindrical (z,r) coordinates.
Axial superficial velocity is quasi-steady pressure-driven Ergun, NOT Navier-Stokes.
Gas/solid separately retained temperature, but solid particles are NOT resolved.
NiO+H2 is a kinetic surrogate, not Li-ion black-mass experimental process.
Strict integral mass/energy ledgers include boundary fluxes and volumetric sources.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import RegularGridInterpolator
from v6.axisymmetric_em import Grid,run as em_run,sample_heating
from v5.thermochemistry import thermo

R_G=8.314462618

def geometry(nz=24,nr=10,z0=.040,z1=.100,rad=.025):
    ze=np.linspace(z0,z1,nz+1);re=np.linspace(0,rad,nr+1)
    z=.5*(ze[:-1]+ze[1:]);r=.5*(re[:-1]+re[1:])
    dz=ze[1]-ze[0];dr=re[1]-re[0]
    A=np.pi*(re[1:]**2-re[:-1]**2)
    V=dz*A[None,:]*np.ones((nz,1))
    return {'z':z,'r':r,'z_edges':ze,'r_edges':re,'dz':dz,'dr':dr,'A':A,'V':V,'nz':nz,'nr':nr}

def ergun_velocity(delta_p,Ts,X,g,epsilon=.45,dp0=.003):
    """Solve integral Ergun pressure-loss equation for constant superficial u>0.
    Gas density fixed 0.15 kg/m3; viscosity grows mildly with T.
    This is pressure-driven Darcy/Ergun closure, not a compressible CFD field.
    """
    mu=3.2e-5*(np.maximum(Ts,450)/850.)**.7
    dp=dp0*(1+.12*np.clip(X,0,1))
    rho_g=.15
    a=150*(1-epsilon)**2*mu/(epsilon**3*dp**2)
    b=1.75*(1-epsilon)*rho_g/(epsilon**3*dp)
    # Area-average cell resistances and axial integration.
    a_int=np.sum(np.mean(a,axis=1)*g['dz'])
    b_int=np.sum(np.mean(b,axis=1)*g['dz'])
    u=(-a_int+np.sqrt(a_int*a_int+4*b_int*delta_p))/(2*b_int)
    pressure_profile=delta_p-np.cumsum((np.mean(a,axis=1)*u+np.mean(b,axis=1)*u*u)*g['dz'])
    return float(u),float(a_int),float(b_int),pressure_profile

def fluxes(field,g,u=0.,inlet=0.,gamma=0.):
    """Outgoing fluxes at z/r faces in scalar units per second; cylindrical areas."""
    nz,nr=g['nz'],g['nr'];A=g['A'];dz=g['dz'];dr=g['dr'];re=g['r_edges']
    zf=np.zeros((nz+1,nr));rf=np.zeros((nz,nr+1))
    if u:
        zf[0]=u*A*inlet
        zf[1:]=u*A[None,:]*field
    if gamma:
        zf[1:-1] += -gamma*A[None,:]*(field[1:]-field[:-1])/dz
        rf[:,1:-1]=-gamma*(2*np.pi*re[None,1:-1]*dz)*(field[:,1:]-field[:,:-1])/dr
    return zf,rf

def divergence_rate(zf,rf):
    return -(zf[1:]-zf[:-1]+rf[:,1:]-rf[:,:-1])

def build_em_temperature_field(Ts,g,em_grid:Grid):
    interp=RegularGridInterpolator((g['z'],g['r']),Ts,bounds_error=False,fill_value=850.)
    zz,rr=em_grid.coords()
    return interp(np.stack([zz.ravel(),rr.ravel()],axis=-1)).reshape(em_grid.nz,em_grid.nr)

def run(nz=24,nr=10,end_s=2.0,dt_s=.004,pressure_Pa=6.,induction=True,reaction=True,
        heat_exchange_W_m3K=600.,em_feedback=True,em_nr=36,em_nz=58,
        out='v6/results/axisymmetric_bed',write=True):
    g=geometry(nz,nr);V=g['V'];epsilon=.45
    Cpg=800. # J/m^3 gas/K, illustrative reference
    Cps=1.05e6 # J/m^3 bed/K effective solid heat capacity
    H2_in=11.; Tg_in=825.
    H2=np.full((nz,nr),H2_in);H2O=np.zeros_like(H2)
    Tg=np.full_like(H2,825.);Ts=np.full_like(H2,850.);X=np.zeros_like(H2)
    solid_site_density=450. # mol NiO/m3 bulk: fictive reactive fraction
    mu_h2=epsilon*2e-5 # effective Fick diffusivity of hydrogen [m2/s]
    k_g=.025;k_s=1.0 # W/m/K
    em_grid=Grid(nr=em_nr,nz=em_nz)
    qem=np.zeros_like(Ts)
    if induction:
        _,field=em_run(em_grid)
        qem=sample_heating(field,g['z'],g['r'])
    nsteps=int(round(end_s/dt_s))
    if nsteps<=0 or abs(nsteps*dt_s-end_s)>1e-9:raise ValueError('end_s must be a multiple of dt_s')
    E0=np.sum((epsilon*Cpg*Tg+Cps*Ts)*V)
    H0=float(np.sum((epsilon*H2+solid_site_density*X)*V))
    H2O0=float(np.sum(epsilon*H2O*V))
    ext_mass=0.; ext_h2o=0.;ext_energy=0.;lastu=0.
    history=[]; q_power_prev=0.
    for step in range(nsteps):
        if induction and em_feedback and step>0 and step % 80==0:
            _,field=em_run(em_grid,bed_temperature_K=build_em_temperature_field(Ts,g,em_grid))
            qem=sample_heating(field,g['z'],g['r'])
        u,aint,bint,press=ergun_velocity(pressure_Pa,Ts,X,g,epsilon)
        # First compute conservative advection/diffusion.
        h2z,h2r=fluxes(H2,g,u,H2_in,mu_h2)
        wtz,wtr=fluxes(H2O,g,u,0.,mu_h2)
        hg=H2+dt_s*divergence_rate(h2z,h2r)/(epsilon*V)
        wg=H2O+dt_s*divergence_rate(wtz,wtr)/(epsilon*V)
        if min(float(hg.min()),float(wg.min())) < -1e-9:
            raise FloatingPointError('gas negative before kinetic source; reduce dt')
        if reaction:
            K=thermo(Ts)['Kp']
            aff=np.maximum(hg-wg/K,0.)
            arr=750.*np.exp(-48000./(R_G*Ts))
            rate=epsilon*arr*aff*np.maximum(1.-X,0.)
            rate=np.minimum(rate,np.minimum(epsilon*np.maximum(hg,0.)/dt_s,
                                            solid_site_density*np.maximum(1.-X,0.)/dt_s))
            dH=thermo(Ts)['delta_H_Jmol']
        else:
            rate=np.zeros_like(H2);dH=np.zeros_like(Ts)
        H2=hg-dt_s*rate/epsilon
        H2O=wg+dt_s*rate/epsilon
        X+=dt_s*rate/solid_site_density
        if np.any(H2 < -1e-9) or np.any(H2O < -1e-9) or np.any(X>1+1e-9):
            raise FloatingPointError('unphysical mass field')
        # Energy is deposited into solid by EM and reduction reaction.
        tgz,tgr=fluxes(Tg,g,u*Cpg,Tg_in,k_g)
        tsz,tsr=fluxes(Ts,g,0.,0.,k_s)
        hgs=heat_exchange_W_m3K*(Ts-Tg)
        qreaction=-dH*rate
        Tg=Tg+dt_s*(divergence_rate(tgz,tgr)/V+hgs)/(epsilon*Cpg)
        Ts=Ts+dt_s*(divergence_rate(tsz,tsr)/V-hgs+qem+qreaction)/Cps
        if min(Tg.min(),Ts.min())<0 or not np.all(np.isfinite(Tg)):
            raise FloatingPointError('unphysical temperature')
        mass_flux_in=float(np.sum(h2z[0]-h2z[-1])); product_flux_in=float(np.sum(wtz[0]-wtz[-1]))
        energy_flux_in=float(np.sum(tgz[0]-tgz[-1]))
        ext_mass+=dt_s*mass_flux_in;ext_h2o+=dt_s*product_flux_in
        ext_energy+=dt_s*(energy_flux_in+float(np.sum((qreaction+qem)*V)))
        Mnow=float(np.sum((epsilon*H2+solid_site_density*X)*V))
        Wnow=float(np.sum((epsilon*H2O-solid_site_density*X)*V))
        Enow=float(np.sum((epsilon*Cpg*Tg+Cps*Ts)*V))
        r_mass=Mnow-H0-ext_mass
        r_product=Wnow-H2O0-ext_h2o
        r_energy=Enow-E0-ext_energy
        if step==0 or (step+1)%max(1,nsteps//30)==0 or step==nsteps-1:
            history.append([dt_s*(step+1),float(np.mean(X)),float(np.max(Ts)),u,
                            r_mass,r_product,r_energy,float(np.sum(qem*V))])
        lastu=u
    report={'model':'Axisymmetric FV plug/Ergun-flow LTNE reacting fixed bed with EM feedback (not Navier-Stokes CFD)',
            'mesh_nz_nr':[nz,nr], 'steps':nsteps,'end_s':end_s,'dt_s':dt_s,
            'em_2d_axisymmetric':induction,'thermo_NiO_H2':reaction,
            'temperature_feedback_on_electrical_conductivity':em_feedback and induction,
            'mean_NiO_conversion':float(np.mean(X)),
            'mean_superficial_U_m_s':lastu,'estimated_pressure_drop_Pa':pressure_Pa,
            'max_solid_temperature_K':float(np.max(Ts)),
            'final_inductive_power_W':float(np.sum(qem*V)),
            'hydrogen_inventory_balance_residual_mol':float(r_mass),
            'water_plus_solid_inventory_balance_residual_mol':float(r_product),
            'sensible_energy_ledger_residual_J':float(r_energy),
            'max_X':float(X.max()),'min_H2_mol_m3':float(H2.min()),
            'qualifications':['surrogate nickel oxide, not NMC/LFP black-mass rate',
                              'superficial flow is solved via Ergun not Navier-Stokes',
                              'illustrative constants; pressure/compressibility approximate',
                              'no resolved particles, plasma chemistry, or reactor experimental validation']}
    if write:
        dest=Path(out);dest.mkdir(parents=True,exist_ok=True)
        (dest/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
        np.savez_compressed(dest/'fields.npz',z=g['z'],r=g['r'],X=X,H2=H2,H2O=H2O,Tg=Tg,Ts=Ts,qEM=qem)
        np.savetxt(dest/'history.csv',np.asarray(history),delimiter=',',comments='',
          header='time_s,mean_NiO_conversion,max_Tsolid_K,superficial_U_m_s,H2_balance_mol,H2O_balance_mol,sensible_E_balance_J,EM_power_W')
    return report
