#!/usr/bin/env python3
"""2-D incompressible pressure-projection MAC Navier-Stokes channel benchmark.

Periodic x, no-slip y; drive with constant x body force. Optional fixed
Brinkman resistance localized in x produces fully 2-D flow. Finite-difference
predictor includes convective terms; FFT(x) + DCT(y) pressure-Poisson solve.
Dimensionless assumptions: Newtonian incompressible constant density/viscosity.
Not a multiphase or reacting Navier-Stokes solver.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.fft import fft, ifft, dct, idct


def poisson_neumann_periodic(rhs,dx,dy):
    nx,ny=rhs.shape
    # Invert the discrete MAC div(grad) operator exactly in Fourier/DCT basis.
    a=dct(fft(rhs,axis=0),axis=1,type=2,norm='ortho')
    lx=2*(np.cos(2*np.pi*np.arange(nx)/nx)-1)/dx**2
    ly=2*(np.cos(np.pi*np.arange(ny)/ny)-1)/dy**2
    eigen=lx[:,None]+ly[None,:]
    a[0,0]=0
    eigen[0,0]=1
    p=idct(ifft(a/eigen,axis=0),axis=1,type=2,norm='ortho').real
    return p


def simulate(nx=48,ny=16,length=0.12,height=0.02,nu=5e-5,
             drive_velocity=.016,dt=.003,end=7.,porous_drag=0.,
             out='results/v3/flow', B_T=0., conductivity_S_m=7e5, density_kg_m3=7000.):
    if min(nx,ny)<6 or min(length,height,nu,drive_velocity,dt,end)<=0 or porous_drag<0 or B_T<0 or conductivity_S_m<=0 or density_kg_m3<=0:
        raise ValueError('invalid flow parameters')
    dx=length/nx; dy=height/ny
    if dt*nu*(2/dx**2+2/dy**2)>=0.4:
        raise ValueError('viscous stability: reduce dt or coarsen mesh')
    u=np.zeros((nx,ny))   # periodic x faces, y centered
    v=np.zeros((nx,ny+1)) # x centered, y faces; v=0 at top/bottom
    p=np.zeros((nx,ny))
    body=12*nu*drive_velocity/height**2
    # Low-magnetic-Re magnetic braking benchmark: B=(0,B,0), E=0.
    # Lorentz x-force = -sigma*B^2*u / rho. For insulating walls, this
    # is an idealized imposed-electric-field assumption, not general MHD.
    mhd_drag=conductivity_S_m*B_T**2/density_kg_m3
    xface=np.arange(nx)*dx
    xcell=(np.arange(nx)+.5)*dx
    ycell=(np.arange(ny)+.5)*dy
    if porous_drag>0:
        # Periodic resistance zone fixed in space, spatially nonuniform.
        zone=(xface[:,None]>0.35*length)&(xface[:,None]<0.65*length)
        mask_u=porous_drag*zone*np.ones((1,ny))
        zone_v=(xcell[:,None]>0.35*length)&(xcell[:,None]<0.65*length)
        mask_v=porous_drag*zone_v*np.ones((1,ny+1))
    else:
        mask_u=np.zeros_like(u);mask_v=np.zeros_like(v)
    mask_u+=mhd_drag
    steps=int(np.ceil(end/dt));dt=end/steps
    log=[]
    for step in range(steps):
        ub=np.concatenate([-u[:,:1],u[:,:-1]],axis=1)
        ut=np.concatenate([u[:,1:],-u[:,-1:]],axis=1)
        lapu=(np.roll(u,-1,axis=0)-2*u+np.roll(u,1,axis=0))/dx**2 + (ut-2*u+ub)/dy**2
        vatu=0.25*(v+np.roll(v,1,axis=0))[:,:-1]+0.25*(v+np.roll(v,1,axis=0))[:,1:]
        dux=np.where(u>=0,(u-np.roll(u,1,axis=0))/dx,(np.roll(u,-1,axis=0)-u)/dx)
        duy=np.where(vatu>=0,(u-ub)/dy,(ut-u)/dy)
        ust=u+dt*(nu*lapu-u*dux-vatu*duy+body-mask_u*u)
        # v is zero at both walls; ghost values only required at x-periodic faces.
        vx_left=np.roll(v,1,axis=0)
        vx_right=np.roll(v,-1,axis=0)
        lapv=(vx_left[:,1:-1]-2*v[:,1:-1]+vx_right[:,1:-1])/dx**2 + (v[:,:-2]-2*v[:,1:-1]+v[:,2:])/dy**2
        uatv=0.25*(u+np.roll(u,-1,axis=0))[:,:-1]+0.25*(u+np.roll(u,-1,axis=0))[:,1:]
        dvx=np.where(uatv>=0,(v[:,1:-1]-vx_left[:,1:-1])/dx,(vx_right[:,1:-1]-v[:,1:-1])/dx)
        dvy=np.where(v[:,1:-1]>=0,(v[:,1:-1]-v[:,:-2])/dy,(v[:,2:]-v[:,1:-1])/dy)
        vst=v.copy()
        vst[:,1:-1]+=dt*(nu*lapv-uatv*dvx-v[:,1:-1]*dvy-mask_v[:,1:-1]*v[:,1:-1])
        div=(np.roll(ust,-1,axis=0)-ust)/dx + (vst[:,1:]-vst[:,:-1])/dy
        pc=poisson_neumann_periodic(div/dt,dx,dy)
        u=ust-dt*(pc-np.roll(pc,1,axis=0))/dx
        v=vst
        v[:,1:-1]-=dt*(pc[:,1:]-pc[:,:-1])/dy
        p+=pc
        # Pressure has an arbitrary constant; recenter for stable output.
        p-=p.mean()
        if step%max(1,steps//100)==0 or step==steps-1:
            d=(np.roll(u,-1,axis=0)-u)/dx+(v[:,1:]-v[:,:-1])/dy
            log.append(((step+1)*dt,float(np.max(np.abs(d))),float(np.mean(u)),
                        float(np.max(u)),float(np.max(np.abs(v)))))
        if not np.all(np.isfinite(u)) or np.max(np.abs(u))>5:
            raise FloatingPointError('unstable MAC flow; decrease dt')
    d=(np.roll(u,-1,axis=0)-u)/dx+(v[:,1:]-v[:,:-1])/dy
    y=(np.arange(ny)+.5)*dy
    exact=body*y*(height-y)/(2*nu)
    if mhd_drag>0:
        k=np.sqrt(mhd_drag/nu)
        exact=body/mhd_drag*(1-np.cosh(k*(y-height/2))/np.cosh(k*height/2))
    err=float(np.linalg.norm(u.mean(axis=0)-exact)/max(np.linalg.norm(exact),1e-30)) if porous_drag==0 else None
    summary={'model':'2-D incompressible MAC NS periodic pressure-driven channel',
             'flow_has_convection':True,'nx':nx,'ny':ny,'end_s':end,'dt_s':dt,
             'porous_Brinkman_drag_per_s':porous_drag,
             'body_force_m_s2':body,'mean_U_m_s':float(u.mean()),
             'max_U_m_s':float(np.max(u)),
             'max_transverse_velocity_m_s':float(np.max(np.abs(v))),
             'max_divergence_1_s':float(np.max(np.abs(d))),
             'relative_l2_vs_exact_channel_profile':err,
             'vs_Poiseuille_profile_relative_l2':err if mhd_drag==0 else None,
             'magnetic_field_T':B_T,'Hartmann_number_half_channel':float(B_T*height/2*np.sqrt(conductivity_S_m/(density_kg_m3*nu))),
             'magnetic_braking_per_s':mhd_drag,
             'Re_based_on_height_and_mean_U':float(np.mean(u)*height/nu),
             'flow_boundary':'periodic in x, no-slip horizontal walls',
             'limitation':'constant-property laminar flow, prescribed B, E=0 and low magnetic Reynolds; not a plasma model or full JxB Maxwell coupling; porous resistance assumed'}
    path=Path(out);path.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path/'velocity.npz',u_xfaces_m_s=np.concatenate([u,u[:1,:]],axis=0),
                        v_yfaces_m_s=v,x_m=(np.arange(nx)+.5)*dx,y_m=(np.arange(ny)+.5)*dy,
                        p_specific_m2_s2=p)
    (path/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    np.savetxt(path/'history.csv',np.array(log),delimiter=',',comments='',
               header='time_s,max_divergence_1_s,mean_u_m_s,max_u_m_s,max_abs_v_m_s')
    print(json.dumps({k:summary[k] for k in ['mean_U_m_s','max_transverse_velocity_m_s','max_divergence_1_s','relative_l2_vs_exact_channel_profile','Hartmann_number_half_channel']},indent=2))
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--nx',type=int,default=48);p.add_argument('--ny',type=int,default=16)
    p.add_argument('--dt',type=float,default=.003);p.add_argument('--end',type=float,default=7.)
    p.add_argument('--drag',type=float,default=0.)
    p.add_argument('--B',type=float,default=0.,help='transverse imposed magnetic field [T]')
    p.add_argument('--sigma-fluid',type=float,default=7e5,help='conducting-fluid sigma [S/m]')
    p.add_argument('--rho-fluid',type=float,default=7000.,help='conducting-fluid density [kg/m3]')
    p.add_argument('--out',default='results/v3/flow')
    a=p.parse_args();simulate(nx=a.nx,ny=a.ny,dt=a.dt,end=a.end,porous_drag=a.drag,out=a.out,B_T=a.B,conductivity_S_m=a.sigma_fluid,density_kg_m3=a.rho_fluid)
