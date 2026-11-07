"""Conservative 1-D weakly-compressible porous-gas transport.

State: gas species in mol/m3 of *pore gas*, T [K], p = sum(c)*R*T [Pa].
Darcy/Ergun face flux is implicitly pressure-corrected at lagged temperature.
Pressure system comes from total molar conservation; species advection uses
exactly the same pressure-corrected total molar flux. This is NOT full NS/PIMPLE.
"""
import numpy as np
from scipy.linalg import solve_banded
R = 8.31446261815324
M = np.array([.0280101,.0440095,.0280134,.00201588,.00100794])

def pressure_corrected_step(gas, Tg, source_bulk, dt, length, porosity,
                            pin, pout, diameter=.003, viscosity=3.5e-5,
                            *, inlet_fraction=None):
    gas=np.asarray(gas,float)
    ns,n=gas.shape
    if ns != 5 or n<3: raise ValueError('gas must have 5 species and >=3 cells')
    if np.min(gas)<0 or np.min(Tg)<=0 or not 0<porosity<1 or dt<=0:raise ValueError('invalid gas state')
    dz=length/n
    c=gas.sum(axis=0)
    p_old=c*R*Tg
    if inlet_fraction is None: inlet_fraction=np.array([0.,.12,.85,.03,0.])
    inlet_fraction=np.asarray(inlet_fraction,float)
    if abs(inlet_fraction.sum()-1)>1e-12 or inlet_fraction.min()<0:raise ValueError('invalid inlet mixture')
    # Freeze mobility at p_old and T_old, then update p using implicit mass balance.
    p_ext=np.r_[pin,p_old,pout]
    temperature_face=np.empty(n+1)
    temperature_face[0]=Tg[0]; temperature_face[-1]=Tg[-1]
    temperature_face[1:-1]=.5*(Tg[:-1]+Tg[1:])
    rho_cell=gas.T@M
    rho_face=np.empty(n+1)
    rho_face[0]=rho_cell[0];rho_face[-1]=rho_cell[-1]
    rho_face[1:-1]=.5*(rho_cell[:-1]+rho_cell[1:])
    dp=np.diff(p_ext)
    distance=np.full(n+1,dz);distance[[0,-1]]=dz/2
    grad=np.abs(dp/distance)
    eps=porosity
    a=150*viscosity*(1-eps)**2/(eps**3*diameter**2)
    b=1.75*rho_face*(1-eps)/(eps**3*diameter)
    # velocity magnitude for a|u|+b*u**2=|grad p|; stable as grad -> 0
    mobility=2./(a+np.sqrt(a*a+4*b*grad))
    p_face_old=np.r_[.5*(pin+p_old[0]), .5*(p_old[:-1]+p_old[1:]), .5*(p_old[-1]+pout)]
    K=(p_face_old/(R*temperature_face))*mobility/distance # mol/(m2 s Pa)
    diag=eps*dz/(R*Tg)+dt*(K[:-1]+K[1:])
    upper=-dt*K[1:-1]
    lower=-dt*K[1:-1]
    rhs=eps*dz*c+dt*dz*source_bulk.sum(axis=0)
    rhs[0]+=dt*K[0]*pin
    rhs[-1]+=dt*K[-1]*pout
    ab=np.zeros((3,n));ab[0,1:]=upper;ab[1,:]=diag;ab[2,:-1]=lower
    p_new=solve_banded((1,1),ab,rhs,check_finite=True)
    if np.min(p_new)<=0:raise FloatingPointError('negative pressure')
    # Face flux uses the identical implicit operator coefficients.
    face_flux=np.empty(n+1)
    face_flux[0]=K[0]*(pin-p_new[0]);face_flux[-1]=K[-1]*(p_new[-1]-pout)
    face_flux[1:-1]=K[1:-1]*(p_new[:-1]-p_new[1:])
    y=gas/np.maximum(c,1e-300)
    species_flux=np.zeros((ns,n+1))
    for face in range(n+1):
        J=face_flux[face]
        if face==0: yf=inlet_fraction if J>=0 else y[:,0]
        elif face==n: yf=y[:,-1] if J>=0 else inlet_fraction
        else:yf=y[:,face-1] if J>=0 else y[:,face]
        species_flux[:,face]=J*yf
    new=gas+dt*(source_bulk-(species_flux[:,1:]-species_flux[:,:-1])/dz)/eps
    if float(new.min()) < -1e-9: raise FloatingPointError(f'negative gas concentration: {new.min()}')
    new=np.maximum(new,0.)
    # Diagnostic: total species follows the implicitly corrected pressure eq.
    relative_closure=np.max(np.abs(new.sum(axis=0)-p_new/(R*Tg)))/max(np.max(new.sum(axis=0)),1e-16)
    return new,face_flux,species_flux,dict(
       pressure_Pa=p_new, face_molar_flux_mol_m2_s=face_flux,
       pressure_closure_relative=float(relative_closure),
       face_superficial_velocity_m_s=face_flux*R*temperature_face/np.maximum(p_face_old,1.),
       max_species_CFL=float(np.max(np.abs(face_flux))*dt/dz/eps/max(float(c.min()),1e-8)),
    )
