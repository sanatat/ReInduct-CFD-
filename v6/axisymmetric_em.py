"""Verified axisymmetric, frequency-domain magnetoquasistatic induction model.

Geometry (r,z); azimuthal current J_theta, magnetic potential A_theta=r*psi.
Equation: -(d2 psi/dr2+3/r*dpsi/dr+d2 psi/dz2)
          + 1j*omega*mu*sigma*psi = mu*J_coil_theta/r.
Axis regularity psi_r(0)=0; A=0 at far radial and axial boundaries.
Finite-volume radial operator integrates r**3 weight. Not a plasma model.
All conductivities and coil drives are *illustrative*, not measured battery data.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from scipy.interpolate import RegularGridInterpolator

MU0 = 4*np.pi*1e-7

@dataclass(frozen=True)
class Grid:
    nr:int=40
    nz:int=64
    radius_m:float=0.065
    length_m:float=0.14
    @property
    def dr(self):return self.radius_m/self.nr
    @property
    def dz(self):return self.length_m/self.nz
    @property
    def r(self):return (np.arange(self.nr)+.5)*self.dr
    @property
    def z(self):return (np.arange(self.nz)+.5)*self.dz
    def coords(self):return np.meshgrid(self.z,self.r,indexing='ij')

def operator(grid:Grid,sigma:np.ndarray,omega:float):
    """Sparse positive elliptic operator for psi, on z/r cell ordering."""
    nr,nz=grid.nr,grid.nz;dr,dz=grid.dr,grid.dz
    sig=np.asarray(sigma,dtype=float)
    if sig.shape!=(nz,nr) or np.any(sig<0):raise ValueError('sigma shape/nonnegative required')
    mat=lil_matrix((nr*nz,nr*nz),dtype=complex)
    for j in range(nz):
        for i in range(nr):
            p=j*nr+i
            rl=i*dr;rh=(i+1)*dr
            vol=(rh**4-rl**4)/4.0
            # effective volume-based radial face stiffness
            west=rl**3/(dr*vol)
            east=rh**3/(dr*vol)
            diag=1j*omega*MU0*sig[j,i]
            if i>0:mat[p,p-1]=-west;diag+=west
            if i<nr-1:mat[p,p+1]=-east;diag+=east
            else:diag+=2*east   # psi(R)=0, half-cell boundary distance
            if j>0:mat[p,p-nr]=-1/dz**2;diag+=1/dz**2
            else:diag+=2/dz**2
            if j<nz-1:mat[p,p+nr]=-1/dz**2;diag+=1/dz**2
            else:diag+=2/dz**2
            mat[p,p]=diag
    return mat.tocsr()

def manufactured(grid:Grid=Grid(), frequency_hz:float=1200., sigma_Sm:float=1e4):
    """Continuous manufactured exact psi=cos(pi*r/2R)*sin(pi*z/L)."""
    z,r=grid.coords(); kr=np.pi/(2*grid.radius_m);kz=np.pi/grid.length_m
    exact=np.cos(kr*r)*np.sin(kz*z)*(1+.15j)
    continuous=(-((kr**2+kz**2)*np.cos(kr*r) +3*kr*np.sin(kr*r)/r)*np.sin(kz*z)*(1+.15j))
    # PDE source is -laplacian4(psi) + i omega mu sigma psi
    source=-continuous+1j*2*np.pi*frequency_hz*MU0*sigma_Sm*exact
    K=operator(grid,np.full((grid.nz,grid.nr),sigma_Sm),2*np.pi*frequency_hz)
    solved=spsolve(K,source.ravel()).reshape(grid.nz,grid.nr)
    return {'grid':[grid.nz,grid.nr],'relative_L2':float(np.linalg.norm(solved-exact)/np.linalg.norm(exact)),
            'algebraic_relative_residual':float(np.linalg.norm(K@solved.ravel()-source.ravel())/np.linalg.norm(source)),
            'exact_psi':exact,'numerical_psi':solved}

def run(grid:Grid=Grid(),frequency_hz:float=5000., coil_current_A_m2:float=1.4e7,
        conductivity_bed_Sm:float=8e4,bed_temperature_K=None, write=False,out='v6/results/axisymmetric'):
    z,r=grid.coords()
    # Finite-thickness azimuthal coil around a finite cylindrical bed.
    # Magnetic vector potential is axisymmetric: this is closer to a real coil than v5 planar legs.
    zr=.070
    coilshape=np.exp(-.5*((r-.036)/.004)**2)*np.exp(-.5*((z-zr)/.034)**2)
    sigma_r=(.5*(1-np.tanh((r-.027)/.0015)))
    sigma_z=(.5*(1+np.tanh((z-.038)/.002)))*(.5*(1-np.tanh((z-.102)/.002)))
    sigma=conductivity_bed_Sm*sigma_r*sigma_z
    if bed_temperature_K is not None:
        tmp=np.asarray(bed_temperature_K,dtype=float)
        if tmp.shape != sigma.shape:raise ValueError('temperature grid mismatch')
        sigma/=np.clip(1+.002*(tmp-850),.25,8.)
    j_coil=coil_current_A_m2*coilshape
    K=operator(grid,sigma,2*np.pi*frequency_hz)
    psi=spsolve(K,(MU0*j_coil/np.maximum(r,1e-9)).ravel()).reshape(grid.nz,grid.nr)
    Atheta=r*psi; Etheta=-1j*2*np.pi*frequency_hz*Atheta
    jed=sigma*Etheta
    q=.5*sigma*np.abs(Etheta)**2  # time-averaged W/m^3
    # B_r=-dA_theta/dz; B_z=(1/r)*d(r A_theta)/dr
    dAdz,dAdr=np.gradient(Atheta,grid.dz,grid.dr)
    Br=-dAdz;Bz=2*psi+r*np.gradient(psi,grid.dr,axis=1)
    vol=np.pi*(((np.arange(grid.nr)+1)*grid.dr)**2-(np.arange(grid.nr)*grid.dr)**2)*grid.dz
    power=float(np.sum(q*vol[None,:]))
    source=(MU0*j_coil/np.maximum(r,1e-9)).ravel()
    residual=float(np.linalg.norm(K@psi.ravel()-source)/np.linalg.norm(source)) if np.linalg.norm(source)>0 else float(np.linalg.norm(K@psi.ravel()))
    summary={'model':'2D axisymmetric azimuthal A_theta=r*psi MQS, cylindrical coil band, no plasma',
             'mesh_nz_nr':[grid.nz,grid.nr], 'frequency_Hz':frequency_hz,
             'drive_Jtheta_A_m2':coil_current_A_m2,'bed_reference_sigma_S_m':conductivity_bed_Sm,
             'joule_power_W':power,'mean_q_bed_W_m3':float(q[(sigma_r*sigma_z)>.5].mean()),
             'peak_B_T':float(np.max(np.hypot(np.abs(Br),np.abs(Bz)))),
             'algebraic_relative_residual':residual,'temperature_feedback':bed_temperature_K is not None,
             'coil_type':'axisymmetric finite Gaussian azimuthal current band; not a CAD winding',
             'measured_properties':False,'experimental_validation':False}
    fields={'r':grid.r,'z':grid.z,'q':q,'sigma':sigma,'Atheta':Atheta,'Br':Br,'Bz':Bz,'Jeddy':jed}
    if write:
        dest=Path(out);dest.mkdir(parents=True,exist_ok=True)
        (dest/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        np.savez_compressed(dest/'fields.npz',**{k: np.real(v) if k not in ('Atheta','Br','Bz','Jeddy') else v for k,v in fields.items()})
    return summary,fields

def sample_heating(fields,z_points,r_points):
    interp=RegularGridInterpolator((fields['z'],fields['r']),fields['q'],bounds_error=False,fill_value=0.)
    zz,rr=np.meshgrid(z_points,r_points,indexing='ij')
    return interp(np.stack([zz.ravel(),rr.ravel()],axis=1)).reshape(len(z_points),len(r_points))
