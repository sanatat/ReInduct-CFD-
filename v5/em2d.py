"""2-D time-harmonic magnetic vector-potential (A_z) induction demonstrator.

Phasors use e^(i omega t): E_z=-i omega A_z, J_eddy=sigma E_z,
(-laplace + i*omega*mu*sigma) A_z=mu J_coil_z.
The geometry is *translationally invariant straight conductors*, NOT a
real axisymmetric coil or plasma discharge. Effective conductivity, coil
currents, and properties are illustrative rather than measured.
Boundary A_z=0 at the outer rectangular box; use a wide air buffer in
industrial models. Heat deposited in conductor: q=0.5*sigma*|E|^2 [W/m3].
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import spsolve

MU0=4e-7*np.pi

def scalar_operator(nx:int, ny:int, dx:float, dy:float, sigma, freq_hz:float):
    """Return sparse complex matrix for A_z with Dirichlet A=0 walls."""
    n=nx*ny
    ex=np.ones(nx); ey=np.ones(ny)
    tx=diags([-ex[:-1],2*ex,-ex[:-1]],[-1,0,1],format='csr') /dx**2
    ty=diags([-ey[:-1],2*ey,-ey[:-1]],[-1,0,1],format='csr') /dy**2
    K=kron(eye(ny,format='csr'),tx,format='csr')+kron(ty,eye(nx,format='csr'),format='csr')
    K=K.astype(complex)+diags((1j*2*np.pi*freq_hz*MU0*np.asarray(sigma).reshape(-1)),format='csr')
    return K

def mesh(nx,ny,Lx,Ly):
    dx=Lx/(nx+1);dy=Ly/(ny+1)
    x=(np.arange(nx)+1)*dx;y=(np.arange(ny)+1)*dy
    return np.meshgrid(x,y),dx,dy

def manufactured(nx=24,ny=18,frequency_hz=25000.,sigma=8e4):
    """Manufactured smooth field; source calculated from continuous PDE."""
    Lx=.12;Ly=.05
    (xx,yy),dx,dy=mesh(nx,ny,Lx,Ly)
    exact=np.sin(np.pi*xx/Lx)*np.sin(np.pi*yy/Ly)*(1.+.25j)
    source=(np.pi**2/Lx**2+np.pi**2/Ly**2 + 1j*2*np.pi*frequency_hz*MU0*sigma)*exact
    op=scalar_operator(nx,ny,dx,dy,np.full((ny,nx),sigma),frequency_hz)
    numerical=spsolve(op,source.reshape(-1)).reshape(ny,nx)
    return {'nx':nx,'ny':ny,'relative_l2':float(np.linalg.norm(numerical-exact)/np.linalg.norm(exact)),
            'abs_residual':float(np.linalg.norm(op@numerical.reshape(-1)-source.reshape(-1))/np.linalg.norm(source))}

def run(nx=51,ny=25,Lx=.12,Ly=.02,freq_hz=30000.,coil_j_amp=7.5e7,
        bed_sigma=1.2e5, temperature=None,out='v5/results/em2d',write=True):
    (xx,yy),dx,dy=mesh(nx,ny,Lx,Ly)
    # Smooth coil legs and conductivity interfaces ensure a convergent mesh limit.
    # The Gaussian current regions approximate finite-width long straight conductors.
    bed_weight=(.5*(np.tanh((xx-.035)/.0015)+1.))*(.5*(np.tanh((.085-xx)/.0015)+1.))
    bed_weight*=.5*(np.tanh((yy-.003)/.001)+1.)*.5*(np.tanh((.017-yy)/.001)+1.)
    bed=bed_weight>.5
    if temperature is None:
        sigma=bed_sigma*bed_weight
    else:
        T=np.asarray(temperature)
        if T.shape!=(ny,nx):raise ValueError('temperature must be (ny,nx)')
        sigma=bed_sigma*bed_weight/np.maximum(.25,1+.003*(T-850.))
    sx=.005;sy=.0035
    lcoil=np.exp(-.5*((xx-.017)/sx)**2-.5*((yy-.01)/sy)**2)
    rcoil=np.exp(-.5*((xx-.103)/sx)**2-.5*((yy-.01)/sy)**2)
    j_source=coil_j_amp*(lcoil-rcoil)
    op=scalar_operator(nx,ny,dx,dy,sigma,freq_hz)
    az=spsolve(op,(MU0*j_source).reshape(-1)).reshape(ny,nx)
    ez=-1j*2*np.pi*freq_hz*az
    jed=ez*sigma
    q=.5*sigma*np.abs(ez)**2
    # Bx=dAz/dy, By=-dAz/dx. Time-averaged Lorentz in the solid.
    da_dy, da_dx=np.gradient(az,dy,dx)
    bx=da_dy;by=-da_dx
    fx=-.5*np.real(jed*np.conj(by));fy=.5*np.real(jed*np.conj(bx))
    residual=np.linalg.norm(op@az.reshape(-1)-(MU0*j_source).reshape(-1))/max(np.linalg.norm((MU0*j_source).reshape(-1)),1e-30)
    p_per_depth=float(q.sum()*dx*dy)
    summary={'model':'2D MQS A_z long-conductor Cartesian approximation',
      'geometry':'two oppositely directed parallel straight coil legs; NOT an axisymmetric induction coil',
      'nx':nx,'ny':ny,'frequency_hz':freq_hz,'coil_current_density_A_m2':coil_j_amp,
      'sigma_bed_S_m_ref':bed_sigma,'conducting_cells':int(bed.sum()),
      'total_joule_power_W_per_m_out_of_plane':p_per_depth,'mean_q_bed_W_m3':float(np.mean(q[bed])),
      'max_q_W_m3':float(q.max()),'max_B_T':float(np.max(np.hypot(np.abs(bx),np.abs(by)))),
      'algebraic_relative_residual':float(residual),'thermal_feedback':temperature is not None,
      'limitations':'illustrative sigma/coil geometry; no thermal/EM material measurements, no plasma, no geometric coil fidelity, no experimental EM validation'}
    if write:
        d=Path(out);d.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(d/'fields.npz',A_z_real=az.real,A_z_imag=az.imag,q=q,sigma=sigma,
                            Bx_abs=np.abs(bx),By_abs=np.abs(by),force_x=fx,force_y=fy,
                            x=xx[0],y=yy[:,0])
        (d/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary,{'q':q,'A':az,'x':xx[0],'y':yy[:,0],'Fx':fx,'Fy':fy,'Bx':bx,'By':by,'sigma':sigma}

def convergence(out='v5/results/em2d'):
    from matplotlib import pyplot as plt
    d=Path(out);d.mkdir(parents=True,exist_ok=True)
    manufactured_rows=[manufactured(n,int(.7*n)) for n in (18,36,72)]
    coils=[]
    for n in (31,51,81):
        s,_=run(nx=n,ny=max(19,int(n*.5)),write=False)
        coils.append({'nx':n,'power_W_per_m':s['total_joule_power_W_per_m_out_of_plane']})
    s,fields=run(out=d,write=True)
    (d/'validation.json').write_text(json.dumps({'manufactured':manufactured_rows,'coil_mesh':coils},indent=2)+'\n')
    fig,axs=plt.subplots(1,2,figsize=(11,3.8))
    im=axs[0].pcolormesh(fields['x']*1e3,fields['y']*1e3,fields['q'],shading='nearest',cmap='inferno')
    axs[0].set_title('Conductive-bed induction Joule heating')
    axs[0].set_xlabel('x [mm]');axs[0].set_ylabel('y [mm]');fig.colorbar(im,ax=axs[0],label='W/m³')
    axs[1].plot([a['nx'] for a in manufactured_rows],[a['relative_l2'] for a in manufactured_rows],'o-',label='Manufactured solution L2')
    axs[1].set_yscale('log');axs[1].set_xlabel('Grid cells in x');axs[1].set_ylabel('Relative L2 error')
    axs[1].grid(alpha=.25);axs[1].legend(fontsize=8)
    fig.tight_layout();fig.savefig(d/'em2d_verification.png',dpi=155);plt.close(fig)
    print(json.dumps({'em2d':s,'verification':manufactured_rows,'coil_mesh':coils},indent=2))
    return s,manufactured_rows,coils

if __name__=='__main__':
    convergence()
