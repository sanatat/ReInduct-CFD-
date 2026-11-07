"""Two-way EM <-> heat fixed-point proof-of-concept with cylindrical coil.

Uses v6 frequency-domain axisymmetric EM solver. On every fixed-point iteration:
Ts -> sigma(Ts) -> A_theta -> Joule q -> Ts; radiation and gas exchange included.
This is a thermal/EM benchmark, NOT full reactive/flow OpenFOAM co-simulation.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from v6.axisymmetric_em import Grid,run as em_run
from v8.physics import volumetric_radiation_loss

def solve(nr=18,nz=30, max_iter=35,relaxation=.55,
          dt_s=80.,ambient_K=700., initial_K=850.,gas_K=825.,
          wall_emissivity=.7,heat_transfer_W_m3K=600.,
          coil_current_A_m2=2.0e6, conductivity_W_mK=1.0, tol_K=2e-5, write=False,out='v8/results/em_thermal'):
    grid=Grid(nr=nr,nz=nz)
    z,r=grid.coords()
    bed=(r>=.005)&(r<=.025)&(z>=.040)&(z<=.100)
    rmax=.025
    Ts=np.full((nz,nr),initial_K)
    wall_area_per_volume=np.zeros_like(Ts)
    radial_edges=np.arange(nr+1)*grid.dr
    # The outermost bed ring represents exposed cylindrical reactor boundary.
    ring=np.flatnonzero(grid.r<rmax)
    if not len(ring):raise ValueError('EM grid lacks reactor cells')
    wall_i=int(ring[-1])
    wall_area_per_volume[:,wall_i]=2*np.pi*radial_edges[wall_i+1]*grid.dz/(np.pi*((wall_i+1)**2-wall_i**2)*grid.dr**2*grid.dz)
    # Expose the last ring only where the bed exists; radiation is NOT counted over internal particles.
    wall_area_per_volume*=bed
    # Assemble conservative cylindrical heat conduction on the actual bed region.
    # Rows are divided by cylindrical cell volume, giving W/(m3 K).
    ids=-np.ones_like(Ts,dtype=int)
    pairs=np.argwhere(bed)
    ids[bed]=np.arange(len(pairs))
    cp=1.05e6
    thermal_matrix=lil_matrix((len(pairs),len(pairs)),dtype=float)
    for idx,(j,i) in enumerate(pairs):
        volume=np.pi*(radial_edges[i+1]**2-radial_edges[i]**2)*grid.dz
        diagonal=cp/dt_s+heat_transfer_W_m3K
        for j2,i2 in [(j-1,i),(j+1,i),(j,i-1),(j,i+1)]:
            if not (0<=j2<nz and 0<=i2<nr and bed[j2,i2]):continue
            if j2!=j:
                area=np.pi*(radial_edges[i+1]**2-radial_edges[i]**2)
                distance=grid.dz
            else:
                face=radial_edges[max(i,i2)]
                area=2*np.pi*face*grid.dz
                distance=grid.dr
            coeff=conductivity_W_mK*area/(distance*volume)
            thermal_matrix[idx,ids[j2,i2]]-=coeff
            diagonal+=coeff
        thermal_matrix[idx,idx]+=diagonal
    thermal_matrix=thermal_matrix.tocsr()
    history=[]
    q=np.zeros_like(Ts)
    for k in range(max_iter):
        em_summary,em=em_run(grid,bed_temperature_K=Ts,coil_current_A_m2=coil_current_A_m2)
        q=em['q']*bed
        rad=wall_emissivity*5.670374419e-8*wall_area_per_volume*(Ts**4-ambient_K**4)
        # Conservative implicit FV conduction + implicit heat exchange, explicit
        # lagged radiation (iterated to convergence), thermal inertia in pseudo-step.
        rhs=(cp*initial_K/dt_s+q+heat_transfer_W_m3K*gas_K-rad)[bed]
        target=np.full_like(Ts,initial_K)
        target[bed]=spsolve(thermal_matrix,rhs)
        if not np.all(np.isfinite(target)) or np.min(target[bed])<1:
            raise RuntimeError('nonphysical temperature in EM thermal iteration')
        new=Ts+relaxation*(target-Ts)
        delta=float(np.max(np.abs(new[bed]-Ts[bed])))
        history.append({'iteration':k+1,'max_temperature_change_K':delta,
                        'peak_bed_temperature_K':float(new[bed].max()),
                        'bed_joule_power_W':float(np.sum(q*2*np.pi*r*grid.dr*grid.dz)),
                        'radiation_loss_W':float(np.sum(rad*2*np.pi*r*grid.dr*grid.dz))})
        Ts=new
        if delta<tol_K:break
    else:
        raise RuntimeError(f'EM-thermal fixed point did not converge: {delta} K')
    _, em_final=em_run(grid,bed_temperature_K=Ts,coil_current_A_m2=coil_current_A_m2)
    q_final=em_final['q']*bed
    _,cold_em=em_run(grid,bed_temperature_K=np.full_like(Ts,initial_K),coil_current_A_m2=coil_current_A_m2)
    v=2*np.pi*r*grid.dr*grid.dz
    p_final=float(np.sum(q_final*v))
    p_cold=float(np.sum(cold_em['q']*bed*v))
    report={'model':'axisymmetric coil MQS + grey-wall radiation + gas/solid exchange, two-way EM-thermal fixed point',
            'grid_nr_nz':[nr,nz],'iterations':len(history),
            'converged':True,'last_delta_K':history[-1]['max_temperature_change_K'],
            'peak_temperature_K':float(Ts[bed].max()),
            'initial_temperature_K':initial_K,
            'power_constant_sigma_W':p_cold,'power_feedback_W':p_final,
            'power_feedback_relative_change_pct':100*(p_final/p_cold-1) if p_cold else 0.,
            'wall_emissivity':wall_emissivity,'solid_conductivity_W_mK':conductivity_W_mK,
            'limitations':['not a reactor CFD simulation','thermal pseudo-steady approximation; conduction included',
                           'illustrative conductivity/coil/heat-transfer properties','no experimental coil validation']}
    if write:
        out=Path(out);out.mkdir(parents=True,exist_ok=True)
        (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
        (out/'iterations.json').write_text(json.dumps(history,indent=2)+'\n')
        np.savez_compressed(out/'fields.npz',r=grid.r,z=grid.z,Ts=Ts,qEM=q_final)
    return report,history

if __name__=='__main__':
    print(json.dumps(solve(write=True)[0],indent=2))
