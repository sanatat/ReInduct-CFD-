"""V11 independent coupled axisymmetric-induction / axial finite-volume battery model.

NOT an OpenFOAM/industrial reactor. A transparent reduced-order verification target.

Three element-balanced *illustrative pathways*, extents R [mol/(m3 bulk s)]:
 A: 2LiCoO2 + C(s) -> Li2O + 2CoO + CO(g)
 B: 2CoO + 2C(s) -> 2Co + 2CO(g)
 C: Li2O + CO2(g) -> Li2CO3.
Published Nuraeni stage-specific isothermal model constants are used as effective
first-order proxies only, NOT validated local heterogeneous reaction rates.

FV conventions: all solid inventories mol per bulk m3, gas concentration mol per
pore-gas m3, convective phi=u_s*c; axial cell width dz, area 1m2 for ledgers.
The gas mixture is weakly compressible and uses a common superficial velocity
computed from Ergun at the instantaneous mean gas state; this is NOT a
Navier-Stokes solver and no claim of local gas-pressure closure is made.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict,replace
from pathlib import Path
import argparse,json
import numpy as np
from v6.axisymmetric_em import Grid,run as solve_em
from v10.lco_published_rates import EARLY,LATE,arrhenius_fit

R_GAS=8.31446261815324
SIGMA_SB=5.670374419e-8
CP_MOLAR=30.0  # placeholder J/mol/K, all gases equally for exact sensible budget
SOLID_CP_BULK=1.05e6
H_A=48.0
# Reaction enthalpies: purely illustrative parameter sensitivity values, not calorimetry.
DELTA_H=np.array([45000.,-30000.,-25000.])
SPECIES=('LiCoO2','C_s','Li2O','CoO','Co','Li2CO3','CO','CO2','N2')
# Each column is one balanced reaction, with signed mole contributions.
STOICH=np.array([
 [-2.,0.,0.],   # LiCoO2
 [-1.,-2.,0.],  # solid carbon
 [1.,0.,-1.],   # Li2O
 [2.,-2.,0.],   # CoO
 [0.,2.,0.],    # Co metal
 [0.,0.,1.],    # Li2CO3
 [1.,2.,0.],    # gas CO
 [0.,0.,-1.],   # gas CO2
 [0.,0.,0.],    # inert N2
])
# Li, Co, C, O, N atoms per molecular formula, rows correspond to SPECIES.
ATOMS=np.array([
 [1,1,0,2,0], [0,0,1,0,0], [2,0,0,1,0],
 [0,1,0,1,0], [0,1,0,0,0], [2,0,1,3,0],
 [0,0,1,1,0], [0,0,1,2,0], [0,0,0,0,2]
],dtype=float)
assert np.max(np.abs(ATOMS.T@STOICH))<1e-14

@dataclass(frozen=True)
class Config:
    nz:int=30
    length_m:float=.06
    dt_s:float=.10
    duration_s:float=240.
    porosity:float=.45
    particle_diameter_m:float=.003
    pressure_gradient_Pa_m:float=8.
    inlet_T_K:float=1080.
    initial_T_K:float=1080.
    ambient_T_K:float=1010.
    initial_LCO_mol_m3:float=2000.
    carbon_LCO_ratio:float=1.8
    inlet_CO2_mole_fraction:float=.12
    inlet_pressure_Pa:float=101325.
    radiation_emissivity:float=.70
    radiation_exposure_factor:float=.15
    conduction_W_m_K:float=.7
    energy_reaction_scale:float=1.0
    coil_current_A_m2:float=2.5e6
    em_updates:int=5
    em_grid_nr:int=22
    em_grid_nz:int=36
    carbonate_k_s:float=8.e-5 # UNCALIBRATED; proxy only
    enable_chemistry:bool=True
    enable_em:bool=True
    enable_radiation:bool=True

    def validated(self):
        if self.nz<8 or self.dt_s<=0 or self.duration_s<=0 or self.length_m<=0:raise ValueError('invalid grid/time')
        if not 0<self.porosity<1 or self.carbon_LCO_ratio<0:raise ValueError('invalid porosity/carbon')
        if not 0<=self.inlet_CO2_mole_fraction<1:raise ValueError('invalid CO2 fraction')
        if self.inlet_T_K<973.15 or self.initial_T_K<973.15:raise ValueError('initial/inlet T outside published 700C+ range')
        if self.em_updates<1:raise ValueError('em_updates must >=1')
        return self


def kinetic_rates_per_s(T):
    """Arrhenius interpolant from 2023 published stage-fitted k/min, 700–1100C only."""
    temp=np.asarray(T,dtype=float)
    if np.any((temp<973.15)|(temp>1373.15)):
        raise ValueError('Temperature outside 700-1100C: no unmarked extrapolation')
    rates=[]
    for rows in (EARLY,LATE):
        fit=arrhenius_fit(rows)
        Ea=fit['E_a_J_per_mol']
        arr=np.array(rows,float)
        b=float(np.mean(np.log(arr[:,2])+Ea/(R_GAS*(arr[:,0]+273.15))))
        rates.append(np.exp(b-Ea/(R_GAS*temp))/60.)
    return np.array(rates)


def ergun_superficial(cfg, temperature_K, gas_conc):
    # constant pressure-gradient closure, not Navier-Stokes; gas state changes drag
    t=max(200.,float(np.mean(temperature_K)))
    y=np.asarray(gas_conc).mean(axis=1)
    rho=max(0.01,float(np.dot(y,np.array([.0280101,.0440095,.0280134]))))
    viscosity=3.5e-5*(t/1080.)**.7
    eps=cfg.porosity;d=cfg.particle_diameter_m
    a=150*viscosity*(1-eps)**2/(eps**3*d*d)
    b=1.75*rho*(1-eps)/(eps**3*d)
    return 2*cfg.pressure_gradient_Pa_m/(a+np.sqrt(a*a+4*b*cfg.pressure_gradient_Pa_m))


def _upwind_rhs(c,u,cin,dz,porosity):
    flux=np.r_[u*cin,u*c]
    return -(flux[1:]-flux[:-1])/(dz*porosity),flux[0],flux[-1]


def _solid_conduction(T,ks,dz):
    flux=np.zeros(T.size+1)
    flux[1:-1]=-ks*(T[1:]-T[:-1])/dz
    return -(flux[1:]-flux[:-1])/dz


def inventory_by_elements(solid,gas,porosity,dz):
    arr=np.concatenate((solid,gas*porosity),axis=0)
    return (ATOMS.T@arr).sum(axis=1)*dz


def _temperature_field_to_em_grid(Ts,cfg,grid):
    bed_z=grid.z
    bed_T=np.interp(bed_z-.04,(np.arange(cfg.nz)+.5)*cfg.length_m/cfg.nz,Ts,
                    left=Ts[0],right=Ts[-1])
    return np.broadcast_to(bed_T[:,None],(grid.nz,grid.nr)).copy()


def calculate_em_heat(Ts,cfg):
    """Axisymmetric EM with material temperature feedback, then conservative annular averaging."""
    z=(np.arange(cfg.nz)+.5)*cfg.length_m/cfg.nz
    grid=Grid(nr=cfg.em_grid_nr,nz=cfg.em_grid_nz)
    temperatures=_temperature_field_to_em_grid(Ts,cfg,grid)
    info,fields=solve_em(grid=grid,coil_current_A_m2=cfg.coil_current_A_m2,
                         bed_temperature_K=temperatures)
    dr=grid.dr; dz=grid.dz
    r=grid.r;zr=grid.z
    in_r=(r>=.005)&(r<=.025)
    in_z=(zr>=.04)&(zr<=.10)
    ring=np.pi*((r+dr/2)**2-(r-dr/2)**2)
    weights=ring[in_r]
    q_z=(fields['q'][:,in_r]*weights[None,:]).sum(axis=1)/weights.sum()
    # cell-center interpolation; rescale to match EM annular power after mapping.
    in_domain_z=zr[in_z]-.04
    in_domain_q=q_z[in_z]
    q_cfd=np.interp(z,in_domain_z,in_domain_q,left=in_domain_q[0],right=in_domain_q[-1])
    bed_A=np.pi*(.025**2-.005**2)
    p_em=float(np.sum((fields['q'][np.ix_(in_z,in_r)])*weights[None,:])*dz)
    p_mapped=float(np.sum(q_cfd)*(cfg.length_m/cfg.nz)*bed_A)
    q_cfd*=p_em/p_mapped
    mapped=float(np.sum(q_cfd)*(cfg.length_m/cfg.nz)*bed_A)
    return q_cfd, {'annulus_EM_power_W':p_em,'mapped_CFD_power_W':mapped,
                   'relative_power_error':abs(mapped-p_em)/max(p_em,1e-16),
                   'em_algebraic_residual':info['algebraic_relative_residual'],
                   'max_q_W_m3':float(q_cfd.max())}


def simulate(cfg=Config(),output_dir=None):
    cfg.validated()
    nz=cfg.nz;eps=cfg.porosity;dz=cfg.length_m/nz
    nstep=round(cfg.duration_s/cfg.dt_s)
    if abs(nstep*cfg.dt_s-cfg.duration_s)>1e-8:raise ValueError('duration/dt nonintegral')
    Ts=np.full(nz,cfg.initial_T_K);Tg=Ts.copy()
    solid=np.zeros((6,nz));solid[0]=cfg.initial_LCO_mol_m3
    solid[1]=cfg.initial_LCO_mol_m3*cfg.carbon_LCO_ratio
    cin0=cfg.inlet_pressure_Pa/(R_GAS*cfg.inlet_T_K)
    cin=np.array([0.,cfg.inlet_CO2_mole_fraction*cin0,
                 (1-cfg.inlet_CO2_mole_fraction)*cin0])
    gas=np.broadcast_to(cin[:,None],(3,nz)).copy()
    Eg=eps*CP_MOLAR*gas.sum(axis=0)*Tg
    Es=SOLID_CP_BULK*Ts
    init_elements=inventory_by_elements(solid,gas,eps,dz)
    init_energy=float((Eg+Es).sum()*dz)
    external_elements=np.zeros(5);external_energy=0.0
    em_updates=[]; snapshots=[]
    qEM=np.zeros(nz);em_schedule=set(np.linspace(0,nstep-1,min(nstep,cfg.em_updates),dtype=int))
    A_wall_over_V=2*.025/(.025**2-.005**2)*cfg.radiation_exposure_factor
    max_cfl=0.
    for step in range(nstep):
        if cfg.enable_em and (step in em_schedule):
            qEM,em_report=calculate_em_heat(Ts,cfg)
            em_updates.append({'time_s':step*cfg.dt_s,'Ts_mean_K':float(Ts.mean()),**em_report})
        # temperature and species dependent fixed pressure-gradient Ergun closure
        u=ergun_superficial(cfg,Tg,gas)
        cfl=u*cfg.dt_s/(dz*eps)
        max_cfl=max(max_cfl,cfl)
        if cfl>.65:raise ValueError(f'Gas convection unstable CFL={cfl:.3f}')
        r1=np.zeros(nz);r2=np.zeros(nz);r3=np.zeros(nz)
        if cfg.enable_chemistry:
            k1,k2=kinetic_rates_per_s(Ts)
            # Extents are in mol reaction/(bulk m3 s)
            r1=np.minimum(.5*k1*np.maximum(solid[0],0.),np.minimum(np.maximum(solid[0],0.)/2,np.maximum(solid[1],0.))/cfg.dt_s)
            # sequential reaction uses current (old) CoO; no artificially created reactant this step
            r2=np.minimum(.5*k2*np.maximum(solid[3],0.),np.minimum(np.maximum(solid[3],0.)/2,np.maximum(solid[1],0.)/2)/cfg.dt_s)
            carbon_rate=r1+2*r2
            carbon_limit=np.maximum(solid[1],0.)/cfg.dt_s
            scale=np.minimum(1.0,carbon_limit/np.maximum(carbon_rate,1e-40))
            r1*=scale; r2*=scale
            r3=np.minimum(cfg.carbonate_k_s*solid[2]*np.minimum(1.,gas[1]/max(cin0,1e-30)),
                          np.minimum(solid[2],eps*gas[1])/cfg.dt_s)
        reaction=np.vstack((r1,r2,r3))
        delta_s=STOICH[:6]@reaction
        delta_g=(STOICH[6:9]@reaction)/eps
        # Gas fluxes per area across inlet/outlet; all fluxes use OLD concentrations
        gas_rhs=np.zeros_like(gas)
        flux_in=np.zeros(3);flux_out=np.zeros(3)
        for k in range(3):
            gas_rhs[k],flux_in[k],flux_out[k]=_upwind_rhs(gas[k],u,cin[k],dz,eps)
        # Element ledger includes net gas transport (solid is fixed)
        external_elements+=cfg.dt_s*(ATOMS[6:9].T@(flux_in-flux_out))
        # Gas energy advection and heat exchange; reaction creates/destroys
        # pore-gas sensible enthalpy at solid temperature, opposite solid term.
        gas_source_bulk=eps*delta_g.sum(axis=0)
        gas_sensible_transfer=CP_MOLAR*gas_source_bulk*Ts
        h_transfer=H_A*(Ts-Tg)
        _,hin,hout=_upwind_rhs(Eg/eps,u,cin.sum()*CP_MOLAR*cfg.inlet_T_K,dz,eps)
        # Flux divergence for bulk energy (do NOT divide by porosity)
        energy_adv=-(np.r_[u*cin.sum()*CP_MOLAR*cfg.inlet_T_K,u*(Eg/eps)][1:]
                     -np.r_[u*cin.sum()*CP_MOLAR*cfg.inlet_T_K,u*(Eg/eps)][:-1])/dz
        qrad=(cfg.radiation_emissivity*SIGMA_SB*A_wall_over_V*
              (Ts**4-cfg.ambient_T_K**4)) if cfg.enable_radiation else np.zeros(nz)
        qs=cfg.energy_reaction_scale*(-DELTA_H@reaction)
        # All energy source/sinks accounted in ledger below
        external_energy+=cfg.dt_s*(hin-hout + float(np.sum(qEM-qrad+qs)*dz))
        Eg+=cfg.dt_s*(energy_adv+h_transfer+gas_sensible_transfer)
        Es+=cfg.dt_s*(_solid_conduction(Ts,cfg.conduction_W_m_K,dz)-h_transfer
                       -gas_sensible_transfer+qEM-qrad+qs)
        solid+=cfg.dt_s*delta_s
        gas+=cfg.dt_s*(gas_rhs+delta_g)
        # no hiding negative inventories behind clipping
        if min(float(solid.min()),float(gas.min()))<-1e-10:raise FloatingPointError(f'negative inventory at {step}: solid {solid.min()} species {np.unravel_index(solid.argmin(),solid.shape)}, gas {gas.min()}')
        if float(gas.sum(axis=0).min())<1e-6:raise FloatingPointError('gas depleted')
        Tg=Eg/(eps*CP_MOLAR*gas.sum(axis=0))
        Ts=Es/SOLID_CP_BULK
        if min(Ts.min(),Tg.min())<200:raise FloatingPointError('unphysical T')
        if (step+1)%max(1,nstep//10)==0 or step==nstep-1:
            snapshots.append({'time_s':(step+1)*cfg.dt_s,'mean_Tsolid_K':float(Ts.mean()),
                'max_Tsolid_K':float(Ts.max()),'mean_Tgas_K':float(Tg.mean()),
                'mean_LCO_conversion':float(1-solid[0].mean()/cfg.initial_LCO_mol_m3),
                'mean_CO_mol_m3_pore':float(gas[0].mean()),
                'mean_Li2CO3_mol_m3_bulk':float(solid[5].mean()),'superficial_velocity_m_s':u})
    inferred_pressure=gas.sum(axis=0)*R_GAS*Tg
    pressure_max_relative_deviation=float(np.max(np.abs(inferred_pressure/cfg.inlet_pressure_Pa-1.)))
    final_elements=inventory_by_elements(solid,gas,eps,dz)
    final_energy=float((Eg+Es).sum()*dz)
    elt_error=final_elements-init_elements-external_elements
    energy_error=final_energy-init_energy-external_energy
    conversion=float(1-solid[0].mean()/cfg.initial_LCO_mol_m3)
    report={'version':'11.0','solver':'conservative 1D axial FV + 2D axisymmetric MQS feedback (not OpenFOAM)',
      'notes':'Stage k/min Nuraeni 2023 used as uncalibrated effective first-order proxies; stoichiometric pathway illustrative.',
      'settings':asdict(cfg),'steps':nstep,'max_gas_CFL':max_cfl,
      'element_names':['Li','Co','C','O','N'],
      'initial_atom_inventories_mol_per_m2':init_elements.tolist(),
      'final_atom_inventories_mol_per_m2':final_elements.tolist(),
      'boundary_net_atom_inventories_mol_per_m2':external_elements.tolist(),
      'element_conservation_residual_mol_per_m2':elt_error.tolist(),
      'max_abs_element_residual':float(np.max(np.abs(elt_error))),
      'initial_total_energy_J_m2':init_energy,'final_total_energy_J_m2':final_energy,
      'external_energy_J_m2':external_energy,'energy_ledger_residual_J_m2':energy_error,
      'relative_energy_residual':abs(energy_error)/max(abs(init_energy),abs(external_energy),1.),
      'mean_conversion_LCO':conversion,'mean_final_Ts_K':float(Ts.mean()),
      'max_final_Ts_K':float(Ts.max()),'outlet_CO_mol_m3_pore':float(gas[0,-1]),
      'mean_lithium_carbonate_mol_m3_bulk':float(solid[5].mean()),
      'max_inferred_ideal_gas_pressure_deviation_fraction':pressure_max_relative_deviation,
      'continuity_limit':'No compressible pressure correction; pressure deviation is monitored, not corrected',
      'em_updates':em_updates,'snapshots':snapshots,
      'experimental_reactor_validated':False,'openfoam_runtime_verified':False,
      'multiphase_Navier_Stokes':False,'full_plasma_MHD':False}
    if output_dir:
        out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
        (out/'evidence.json').write_text(json.dumps(report,indent=2)+'\n')
        np.savez_compressed(out/'fields.npz', z_m=(np.arange(nz)+.5)*dz,Ts_K=Ts,Tg_K=Tg,
                            solid_inventories=solid,gas_concentrations=gas,qEM_W_m3=qEM)
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='v11/results/base');
    p.add_argument('--duration',type=float,default=240.);p.add_argument('--nz',type=int,default=30)
    a=p.parse_args();result=simulate(Config(duration_s=a.duration,nz=a.nz),a.output)
    print(json.dumps({k:result[k] for k in ('mean_conversion_LCO','mean_final_Ts_K','max_abs_element_residual','relative_energy_residual','em_updates')},indent=2))
