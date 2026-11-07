"""V12 conservative low-dimensional thermo-reactive column / MQS induction surrogate.

Gas: CO, CO2, N2, H2, H. Solids: LCO, C, Li2O, CoO, Co, Li2CO3.
The radical kinetics is an illustrative power-limited H2 <-> 2H model.
No ions, plasma field, electron kinetic equation, or direct hydrogen reduction
of LCO is claimed. This is not a full plasma/MHD/OpenFOAM implementation.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict
from pathlib import Path
import json
import numpy as np
from v11.integrated_reactor import (Config as LegacyConfig, STOICH, kinetic_rates_per_s,
                                     calculate_em_heat, _solid_conduction)
from v13.compressible_transport import R,pressure_corrected_step
from v13.thermo_pressure import advance as coupled_pressure_energy_advance
from v13.em_feedback import solve_em_thermal_feedback
from v13.gas_properties import gas_energy,enthalpy_molar
CP=30.0 # J/mol/K representative; verified sensible-energy balance, NOT measured mixture data
SOLID_CP=1.05e6 # J/m3/K representative
H_TRANSFER=48. # W/m3/K
BOND_H2=436000. # J/mol H2 dissociated, illustrative bond enthalpy
SB=5.670374419e-8
# Li, Co, C, O, N, H atom counts, corresponding to 6 solid + 5 gas species
ATOMS=np.array([
[1,1,0,2,0,0], [0,0,1,0,0,0], [2,0,0,1,0,0],
[0,1,0,1,0,0], [0,1,0,0,0,0],[2,0,1,3,0,0],
[0,0,1,1,0,0],[0,0,1,2,0,0],[0,0,0,0,2,0],
[0,0,0,0,0,2],[0,0,0,0,0,1]],dtype=float)
S=np.zeros((11,4));S[:9,:3]=STOICH;S[9,3]=-1.;S[10,3]=2.
assert np.max(np.abs(ATOMS.T@S))<1e-14
@dataclass(frozen=True)
class Config:
    nz:int=18
    length_m:float=.06
    dt_s:float=.025
    duration_s:float=5.0
    porosity:float=.45
    pressure_in_Pa:float=101330.
    pressure_out_Pa:float=101329.5
    T_in_K:float=1080.
    T_initial_K:float=1080.
    T_ambient_K:float=1010.
    initial_LCO_mol_m3:float=2000.
    carbon_LCO_ratio:float=1.8
    inlet_CO2_fraction:float=.12
    inlet_H2_fraction:float=.03
    # Illustrative power-limited radical chemistry, explicitly not a fitted plasma law
    radical_dissociation_rate_s:float=0.04
    radical_recombination_m3_mol_s:float=0.001
    plasma_electron_power_W_m3:float=4000.
    # Uncalibrated surface chemistry switch; heat of solid reactions is intentionally DISABLED
    carbonate_rate_s:float=8e-5
    radiation_emissivity:float=.7
    wall_exposure_factor:float=.15
    coil_current_A_m2:float=2.5e6
    em_nr:int=16
    em_nz:int=24
    em_update_interval_s:float=1.
    em_feedback_tolerance:float=1e-7
    enable_em:bool=True
    enable_solid_reactions:bool=True
    enable_radicals:bool=True
    enable_radiation:bool=True
    def check(self):
        if self.nz<5 or self.dt_s<=0 or self.length_m<=0 or self.duration_s<=0: raise ValueError('grid/time')
        if self.pressure_in_Pa<=self.pressure_out_Pa or self.pressure_out_Pa<1e3:raise ValueError('pressure BC')
        if not 0<self.porosity<1 or not 0<=self.inlet_H2_fraction<=1 or not 0<=self.inlet_CO2_fraction<=1 or self.inlet_H2_fraction+self.inlet_CO2_fraction>1:raise ValueError('mixture')
        if self.em_update_interval_s<=0 or not 0<self.em_feedback_tolerance<1: raise ValueError('EM coupling interval/tolerance')
        if self.coil_current_A_m2<0 or self.radiation_emissivity<0 or self.radiation_emissivity>1:raise ValueError('heat source parameters')
        if not 973.15<=self.T_initial_K<=1373.15:raise ValueError('outside kinetic domain')
        if self.plasma_electron_power_W_m3<0 or self.radical_dissociation_rate_s<0 or self.radical_recombination_m3_mol_s<0:raise ValueError('negative plasma/radical parameters')
        return self

def run(cfg=Config(),output=None):
    cfg.check();dz=cfg.length_m/cfg.nz;n=cfg.nz;eps=cfg.porosity
    steps=round(cfg.duration_s/cfg.dt_s)
    if abs(steps*cfg.dt_s-cfg.duration_s)>1e-9:raise ValueError('time step must exactly divide duration')
    species_names=['CO','CO2','N2','H2','H'];element_names=['Li','Co','C','O','N','H']
    y_in=np.array([0.,cfg.inlet_CO2_fraction,1.-cfg.inlet_CO2_fraction-cfg.inlet_H2_fraction,cfg.inlet_H2_fraction,0.])
    Tgas=np.full(n,cfg.T_initial_K);Tsolid=Tgas.copy()
    inlet_c=cfg.pressure_in_Pa/(R*cfg.T_in_K)
    gas=np.repeat((y_in*inlet_c)[:,None],n,axis=1)
    solid=np.zeros((6,n));solid[0]=cfg.initial_LCO_mol_m3;solid[1]=cfg.carbon_LCO_ratio*cfg.initial_LCO_mol_m3
    Eg=gas_energy(gas,Tgas,eps);Es=SOLID_CP*Tsolid
    chemE=.5*BOND_H2*eps*gas[4] # H radical potential energy relative to molecular hydrogen
    atom_initial=ATOMS.T@np.r_[solid,gas*eps]
    atom_initial=atom_initial.sum(axis=1)*dz
    energy_initial=(Eg+Es+chemE).sum()*dz
    boundary_atoms=np.zeros(6);external_energy=0.;max_closure=0.;min_p=np.inf;max_p=0.
    qEM=np.zeros(n);em_reports=[];history=[];max_cfl=0.;pressure_profiles=[];max_thermo_pressure_iter=0;max_eos_residual_Pa=0.;max_em_feedback_iters=0;max_em_feedback_err=0.
    update_steps=max(1,round(cfg.em_update_interval_s/cfg.dt_s))
    for step in range(steps):
        em_refresh=cfg.enable_em and step % update_steps==0
        extents=np.zeros((4,n))
        if cfg.enable_solid_reactions:
            k1,k2=kinetic_rates_per_s(Tsolid)
            x1=np.minimum(0.5*k1*np.maximum(solid[0],0),np.minimum(solid[0]/2,solid[1])/cfg.dt_s)
            x2=np.minimum(0.5*k2*np.maximum(solid[3],0),np.minimum(solid[3]/2,solid[1]/2)/cfg.dt_s)
            scale=np.minimum(1.,np.maximum(solid[1],0.)/np.maximum((x1+2*x2)*cfg.dt_s,1e-30))
            extents[0]=x1*scale;extents[1]=x2*scale
            extents[2]=np.minimum(cfg.carbonate_rate_s*solid[2]*np.minimum(1.,gas[1]/max(inlet_c,1e-30)),
                  np.minimum(solid[2],eps*gas[1])/cfg.dt_s)
        if cfg.enable_radicals:
            # mol bulk/m3/s for H2 -> 2H, limited by electron energy available
            d=eps*cfg.radical_dissociation_rate_s*np.maximum(gas[3],0.)
            r=eps*cfg.radical_recombination_m3_mol_s*np.maximum(gas[4],0.)**2/2.
            d=np.minimum(d,eps*gas[3]/cfg.dt_s)
            d=np.minimum(d,(cfg.plasma_electron_power_W_m3/BOND_H2)+r)
            r=np.minimum(r,eps*gas[4]/(2*cfg.dt_s))
            extents[3]=d-r
        delta=S@extents
        species_source=delta[6:]
        # Conservative pressure, species and gas-sensible-energy correction.
        # Iterate thermal EOS rather than stopping at the lagged-T pressure solution.
        solid_gas_moles=species_source[:3].sum(axis=0)
        # Conservatively transfer sensible energy carried by gas species
        # generated/consumed by heterogeneous solid pathways (CO/CO2).
        released_sensible=np.sum(species_source[:3]*enthalpy_molar(Tsolid)[:3],axis=0)
        gas_solid=H_TRANSFER*(Tsolid-Tgas)
        area_over_vol=2*.025/(.025**2-.005**2)*cfg.wall_exposure_factor
        qrad=cfg.radiation_emissivity*SB*area_over_vol*(Tsolid**4-cfg.T_ambient_K**4) if cfg.enable_radiation else np.zeros(n)
        qplasma=cfg.plasma_electron_power_W_m3 if cfg.enable_radicals else 0.
        dissociation_energy=BOND_H2*extents[3]
        energy_other=gas_solid+released_sensible+qplasma-dissociation_energy
        (newgas,flux,sflux,pinfo,Eg_next,Tgas_next,flux_energy,pressure_info)=coupled_pressure_energy_advance(
             gas,Tgas,Eg,species_source,energy_other,dt=cfg.dt_s,length=cfg.length_m,
             porosity=eps,pin=cfg.pressure_in_Pa,pout=cfg.pressure_out_Pa,
             inlet_fraction=y_in,T_in=cfg.T_in_K)
        gas_dE=(Eg_next-Eg)/cfg.dt_s
        max_thermo_pressure_iter=max(max_thermo_pressure_iter,pressure_info['iterations'])
        max_eos_residual_Pa=max(max_eos_residual_Pa,pressure_info['pressure_EOS_residual_Pa'])
        max_closure=max(max_closure,pinfo['pressure_closure_relative'])
        max_cfl=max(max_cfl,pinfo['max_species_CFL'])
        min_p=min(min_p,float(np.min(pinfo['pressure_Pa'])));max_p=max(max_p,float(np.max(pinfo['pressure_Pa'])))
        if pinfo['max_species_CFL']>0.7:raise FloatingPointError('advective CFL too large')
        boundary_atoms+=cfg.dt_s*ATOMS[6:].T@(sflux[:,0]-sflux[:,-1])
        solid_rhs_without_em=_solid_conduction(Tsolid,0.7,dz)-gas_solid-released_sensible-qrad
        if em_refresh:
            lc=LegacyConfig(nz=n,length_m=cfg.length_m,em_grid_nr=cfg.em_nr,
                            em_grid_nz=cfg.em_nz,coil_current_A_m2=cfg.coil_current_A_m2)
            qEM,er=solve_em_thermal_feedback(
                Es,solid_rhs_without_em,cfg.dt_s,SOLID_CP,Tsolid,
                lambda next_T: calculate_em_heat(next_T,lc),
                tolerance=cfg.em_feedback_tolerance)
            em_reports.append({'time_s':step*cfg.dt_s,**er})
            max_em_feedback_iters=max(max_em_feedback_iters,er['feedback_iterations'])
            max_em_feedback_err=max(max_em_feedback_err,er['feedback_relative_power_density_error'])
        sol_dE=solid_rhs_without_em+qEM
        external_energy+=cfg.dt_s*(float(flux_energy[0]-flux_energy[-1])
             +.5*BOND_H2*float(sflux[4,0]-sflux[4,-1])
             +float(np.sum((qEM-qrad+qplasma)*dz)))
        gas=newgas
        solid+=cfg.dt_s*delta[:6]
        Eg+=cfg.dt_s*gas_dE;Es+=cfg.dt_s*sol_dE
        Tgas=Tgas_next;Tsolid=Es/SOLID_CP
        if min(float(gas.min()),float(solid.min()))<-1e-9:raise FloatingPointError('negative inventory')
        if np.min(Tgas)<200 or np.min(Tsolid)<200:raise FloatingPointError('unphysical temperature')
        if np.any(Tsolid<973.15) or np.any(Tsolid>1373.15):raise FloatingPointError('reaction temperature outside data range')
        if step%max(1,steps//10)==0 or step==steps-1:
            pressure_profiles.append(gas.sum(axis=0)*R*Tgas)
            history.append({'time_s':(step+1)*cfg.dt_s,'Tsolid_mean_K':float(Tsolid.mean()),
               'LCO_conversion':float(1.-solid[0].mean()/cfg.initial_LCO_mol_m3),
               'H_radical_mean_mol_m3_pore':float(gas[4].mean()),
               'p_mean_Pa':float((gas.sum(axis=0)*R*Tgas).mean())})
    atom_final=(ATOMS.T@np.r_[solid,gas*eps]).sum(axis=1)*dz
    chemE=.5*BOND_H2*eps*gas[4]
    energy_final=(Eg+Es+chemE).sum()*dz
    atom_err=atom_final-atom_initial-boundary_atoms
    energy_err=energy_final-energy_initial-external_energy
    output_report={'version':'13.0','type':'research reduced-order finite-volume surrogate (NOT OpenFOAM CFD)',
         'time_s':cfg.duration_s,'steps':steps,'config':asdict(cfg),
         'mean_LCO_conversion':float(1-solid[0].mean()/cfg.initial_LCO_mol_m3),
         'H_radical_mean_mol_m3_pore':float(gas[4].mean()),
         'max_element_balance_mol_m2':float(np.max(np.abs(atom_err))),
         'element_residual_by_species_mol_m2':dict(zip(element_names,atom_err.tolist())),
         'relative_energy_ledger_error':float(abs(energy_err)/max(1.,abs(energy_initial),abs(external_energy))),
         'max_implicit_pressure_closure_relative':max_closure,
         'max_advective_CFL':max_cfl,
         'max_pressure_energy_Picard_iterations':max_thermo_pressure_iter,
         'max_end_step_EOS_residual_Pa':max_eos_residual_Pa,
         'max_EM_feedback_iterations_per_refresh':max_em_feedback_iters,
         'max_EM_feedback_relative_residual':max_em_feedback_err,
         'min_pressure_Pa':min_p,'max_pressure_Pa':max_p,
         'end_pressure_profile_Pa':(gas.sum(axis=0)*R*Tgas).tolist(),
         'electromagnetic_updates':em_reports,
         'history':history,
         'limitations':['1D porous low-Mach / Ergun, not Navier-Stokes','pressure and gas sensible enthalpy Picard corrected using NIST Shomate Cp; reaction formation enthalpies still omitted',
           'Solid reaction formation enthalpies deliberately omitted: not thermochemically calibrated',
           'Radical dissociation/recombination illustrative: no electron Boltzmann or sheath physics',
           'EM heat feedback converged only on refresh steps, gas/solid sources lagged; no full monolithic MHD coupling',
           'No OpenFOAM runtime verification','No independent reactor experimental validation'],
         'verified_independent_reactor_data':False,'verified_openfoam':False}
    if output:
        p=Path(output);p.mkdir(parents=True,exist_ok=True)
        (p/'evidence.json').write_text(json.dumps(output_report,indent=2)+'\n')
        np.savez_compressed(p/'fields.npz',solid=solid,gas=gas,Tgas=Tgas,Tsolid=Tsolid,qEM=qEM,
                pressure=gas.sum(axis=0)*R*Tgas,z_m=(np.arange(n)+.5)*dz)
    return output_report
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='v13/results/base');ap.add_argument('--duration',type=float,default=5.)
    args=ap.parse_args()
    print(json.dumps(run(Config(duration_s=args.duration),args.output),indent=2))
