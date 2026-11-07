"""Energy-accounted 0D H2 <-> 2H plasma-chemistry *surrogate*.

Electron impact is represented by an effective first-order H2 dissociation rate,
NOT by electron Boltzmann kinetics, ion transport, Poisson's equation, or a
resolved plasma discharge. The applied electron power is an external input.
Molar inventories are per m3. E_e and E_g are in J/m3, dissociation energy is
J/mol H2. The conserved inventory is (2*H2+H) mol H atoms/m3; total energy
includes radical chemical potential (0.5*H*D) and external power.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import json
from pathlib import Path
import math

DISSOCIATION_J_PER_MOL_H2=436000.0

@dataclass(frozen=True)
class PlasmaCase:
    h2_mol_m3:float=1.0
    h_mol_m3:float=0.0
    electron_energy_J_m3:float=30000.0
    neutral_sensible_energy_J_m3:float=0.0
    electron_heating_W_m3:float=3000.0
    k_diss_s:float=0.035
    k_recomb_m3_mol_s:float=0.012
    dt_s:float=0.005
    duration_s:float=10.0


def run(cfg:PlasmaCase=PlasmaCase(), save_path=None):
    pars=asdict(cfg)
    if not all(math.isfinite(x) for x in pars.values()):raise ValueError('nonfinite input')
    if any(x<0 for x in pars.values()) or cfg.dt_s<=0:raise ValueError('nonphysical input')
    N=round(cfg.duration_s/cfg.dt_s)
    if N<1 or abs(N*cfg.dt_s-cfg.duration_s)>1e-9:raise ValueError('duration/dt not integral')
    h2=cfg.h2_mol_m3; h=cfg.h_mol_m3
    ee=cfg.electron_energy_J_m3; eg=cfg.neutral_sensible_energy_J_m3
    h_atoms0=2*h2+h
    E0=ee+eg+0.5*h*DISSOCIATION_J_PER_MOL_H2
    history=[]; peak_atom_residual=0.; peak_energy_residual=0.
    for n in range(N):
        # External energy deposited at the start of the step. Dissociation energy
        # is debited from the electron reservoir; recombination heats neutrals.
        ee+=cfg.electron_heating_W_m3*cfg.dt_s
        rd=min(cfg.k_diss_s*h2,h2/cfg.dt_s,ee/(DISSOCIATION_J_PER_MOL_H2*cfg.dt_s))
        rr=min(cfg.k_recomb_m3_mol_s*h*h,h/(2*cfg.dt_s))
        delta=(rd-rr)*cfg.dt_s
        h2-=delta;h+=2*delta
        ee-=rd*DISSOCIATION_J_PER_MOL_H2*cfg.dt_s
        eg+=rr*DISSOCIATION_J_PER_MOL_H2*cfg.dt_s
        ea=abs(2*h2+h-h_atoms0)
        e_total=ee+eg+0.5*h*DISSOCIATION_J_PER_MOL_H2
        eb=abs(e_total-E0-cfg.electron_heating_W_m3*cfg.dt_s*(n+1))
        peak_atom_residual=max(peak_atom_residual,ea)
        peak_energy_residual=max(peak_energy_residual,eb)
        if min(h2,h,ee,eg)<-1e-7:raise FloatingPointError('negative inventory or energy')
        if n==0 or (n+1)%max(1,N//50)==0 or n==N-1:
            history.append({'time_s':(n+1)*cfg.dt_s,'H2_mol_m3':h2,'H_mol_m3':h,
                            'electron_energy_J_m3':ee,'neutral_energy_J_m3':eg,
                            'H_atoms_residual_mol_m3':ea,'energy_residual_J_m3':eb})
    report={'model':'Zero-dimensional effective H2/H dissociation/recombination with electron-energy reservoir',
            'initial_condition':pars,'final_H2_mol_m3':h2,'final_H_mol_m3':h,
            'atomic_H_balance_abs_max_mol_m3':peak_atom_residual,
            'total_energy_balance_abs_max_J_m3':peak_energy_residual,
            'final_electron_energy_J_m3':ee,'final_neutral_energy_J_m3':eg,
            'claim_limit':'Reduced-rate kinetics with imposed power; NOT a plasma discharge, MHD, or validated black-mass chemistry'}
    if save_path is not None:
        p=Path(save_path);p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(json.dumps({'report':report,'history':history},indent=2)+'\n')
    return report,history

if __name__=='__main__':
    r,_=run(save_path='v9/results/plasma_energy.json');print(json.dumps(r,indent=2))
