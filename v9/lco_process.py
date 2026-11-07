"""Kinetics-informed, carbon-constrained LiCoO2/C batch scenarios.

This is NOT independent validation. The time-dependent conversion is supplied
by v8's stage-integral model; its prefactors were estimated using 4 approximate
endpoints of ONE publication, not raw measured time histories.

Elemental accounting uses low-T Li2CO3 / Co / CO2 route below 900 C and high-T
Li2O / Co / CO route above 1000 C. The high-T gas product is an assumption.
"""
from __future__ import annotations
import json
from pathlib import Path
from v7.validation.battery_phase_envelope import LOW,HIGH
from v4.battery_chemistry import formula_mass,atomic_residual
from v8.lco_holdout import alpha
from v8.lco_kinetics import KREF_GB_PER_S,KREF_A3_PER_S

def scenario(temp_c=800.,time_min=120.,lco_mol=10.,carbon_mol=7.5,
             kref_gb=KREF_GB_PER_S,kref_a3=KREF_A3_PER_S):
    if lco_mol<=0 or carbon_mol<0 or time_min<0:raise ValueError('Invalid inputs')
    if 700<=temp_c<=900:stoich=LOW;regime='low_T_carbonate_CO2'
    elif 1000<=temp_c<=1100:stoich=HIGH;regime='high_T_oxide_CO_assumed'
    else:raise ValueError('No phase data for 900–1000 C gap or outside tested range')
    kinetic_conversion=alpha(time_min*60,temp_c+273.15,kref_gb,kref_a3)
    max_lco_reacted=min(lco_mol,carbon_mol*(-stoich['LiCoO2'])/(-stoich['C']))
    # Carbon shortage caps conversion, rather than creating negative carbon.
    lco_reacted=min(lco_mol*kinetic_conversion,max_lco_reacted)
    ext=lco_reacted/(-stoich['LiCoO2'])
    inventories={sp:nu*ext for sp,nu in stoich.items()}
    inventories['LiCoO2']+=lco_mol;inventories['C']+=carbon_mol
    if min(inventories.values())<-1e-10:raise RuntimeError('Negative inventory')
    original=lco_mol*formula_mass('LiCoO2')+carbon_mol*formula_mass('C')
    final=sum(formula_mass(sp)*n for sp,n in inventories.items())
    return {'temperature_C':temp_c,'time_min':time_min,'phase_regime':regime,
            'unconstrained_kinetic_conversion':kinetic_conversion,
            'carbon_limited_conversion_LCO':lco_reacted/lco_mol,
            'carbon_available_mol':carbon_mol,'LCO_initial_mol':lco_mol,
            'inventories_mol':inventories,'elemental_balance':atomic_residual(stoich),
            'input_mass_g':original,'final_mass_g':final,'mass_residual_g':final-original,
            'interpretation':'Kinetic scenario with assumed prefactors and gas speciation; not measured recovery'}

def sweep(out='v9/results/lco_scenarios.json'):
    rows=[]
    for t in (700,800,900,1000,1100):
        for tm in (0,60,120,240):
            for carbon in (2.5,7.5):
                rows.append(scenario(t,tm,10,carbon))
    p=Path(out);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(rows,indent=2)+'\n')
    return rows

if __name__=='__main__':print(json.dumps({'scenarios':len(sweep())},indent=2))
