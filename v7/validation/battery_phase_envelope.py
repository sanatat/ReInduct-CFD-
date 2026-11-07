"""Literature-compatible *qualitative* LiCoO2/C phase envelope.

Source: carbothermic LiCoO2 study (2023) reports Li2CO3+Co at 700-900C
and Li2O+Co at 1000-1100C. Balancing a representative high-T route
with CO rather than CO2 is a MODEL ASSUMPTION, NOT published gas speciation.
The scenarios assume conversion and are NOT kinetic predictions or validations.
"""
from __future__ import annotations
import json
from pathlib import Path
from v4.battery_chemistry import atomic_residual,formula_mass
LOW={'LiCoO2':-4,'C':-3,'Li2CO3':2,'Co':4,'CO2':1}
HIGH={'LiCoO2':-2,'C':-3,'Li2O':1,'Co':2,'CO':3}

def scenario(T_C,LCO_mol=10,carbon_mol=15,conversion_fraction=.60):
    if not 0<=conversion_fraction<=1 or LCO_mol<0 or carbon_mol<0:
        raise ValueError('invalid nonnegative inputs')
    if 700<=T_C<=900:reaction=LOW;observed='Li2CO3 + Co';gas_assumption='CO2 in representative balance'
    elif 1000<=T_C<=1100:reaction=HIGH;observed='Li2O + Co';gas_assumption='CO in representative balance (assumed gas stoichiometry)'
    else:raise ValueError('Source does not establish product phase within 900-1000C gap or outside 700-1100C')
    assert all(x==0 for x in atomic_residual(reaction).values())
    ext=conversion_fraction*min(LCO_mol/-reaction['LiCoO2'],carbon_mol/-reaction['C'])
    amount={sp:coef*ext for sp,coef in reaction.items()}
    amount['LiCoO2']+=LCO_mol;amount['C']+=carbon_mol
    mi=LCO_mol*formula_mass('LiCoO2')+carbon_mol*formula_mass('C')
    mf=sum(formula_mass(sp)*n for sp,n in amount.items())
    return {'T_C':T_C,'observed_solid_phase_class':observed,'gas_product_assumption':gas_assumption,
            'assumed_conversion_fraction_of_stoichiometric_limit':conversion_fraction,
            'remaining_and_products_mol':amount,'mass_balance_error_g':mf-mi,
            'scope':'Literature-informed phase selection and theoretical mass balance; no kinetic or process-yield validation'}

def main(out='v7/results/battery_phase'):
    out=Path(out);out.mkdir(exist_ok=True,parents=True)
    results=[scenario(t) for t in (700,800,900,1000,1100)]
    (out/'phase_scenarios.json').write_text(json.dumps(results,indent=2)+'\n')
    return results
if __name__=='__main__':print(json.dumps(main(),indent=2))
