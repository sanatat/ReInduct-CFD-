#!/usr/bin/env python3
"""Literature-reported LiCoO2 + C overall carbon-reduction stoichiometry.

The model provides atomic/mass conservation, ideal stoichiometric phase yields,
not thermodynamic equilibrium or measured recovery performance. See
https://doi.org/10.1039/D3CS00898C for the reported overall reaction:
4 LiCoO2 + 3 C -> 2 Li2CO3 + 4 Co + CO2.
"""
from __future__ import annotations
import json
import re
from collections import defaultdict
from pathlib import Path
ELEMENT_MASS={'Li':6.94,'Co':58.933194,'C':12.011,'O':15.999}
REACTION={'LiCoO2':-4,'C':-3,'Li2CO3':2,'Co':4,'CO2':1}

def formula_atoms(formula):
    atoms=defaultdict(int)
    for symbol,n in re.findall(r'([A-Z][a-z]?)(\d*)',formula):atoms[symbol]+=int(n or 1)
    if ''.join(x+y for x,y in re.findall(r'([A-Z][a-z]?)(\d*)',formula))!=formula:
        raise ValueError(f'Unsupported chemical formula {formula}')
    return dict(atoms)

def atomic_residual(reaction=REACTION):
    residual=defaultdict(int)
    for compound,coef in reaction.items():
        for element,num in formula_atoms(compound).items():residual[element]+=coef*num
    return dict(residual)

def formula_mass(compound):return sum(ELEMENT_MASS[e]*n for e,n in formula_atoms(compound).items())

def batch(feed_mol=10.,carbon_mol=None,conversion=0.6):
    if carbon_mol is None:carbon_mol=.75*feed_mol
    if feed_mol<=0 or carbon_mol<0 or not 0<=conversion<=1:raise ValueError('invalid batch inputs')
    extent=min(feed_mol/4.,carbon_mol/3.)*conversion
    amounts={'LiCoO2':feed_mol-4*extent,'C':carbon_mol-3*extent,
             'Li2CO3':2*extent,'Co':4*extent,'CO2':extent}
    start_mass=feed_mol*formula_mass('LiCoO2')+carbon_mol*formula_mass('C')
    end_mass=sum(n*formula_mass(s) for s,n in amounts.items())
    return {'reaction':'4LiCoO2 + 3C -> 2Li2CO3 + 4Co + CO2',
            'reported_in':'DOI 10.1039/D3CS00898C (reported overall pathway, not single-step kinetics)',
            'feed_LiCoO2_mol':feed_mol,'feed_C_mol':carbon_mol,
            'conversion_of_stoichiometrically_accessible_reactants':conversion,
            'molar_amounts':amounts,'total_input_mass_g':start_mass,'total_output_mass_g':end_mass,
            'absolute_mass_balance_error_g':abs(start_mass-end_mass),
            'atomic_balance':atomic_residual(),
            'co_metal_yield_g':amounts['Co']*formula_mass('Co'),
            'lithium_carbonate_yield_g':amounts['Li2CO3']*formula_mass('Li2CO3'),
            'measured_recovery_percent':None,
            'disclaimer':'illustrative stoichiometric calculation, not calibrated process recovery'}

if __name__=='__main__':
    path=Path('v4/results/lco_stoichiometry.json');path.parent.mkdir(exist_ok=True,parents=True)
    r=batch();path.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
