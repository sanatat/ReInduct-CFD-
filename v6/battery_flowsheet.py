"""LiCoO2/carbon recycling: strictly stoichiometric scenario MATERIAL BALANCE.

4 LiCoO2 + 3 C -> 2 Li2CO3 + 4 Co + CO2 (overall proposed pathway).
No equilibrium model, measured process recovery, or kinetic prediction.
Assumed extent is provided by the user/scenario, not fitted to experiments.
"""
from pathlib import Path
import csv
import json
import numpy as np
from v4.battery_chemistry import formula_mass,atomic_residual

STOICH={'LiCoO2':-4.,'C':-3.,'Li2CO3':2.,'Co':4.,'CO2':1.}

def simulate(feed_kg=1.0,lco_mass_fraction=.45,carbon_mass_fraction=.12,conversion=.6,
             extra_carbon_g=0.,recovery_separation_fraction=1.0):
    """Track reactor and separation inventories separately to avoid false recovery claims.
    recovery_separation_fraction applies to phase separation, is an assumption.
    """
    if not (feed_kg>0 and 0<=lco_mass_fraction<=1 and 0<=carbon_mass_fraction<=1
        and lco_mass_fraction+carbon_mass_fraction<=1 and 0<=conversion<=1
        and extra_carbon_g>=0 and 0<=recovery_separation_fraction<=1):
        raise ValueError('invalid scenario fractions/amounts')
    feed_g=feed_kg*1000
    moles_lco=feed_g*lco_mass_fraction/formula_mass('LiCoO2')
    moles_carbon=(feed_g*carbon_mass_fraction+extra_carbon_g)/formula_mass('C')
    inert_g=feed_g*(1-lco_mass_fraction-carbon_mass_fraction)
    max_ext=min(moles_lco/4,moles_carbon/3)
    ext=conversion*max_ext
    amounts={'LiCoO2':moles_lco-4*ext,'C':moles_carbon-3*ext,'Li2CO3':2*ext,
             'Co':4*ext,'CO2':ext}
    solids=['LiCoO2','C','Li2CO3','Co']
    final_total=sum(amounts[p]*formula_mass(p) for p in amounts)+inert_g
    initial_total=feed_g+extra_carbon_g
    reco=recovery_separation_fraction
    returned={'feed_mass_g':feed_g,'extra_C_g':extra_carbon_g,
       'assumed_reaction_conversion_of_limiting_stoichiometric_extent':conversion,
       'assumed_separation_recovery_fraction':reco,
       'theoretical_Co_formed_g':amounts['Co']*formula_mass('Co'),
       'theoretical_Li2CO3_formed_g':amounts['Li2CO3']*formula_mass('Li2CO3'),
       'assumed_collected_Co_g':amounts['Co']*formula_mass('Co')*reco,
       'assumed_collected_Li2CO3_g':amounts['Li2CO3']*formula_mass('Li2CO3')*reco,
       'CO2_generated_g':amounts['CO2']*formula_mass('CO2'),
       'unreacted_LiCoO2_g':amounts['LiCoO2']*formula_mass('LiCoO2'),
       'carbon_remaining_g':amounts['C']*formula_mass('C'),
       'inert_remaining_g':inert_g,
       'phase_amounts_mol':amounts,
       'mass_balance_residual_g':final_total-initial_total,
       'elemental_balance_reaction':atomic_residual(),
       'process_validation':'none - theoretical stoichiometric and assumed separation/conversion only',
       'source':'overall reaction from review DOI:10.1039/D3CS00898C'}
    return returned

def study(out='v6/results/battery'):
    d=Path(out);d.mkdir(parents=True,exist_ok=True)
    rows=[]
    for carbon in [0.,25.,75.]:
        for c in [.2,.4,.6,.8,1.]:
            r=simulate(extra_carbon_g=carbon,conversion=c,recovery_separation_fraction=.8)
            rows.append({'extra_C_g':carbon,'assumed_conversion':c,
              'assumed_Co_collected_g':r['assumed_collected_Co_g'],
              'assumed_Li2CO3_collected_g':r['assumed_collected_Li2CO3_g'],
              'CO2_generated_g':r['CO2_generated_g'],'mass_balance_error_g':r['mass_balance_residual_g']})
    with (d/'hypothetical_scenarios.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    baseline=simulate()
    (d/'baseline_theoretical.json').write_text(json.dumps(baseline,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(6.5,4))
    for mass in [0,25,75]:
        sub=[r for r in rows if r['extra_C_g']==mass]
        ax.plot([r['assumed_conversion'] for r in sub],
                [r['assumed_Co_collected_g'] for r in sub],marker='o',label=f'Added C {mass} g')
    ax.set_title('Theoretical cobalt collection, assumed 80% separation')
    ax.set_xlabel('Assumed fraction of stoichiometrically possible conversion')
    ax.set_ylabel('Assumed collected Co [g] per 1 kg illustrative feed')
    ax.legend();ax.grid(alpha=.3);fig.tight_layout()
    fig.savefig(d/'assumed_recovery_scenarios.png',dpi=170);plt.close(fig)
    return baseline,rows

if __name__=='__main__':
    print(json.dumps(study()[0],indent=2))
