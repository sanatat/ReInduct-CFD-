#!/usr/bin/env python3
"""Small factorial design-of-experiments on ILLUSTRATIVE model parameters.

Screen conversion/peak temperature/total induced heat, not a calibrated
industrial reactor optimization. Operating ranges and temperature constraint
are illustrative, chosen to demonstrate transparent process screening.
"""
import csv,json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from reactor2d.reactor import simulate

def main():
    out=ROOT/'results/v3/DOE';out.mkdir(parents=True,exist_ok=True)
    table=[]
    for U in [.008,.016,.024]:
        for h0 in [4000.,9000.,14000.]:
            s=simulate(nx=30,ny=10,dt=.003,end=8.,Umean=U,h0=h0,
                       out=out/f'U{U:g}_H{int(h0)}')
            table.append({'U_m_s':U,'H0_A_m':h0,'conversion':s['mean_conversion'],
                          'max_T_K':s['max_solid_T_K'],
                          'Joule_heat_J_per_m':s['cumulative_EM_heat_J_per_m']})
    with (out/'screening.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
    # Demonstration constraint, not a real equipment safety limit.
    temp_threshold=920.
    feasible=[r for r in table if r['max_T_K']<=temp_threshold]
    # Pareto: maximize conversion, minimize energy input (assuming positive values).
    pareto=[]
    for v in feasible:
        dominated=any((w['conversion']>=v['conversion'] and
                       w['Joule_heat_J_per_m']<=v['Joule_heat_J_per_m'] and
                      (w['conversion']>v['conversion'] or
                       w['Joule_heat_J_per_m']<v['Joule_heat_J_per_m'])) for w in feasible)
        if not dominated:pareto.append(v)
    fig,ax=plt.subplots(figsize=(7.2,4.6))
    for U in sorted({v['U_m_s'] for v in table}):
        subset=[v for v in table if v['U_m_s']==U]
        ax.plot([v['Joule_heat_J_per_m']/1000 for v in subset],
                [100*v['conversion'] for v in subset], 'o-',label=f'U = {U:.3f} m/s')
    ax.set(xlabel='Accumulated induction heat [kJ/m depth]',
           ylabel='Simulated conversion [%]',
           title='Illustrative process design screening — not experimental')
    ax.legend();ax.grid(alpha=.25);fig.tight_layout()
    fig.savefig(out/'conversion_vs_heat.png',dpi=160);plt.close(fig)
    kinetics=[]
    for Ea in [75000.,85000.,95000.]:
        s=simulate(nx=30,ny=10,dt=.003,end=8.,kinetic_Ea=Ea,out=out/f'kinetics_Ea{int(Ea)}')
        kinetics.append({'Ea_J_mol':Ea,'conversion':s['mean_conversion']})
    with (out/'kinetics_sensitivity.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['Ea_J_mol','conversion']);w.writeheader();w.writerows(kinetics)
    (out/'summary.json').write_text(json.dumps({
        'assumed_max_T_K_constraint':temp_threshold,
        'number_of_cases':len(table),'feasible_cases':len(feasible),
        'non_dominated_feasible':pareto,'Ea_screening':kinetics,
        'caution':'engineering screening with unvalidated surrogate and illustrative property values'
    },indent=2)+'\n')
    print('Complete: DOE',len(table),'cases; feasible',len(feasible),
          '; Pareto feasible',len(pareto),'\n',out)
if __name__=='__main__':main()
