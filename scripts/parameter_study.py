#!/usr/bin/env python3
"""Sensitivity to temperature-independent prescribed induction heating and gas flow."""
import argparse
import csv
from pathlib import Path
import subprocess

p=argparse.ArgumentParser()
p.add_argument('--exe',default='./build/reactor1d')
a=p.parse_args()
root=Path('results/parameter_study');root.mkdir(parents=True,exist_ok=True)
rows=[]
for U in (0.008,0.016,0.032):
    for heating in (300000,600000,900000):
        dest=root/f'U{U:g}_q{heating}'
        subprocess.run([a.exe,'--velocity',str(U),'--heat',str(heating),'--out',str(dest)],
                       check=True,stdout=subprocess.DEVNULL)
        data={}
        for line in (dest/'summary.txt').read_text().splitlines():
            if '=' in line:
                key,v=line.split('=',1)
                try:data[key]=float(v)
                except ValueError:pass
        rows.append({'superficial_velocity_m_s':U,
                     'heat_peak_W_m3':heating,
                     'mean_conversion':data['mean_conversion'],
                     'max_solid_temperature_K':data['maximum_solid_temperature_K'],
                     'H2_mass_balance_abs_mol_m2':data['max_absolute_H2_balance_error_mol_m2']})
with (root/'comparison.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('Wrote',root/'comparison.csv')
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(6.7,4.2))
    for U in (0.008,0.016,0.032):
        rr=[x for x in rows if x['superficial_velocity_m_s']==U]
        ax.plot([r['heat_peak_W_m3']/1000 for r in rr],
                [100*r['mean_conversion'] for r in rr],marker='o',label=f'U={U:g} m/s')
    ax.set(xlabel='Prescribed peak heat source (kW/m³)',
           ylabel='Mean NiO conversion (%)',
           title='Unvalidated reduced-order sensitivity study')
    ax.legend();fig.tight_layout();fig.savefig(root/'conversion_sensitivity.png',dpi=170)
    plt.close(fig)
except ImportError:
    pass
