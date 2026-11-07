#!/usr/bin/env python3
"""Run three benchmark resolutions, and quantify differences (no exact-solution claim)."""
from pathlib import Path
import subprocess, argparse, csv
p=argparse.ArgumentParser()
p.add_argument('--exe', default='./build/reactor1d')
a=p.parse_args()
root=Path('results/grid_study');root.mkdir(parents=True,exist_ok=True)
rows=[]
for n in (30,60,120):
    target=root/f'n{n}'
    subprocess.run([a.exe,'--nx',str(n),'--dt','0.005','--out',str(target)],check=True)
    data={}
    for line in (target/'summary.txt').read_text().splitlines():
        if '=' in line:
            k,v=line.split('=',1)
            try:data[k]=float(v)
            except ValueError:pass
    rows.append({'cells':n,'conversion':data['mean_conversion'],
                 'Tmax_K':data['maximum_solid_temperature_K'],
                 'H2_balance_error_mol_m2':data['max_absolute_H2_balance_error_mol_m2']})
with (root/'comparison.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('Grid study complete. Differences relative to finest grid:')
for r in rows:
    print(f"N={r['cells']:3d} Xbar={r['conversion']:.6f} Tmax={r['Tmax_K']:.3f}K "
          f"DeltaX_from_N120={abs(r['conversion']-rows[-1]['conversion']):.3g}")
