"""Generate explicitly synthetic parameter studies; not empirical validation."""
import csv,json
from dataclasses import replace
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from .reactor import Config,run

def main():
    dest=Path('v12/results/study');dest.mkdir(parents=True,exist_ok=True)
    base=Config(duration_s=2.)
    scenarios={
      'base':base,'grid12':replace(base,nz=12),'grid24':replace(base,nz=24),
      'dt0.05':replace(base,dt_s=.05),'dt0.0125':replace(base,dt_s=.0125),
      'no_em':replace(base,enable_em=False),'no_radicals':replace(base,enable_radicals=False),
      'H2_1pct':replace(base,inlet_H2_fraction=.01),
      'H2_6pct':replace(base,inlet_H2_fraction=.06),
      'low_carbon':replace(base,carbon_LCO_ratio=.8),
      'high_carbon':replace(base,carbon_LCO_ratio=2.5),
    }
    vals=[]
    for name,cfg in scenarios.items():
        r=run(cfg,dest/name)
        vals.append({'scenario':name,'mean_LCO_conversion':r['mean_LCO_conversion'],
          'H_radical_mean_mol_m3_pore':r['H_radical_mean_mol_m3_pore'],
          'max_element_balance_mol_m2':r['max_element_balance_mol_m2'],
          'relative_energy_ledger_error':r['relative_energy_ledger_error'],
          'max_pressure_closure_relative':r['max_implicit_pressure_closure_relative']})
        print(name, r['mean_LCO_conversion'],r['relative_energy_ledger_error'])
    with (dest/'scenarios.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(vals[0]));w.writeheader();w.writerows(vals)
    (dest/'summary.json').write_text(json.dumps({'scenarios':vals,
        'important':'Synthetic numerical sensitivity; no independent experimental validation.'},indent=2))
    fig,ax=plt.subplots(figsize=(7.6,4.0))
    ax.bar([v['scenario'] for v in vals],[100*v['mean_LCO_conversion'] for v in vals])
    ax.set_ylabel('Mean LiCoO2 conversion (%)')
    ax.tick_params(axis='x',rotation=45)
    ax.set_title('v12: synthetic parameter-sensitivity study (2 s)')
    ax.grid(axis='y',alpha=.2)
    fig.tight_layout();fig.savefig(dest/'conversion_sensitivity.png',dpi=170);plt.close(fig)
    z=(np.arange(base.nz)+.5)*base.length_m/base.nz
    fields=np.load(dest/'base'/'fields.npz')
    fig,ax=plt.subplots(figsize=(7.2,4))
    ax.plot(z,fields['pressure'],label='EOS pressure (Pa)')
    ax.set_xlabel('Axial position (m)');ax.set_ylabel('Pressure (Pa)')
    ax.set_title('v12: corrected gas-pressure profile')
    ax.legend();fig.tight_layout();fig.savefig(dest/'pressure_profile.png',dpi=170);plt.close(fig)
    return vals
if __name__=='__main__':main()
