"""Grid, time-step, boundary-condition and model-switch sensitivity.

This is an uncertainty *screen*, NOT experimental validation or a calibrated
uncertainty quantification of battery recycling performance.
"""
import csv, json
from dataclasses import replace
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from v13.reactor import Config,run


def study(out_dir='v13/results/study'):
    out=Path(out_dir);out.mkdir(parents=True,exist_ok=True)
    base=Config(duration_s=2.,em_nr=16,em_nz=24)
    scenarios={
      'baseline':{},
      'grid12':{'nz':12},
      'grid24':{'nz':24},
      'dt_0125':{'dt_s':.0125},
      'dt_050':{'dt_s':.05},
      'EM_off':{'enable_em':False},
      'radicals_off':{'enable_radicals':False},
      'H2_1pct':{'inlet_H2_fraction':.01},
      'H2_6pct':{'inlet_H2_fraction':.06},
      'pressure_drop_1Pa':{'pressure_out_Pa':base.pressure_in_Pa-1.0},
    }
    rows=[]
    for name,kw in scenarios.items():
        conf=replace(base,**kw)
        res=run(conf)
        rows.append({'scenario':name,'nz':conf.nz,'dt_s':conf.dt_s,
          'mean_LCO_conversion':res['mean_LCO_conversion'],
          'mean_Tsolid_K':res['history'][-1]['Tsolid_mean_K'],
          'max_abs_element_error_mol_m2':res['max_element_balance_mol_m2'],
          'relative_energy_error':res['relative_energy_ledger_error'],
          'max_EOS_residual_Pa':res['max_end_step_EOS_residual_Pa'],
          'max_pressure_Picard_iterations':res['max_pressure_energy_Picard_iterations'],
          'max_EM_feedback_iterations':res['max_EM_feedback_iterations_per_refresh']})
        print(name,rows[-1]['mean_LCO_conversion'],res['relative_energy_ledger_error'],flush=True)
    with (out/'scenarios.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
    (out/'scenarios.json').write_text(json.dumps(rows,indent=2)+'\n')
    plt.figure(figsize=(11,5))
    plt.bar(range(len(rows)),[r['mean_LCO_conversion']*100 for r in rows])
    plt.xticks(range(len(rows)),[r['scenario'] for r in rows],rotation=35,ha='right')
    plt.ylabel('LiCoO2 conversion (%)')
    plt.title('V13 sensitivity, 2 seconds (uncalibrated reactor surrogate)')
    plt.tight_layout();plt.savefig(out/'conversion_sensitivity.png',dpi=160);plt.close()
    baseline=next(r['mean_LCO_conversion'] for r in rows if r['scenario']=='baseline')
    plt.figure(figsize=(11,5))
    plt.bar(range(len(rows)),[100*(r['mean_LCO_conversion']/baseline-1) for r in rows])
    plt.xticks(range(len(rows)),[r['scenario'] for r in rows],rotation=35,ha='right')
    plt.ylabel('Relative difference from baseline (%)')
    plt.title('Conversion sensitivity to selected assumptions (2-second illustrative case)')
    plt.tight_layout();plt.savefig(out/'relative_sensitivity.png',dpi=160);plt.close()
    plt.figure(figsize=(9,5))
    plt.plot([12,18,24],[next(r for r in rows if r['scenario']==s)['mean_LCO_conversion'] for s in ('grid12','baseline','grid24')],marker='o')
    plt.xlabel('Axial cells');plt.ylabel('LiCoO2 conversion (mole fraction)')
    plt.title('Spatial resolution screen (not GCI or full grid independence)')
    plt.tight_layout();plt.savefig(out/'grid_screen.png',dpi=160);plt.close()
    return rows

if __name__=='__main__':study()
