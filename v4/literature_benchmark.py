#!/usr/bin/env python3
"""External published NiO/H2 data bench, NOT validation of the 2D reactor.

Rates [um/s] were *manually transcribed* from Hidayat et al. (2009),
Metallurgical and Materials Transactions B 40, Table III, 100 vol% hydrogen,
DOI 10.1007/s11663-008-9212-0. Uncertainties listed in Table III.
A 500C/675C two-point Arrhenius calibration predicts held-out 600C and the
high-temperature 800/900/1000C rates. 800C violates monotonic Arrhenius
trend; intentionally show cross-regime failure, not hide it.
"""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
T_C=np.array([500.,600.,675.,800.,900.,1000.])
RATE=np.array([0.61,2.08,4.04,2.80,3.82,5.25])
UNC=np.array([.09,.05,.12,.10,.10,.40])
R=8.314462618

def fit_and_score(out='v4/results/literature'):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    T_K=T_C+273.15
    calibration=np.array([0,2]);heldout=np.array([1,3,4,5])
    slope,intercept=np.polyfit(1./T_K[calibration],np.log(RATE[calibration]),1)
    Ea=-slope*R
    pred=np.exp(intercept+slope/T_K)
    errors=pred-RATE
    metrics={
      'title':'NiO initial Ni product layer growth under pure H2: external literature benchmark',
      'reference':'T. Hidayat et al. (2009), Table III, 100% H2, DOI 10.1007/s11663-008-9212-0',
      'data_kind':'measured dense synthetic NiO layer growth, um/s; NOT reactor conversion',
      'method':'Arrhenius log-linear fit at exactly 500 C and 675 C, holdout 600/800/900/1000 C',
      'calibration_temperatures_C':T_C[calibration].tolist(),
      'heldout_temperatures_C':T_C[heldout].tolist(),
      'derived_Ea_J_mol':float(Ea),'derived_A_um_s':float(np.exp(intercept)),
      'heldout_600C_absolute_error_um_s':float(abs(errors[1])),
      'heldout_high_T_MAE_um_s':float(np.mean(np.abs(errors[3:]))),
      'heldout_all_MAE_um_s':float(np.mean(np.abs(errors[heldout]))),
      'predictions_um_s':pred.tolist(),
      'observations_um_s':RATE.tolist(),
      'measurement_sd_um_s':UNC.tolist(),
      'caveat':'A fit to two points gives 0 fitting error by construction. The model cannot reproduce nonmonotonic high-T data and is NOT a calibrated 2D volumetric kinetic law.'}
    with (out/'data.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['T_C','measured_growth_um_per_s','std_um_per_s','arrhenius_prediction_um_per_s','set'])
        for i in range(len(T_C)):
            w.writerow([T_C[i],RATE[i],UNC[i],pred[i],'fit' if i in calibration else 'holdout'])
    (out/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(7,4.5))
    ax.errorbar(T_C,RATE,yerr=UNC,marker='o',linestyle='none',capsize=3,label='Hidayat et al. (2009), Table III')
    ax.plot(T_C,pred,'--',label='2-point Arrhenius model (extrapolation)')
    ax.scatter(T_C[calibration],RATE[calibration],s=85,facecolor='none',edgecolor='black',label='Calibration temperatures')
    ax.set_xlabel('Temperature (°C)');ax.set_ylabel('Initial Ni-layer growth (µm/s)')
    ax.set_title('Measured non-Arrhenius response: extrapolation fails');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'NiO_external_benchmark.png',dpi=170);plt.close(fig)
    print(json.dumps(metrics,indent=2))
    return metrics

if __name__=='__main__':fit_and_score()
