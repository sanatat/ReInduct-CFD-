"""Paper-table-III initial Ni-layer growth rate holdouts. Not reactor calibration.

Uses exactly 19 means and standard deviations from Hidayat et al. (2009),
Table III. A temperature-blind holdout probes the key nonmonotonic behavior;
leave-one-composition-out tests transfer across H2 dilution. None are external-lab validation.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/nio_hidayat_2009_table_III.csv'
R=8.314462618

def read_data():
    with DATA.open() as f: rows=list(csv.DictReader(f))
    arr=np.array([[float(r[k]) for k in ['temperature_C','h2_mole_fraction','rate_um_per_s','sd_um_per_s']] for r in rows]);
    assert arr.shape==(19,4)
    return arr

def fit_global_arrhenius(rows):
    # k(T) y_H2, log-rate regression weighted by reported approx relative uncertainty
    T=rows[:,0]+273.15; x=1/T; y=np.log(rows[:,2]/rows[:,1]); weights=(rows[:,2]/rows[:,3])**2
    mat=np.column_stack((np.ones(len(rows)),x)); coef=np.linalg.lstsq(mat*np.sqrt(weights[:,None]),y*np.sqrt(weights),rcond=None)[0]
    return lambda T_C,h2: np.exp(coef[0]+coef[1]/(np.asarray(T_C)+273.15))*np.asarray(h2)

def fit_temperature_interpolator(rows):
    temps=np.unique(rows[:,0]); rates=[]
    for t in temps:
        s=rows[rows[:,0]==t]
        w=1/s[:,3]**2
        # Explicit first-order in p_H2: minimize weighted SSE for rate = k(T)*y
        k=float(np.sum(w*s[:,1]*s[:,2])/np.sum(w*s[:,1]**2))
        rates.append(k)
    interpolation=PchipInterpolator(temps,np.asarray(rates),extrapolate=False)
    return lambda T_C,h2: interpolation(T_C)*np.asarray(h2)

def metrics(obs,pred):
    z=np.asarray(obs)-np.asarray(pred)
    return {'n':int(len(z)),'MAE_um_s':float(np.mean(abs(z))),'RMSE_um_s':float(np.sqrt(np.mean(z*z))),
            'bias_um_s':float(np.mean(-z))}

def run(out=ROOT/'results'):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    table=read_data()
    # Prespecified holdout: T=800C; all concentrations excluded from fitting.
    train=table[table[:,0]!=800];test=table[table[:,0]==800]
    ar=fit_global_arrhenius(train);interp=fit_temperature_interpolator(train)
    predictions={'global_arrhenius':ar(test[:,0],test[:,1]),'PCHIP_temperature':interp(test[:,0],test[:,1])}
    holdout={name:metrics(test[:,2],vals) for name,vals in predictions.items()}
    # Secondary leave-one-pressure-fraction-out; temperatures included in train, hence easier task.
    leave_h2={}; target=.26
    tr=table[table[:,1]!=target];te=table[table[:,1]==target]
    for name,fit in [('global_arrhenius',fit_global_arrhenius),('PCHIP_temperature',fit_temperature_interpolator)]:
        m=fit(tr);leave_h2[name]=metrics(te[:,2],m(te[:,0],te[:,1]))
    # Leave one *interior* temperature completely unseen from fit.
    folds=[]
    for t in [600.,675.,800.,900.]:
        tr=table[table[:,0]!=t];te=table[table[:,0]==t]
        fold={'temperature_C':t,'n_test':int(len(te))}
        for name,fit in [('arrhenius',fit_global_arrhenius),('PCHIP',fit_temperature_interpolator)]:
            m=fit(tr);fold[name]=metrics(te[:,2],m(te[:,0],te[:,1]))
        folds.append(fold)
    report={'dataset':'Hidayat 2009 Table III; dense NiO initial product-layer advancement [um/s]',
            'n_total':len(table),'holdout_800C':holdout,'leave_H2_26pct_out':leave_h2,
            'leave_one_temperature_out_interior':folds,'scope':'internal paper holdouts; independent reactor validation NOT performed',
            'units':'rates [um/s]; do NOT interpret as mol/m3/s without interfacial-area and solid-density assumptions'}
    (out/'kinetics_holdout.json').write_text(json.dumps(report,indent=2)+'\n')
    with (out/'holdout_800C.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['temperature_C','y_H2','observed_um_s','sd_um_s','global_arrhenius_pred_um_s','PCHIP_pred_um_s'])
        for row,a,b in zip(test,predictions['global_arrhenius'],predictions['PCHIP_temperature']):
            w.writerow([row[0],row[1],row[2],row[3],a,b])
    # Main plot independently held-out 800 C measurements.
    fig,ax=plt.subplots(figsize=(8,5))
    ax.errorbar(test[:,1],test[:,2],yerr=test[:,3],fmt='o',label='Published observation (holdout)')
    for name,v in predictions.items():ax.plot(test[:,1],v,marker='s',label=name)
    ax.set(xlabel='H2 mole fraction [-]',ylabel='Initial Ni layer growth [um/s]',title='800°C held out entirely from fitting')
    ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(out/'holdout_800C.png',dpi=160);plt.close(fig)
    return report

if __name__=='__main__':
    print(json.dumps(run(),indent=2))
