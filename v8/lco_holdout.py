"""Transparent, very sparse LCO holdout experiment using *approximate prose* endpoints.
This is a SAME-PUBLICATION temperature holdout, not independent reactor validation.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from v8.physics import R
from v8.lco_kinetics import OBSERVED_T240_APPROX,EA_GB_J_MOL,EA_A3_J_MOL,TREF_K,ginstling_brounshtein,gb_inverse

def alpha(t_s,T_K,kref_gb,kref_a3,transition_s=1500):
    def arr(k,ea):return k*np.exp(-ea/R*(1/T_K-1/TREF_K))
    k1,k2=arr(kref_gb,EA_GB_J_MOL),arr(kref_a3,EA_A3_J_MOL)
    if t_s<=transition_s:return gb_inverse(min(k1*t_s,1/3-1e-12))
    x0=gb_inverse(min(k1*transition_s,1/3-1e-12))
    a0=(-np.log1p(-x0))**(1/3)
    return float(1-np.exp(-(a0+k2*(t_s-transition_s))**3))

def estimate(train_C=(800,1000),test_C=(700,1100)):
    if set(train_C)&set(test_C):raise ValueError('Train/test temperature leakage')
    def predict(tc,params):return alpha(240*60,tc+273.15,np.exp(params[0]),np.exp(params[1]))
    y=np.asarray([OBSERVED_T240_APPROX[x] for x in train_C])
    fit=least_squares(lambda par:np.array([predict(x,par) for x in train_C])-y,
        np.log([2e-5,2e-5]),bounds=(np.log([1e-10,1e-10]),np.log([1e-2,1e-2])),max_nfev=200)
    results={}
    for title,temps in [('train',train_C),('held_out',test_C)]:
        rows=[{'T_C':tc,'observed_alpha_240min_approx':OBSERVED_T240_APPROX[tc],
               'predicted_alpha_240min':predict(tc,fit.x)} for tc in temps]
        results[title]={'rows':rows,'MAE_conversion_fraction':float(np.mean([abs(q['predicted_alpha_240min']-q['observed_alpha_240min_approx']) for q in rows]))}
    return {'prefactors_kref_GB_and_A3_s_inverse_assumed_vs_calibrated':np.exp(fit.x).tolist(),
            'published_Ea_kJ_mol':[121,95],
            'transition_min_assumed':25,'fit_success':bool(fit.success),
            'source': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9821372/',
            'caveat':'only 4 approximate 240-min endpoints from one study; 2 fitted, 2 held out; not full experimental kinetic validation',
            **results}

def main(out='v8/results/lco_holdout.json'):
    res=estimate();p=Path(out);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(res,indent=2)+'\n');return res
if __name__=='__main__':print(json.dumps(main(),indent=2))
