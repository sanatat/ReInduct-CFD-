"""Group-heldout chemical-kinetic comparisons without model selection leakage.

Published Hidayat 2009 Table III source is embedded in v5/kinetics_validation.py.
Measurement target is early Ni product-layer growth [micrometre/s], NOT bed rate.
For each of four H2 fractions train on all other fractions; report every fold.
All three models declared before looking at folds; no post-hoc 'best' claim.
"""
import json
from pathlib import Path
import numpy as np
from v5.kinetics_validation import ROWS,fit,predict

MODELS=[('single_arrhenius','arrhenius',None),('two_regime_free_order','segmented',None),
        ('two_regime_first_order','segmented',1.)]

def crossvalidate(out='v6/results/kinetics',write=True):
    a=np.asarray(ROWS,float);T,p,y,sd=a.T
    results={}
    for label,mode,order in MODELS:
        folds=[];all_abs=[]
        for held in sorted(np.unique(p)):
            train=~np.isclose(p,held); test=~train
            coef=fit(T[train],y[train],sd[train],p[train],mode,order)
            yh=predict(coef,T[test],p[test],mode,order)
            err=(yh-y[test]);all_abs.extend(np.abs(err))
            folds.append({'heldout_H2_fraction':float(held),'n_train':int(np.sum(train)),
                          'n_test':int(np.sum(test)),'MAE_um_s':float(np.mean(np.abs(err))),
                          'measured_um_s':y[test].tolist(),'pred_um_s':yh.tolist(),
                          'temperatures_C':T[test].tolist()})
        results[label]={'folds':folds,'pooled_MAE_um_s':float(np.mean(all_abs))}
    results['protocol']='Leave-one-H2-fraction-out at 0.26,0.53,0.75,1.00; candidate models fixed a priori for THIS assessment'
    results['provenance']='T. Hidayat et al. (2009), Metall Mater Trans B, Table III, DOI 10.1007/s11663-008-9212-0'
    results['measurement']='Initial Ni layer growth in dense synthetic NiO, micrometre per second; not porous black mass conversion'
    results['interpretation']='Within-one-publication composition generalization, no independent laboratory/reactor validation; segmented boundary fixed at 725 C, not tuned by fold'
    if write:
        d=Path(out);d.mkdir(parents=True,exist_ok=True)
        (d/'group_holdout.json').write_text(json.dumps(results,indent=2)+'\n')
    return results

if __name__=='__main__':
    r=crossvalidate()
    print(json.dumps({k:v['pooled_MAE_um_s'] for k,v in r.items() if isinstance(v,dict)},indent=2))
