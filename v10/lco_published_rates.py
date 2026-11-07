"""Reproducible published *kinetic coefficient* benchmark, not reactor validation.

Values transcribed from Nuraeni et al., Metall. Mater. Trans. B (2023)
DOI 10.1007/s11663-022-02712-1, Tables II & III. The late high-T
regimes stop when a third stage/deceleration begins, so the k values are
NOT total conversion curves at 240 min. All are from a single publication.
"""
from __future__ import annotations
import csv,json
from pathlib import Path
import numpy as np

SOURCE='10.1007/s11663-022-02712-1'
GAS_R=8.31446261815324
# T_C, stage limit minutes, rate constant min-1, R^2
EARLY=[(700,16,.00017,.998),(800,23,.00115,.996),(900,25,.00139,.976),
       (1000,19,.00605,.998),(1100,10,.01744,.987)]
LATE=[(700,240,.00017,.989),(800,240,.00041,.996),(900,126,.00102,.999),
      (1000,44,.002515,.983),(1100,28,.00505,.972)]

# Conditions and reported *recovered product yields* from a DIFFERENT study.
# Explicitly never compare recovered Li2CO3 mass yield to kinetic alpha.
EXTERNAL_QUALITATIVE={
 'source_doi':'10.1177/17475198211066533',
 'study':'She et al., Journal of Chemical Research (2022)',
 'atmosphere':'air, unlike 2023 argon experiment',
 'temperature_C':800,
 'hold_min':45,
 'carbon_to_LCO_mol_ratio':1.25,
 'reported_Li2CO3_recovery_fraction':.84,
 'reported_cobalt_oxide_recovery_fraction':.89,
 'observed_product_classes':['Li2CO3','cobalt_oxide'],
 'claim':'Independent study qualitative phase and recovery cross-check only; product recovery is NOT LiCoO2 reaction conversion.',
}

def arrhenius_fit(rows):
    temp=np.array([r[0]+273.15 for r in rows])
    kval=np.array([r[2] for r in rows])
    slope, intercept=np.polyfit(1/temp,np.log(kval),1)
    pred=np.exp(intercept+slope/temp)
    return {'E_a_J_per_mol':float(-slope*GAS_R),'predicted_k_min_1':pred.tolist(),
            'logRMSE':float(np.sqrt(np.mean((np.log(pred)-np.log(kval))**2)))}

def held_out_rate_test(rows):
    predictions=[]
    for idx,row in enumerate(rows):
        fit=arrhenius_fit([r for j,r in enumerate(rows) if j!=idx])
        t=row[0]+273.15
        train=np.array([r for j,r in enumerate(rows) if j!=idx])
        x=1/(train[:,0]+273.15);y=np.log(train[:,2]);slope,intercept=np.polyfit(x,y,1)
        kp=float(np.exp(intercept+slope/t))
        predictions.append({'held_out_temperature_C':row[0],
                            'measured_published_k_min_1':row[2],
                            'prediction_train_other_4_k_min_1':kp,
                            'log_error':float(np.log(kp/row[2])),
                            'relative_error':float(kp/row[2]-1)})
    return {'rows':predictions,'logMAE':float(np.mean([abs(r['log_error']) for r in predictions])),
            'median_abs_relative_error':float(np.median([abs(r['relative_error']) for r in predictions]))}

def build(outdir='v10/results'):
    p=Path(outdir);p.mkdir(parents=True,exist_ok=True)
    rows=[{'source_doi':SOURCE,'stage':stage,'temperature_C':t,'stage_end_min':end,
           'published_k_min_1':k,'published_fit_R2':r2,'measurement_type':'published model-fit constant (NOT raw conversion)'}
          for stage, seq in [('early_GB',EARLY),('late_A3',LATE)] for t,end,k,r2 in seq]
    with (p/'published_rate_table.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    result={'data_source':SOURCE,'data_type':'published 2023 stage-specific fitted rate constants; same-study LOOCV only',
            'early_fit':arrhenius_fit(EARLY),'late_fit':arrhenius_fit(LATE),
            'early_holdouts':held_out_rate_test(EARLY),'late_holdouts':held_out_rate_test(LATE),
            'independent_qualitative_study':EXTERNAL_QUALITATIVE,
            'external_reactor_conversion_validated':False,
            'validation_limitation':'Only rate constants from publication tables; second study is air-atmosphere product recovery not reaction conversion.'}
    (p/'rate_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    x=build();print(json.dumps({k:v for k,v in x.items() if k!='independent_qualitative_study'},indent=2))
