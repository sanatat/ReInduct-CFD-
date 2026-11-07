"""Fail-closed *independent* reactor validation harness.

No measured reactor time series were provided. This script will NOT create
synthetic observations or count the original fitting study as independent.
It requires independently sourced experimental measurements for the SAME
LiCoO2 + C chemical system; data from NMC + H2 must not be mixed.
"""
from __future__ import annotations
import csv,json,math
from pathlib import Path
import numpy as np
from v9.lco_process import scenario
SOURCE_DOI='10.1007/s11663-022-02712-1'
FIELDS=['source_doi','experiment_id','reactor_id','chemistry','temperature_C',
        'time_min','carbon_mol_per_lco_mol','observed_conversion','measurement_uncertainty_abs']

def evaluate(path,output=None):
    with open(path,newline='') as f:
        reader=csv.DictReader(f)
        if not set(FIELDS).issubset(reader.fieldnames or []):
            raise ValueError('Missing required experimental provenance, reactor, or measurements columns')
        rows=list(reader)
    if len(rows)<3:raise ValueError('Require >=3 real measured points from an independent reactor study')
    if len({r['experiment_id'] for r in rows})<2:
        raise ValueError('Require >=2 independent experiment IDs')
    predictions=[]
    for row in rows:
        doi=row['source_doi'].lower().strip().replace('https://doi.org/','').replace('doi:','')
        if not doi or doi==SOURCE_DOI:
            raise ValueError('Original kinetic source is not independent external validation')
        if row['chemistry'].strip()!='LiCoO2+C':
            raise ValueError('Mismatched chemistry: only LiCoO2+C is supported')
        t=float(row['temperature_C']);time=float(row['time_min'])
        y=float(row['observed_conversion']);unc=float(row['measurement_uncertainty_abs'])
        ratio=float(row['carbon_mol_per_lco_mol'])
        if not math.isfinite(ratio) or ratio<0:raise ValueError('Invalid carbon-to-LCO molar ratio')
        if not all(math.isfinite(x) for x in [t,time,y,unc]):raise ValueError('Invalid data')
        if not 0<=y<=1 or unc<=0:raise ValueError('Invalid measured conversion or uncertainty')
        s=scenario(t,time,10.,10.*ratio)
        yp=s['carbon_limited_conversion_LCO']
        predictions.append({'source_doi':doi,'experiment_id':row['experiment_id'],
                            'reactor_id':row['reactor_id'],
                            'observed':y,'predicted':yp,'abs_error':abs(yp-y),
                            'within_observed_uncertainty':abs(yp-y)<=unc})
    errors=np.array([p['predicted']-p['observed'] for p in predictions])
    result={'test_type':'independent source reactor time-series vs previously fixed assumed kinetic parameters',
            'n_measured':len(predictions),'MAE_conversion_fraction':float(np.abs(errors).mean()),
            'RMSE_conversion_fraction':float(np.sqrt(np.mean(errors**2))),
            'mean_bias':float(errors.mean()),
            'fraction_within_measurement_uncertainty':float(np.mean([p['within_observed_uncertainty'] for p in predictions])),
            'limitations':'not comparable if carbon feed ratio/reactor conditions differ from illustrative defaults; extend metadata before publication',
            'individual':predictions}
    if output:
        p=Path(output);p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--csv',required=True);ap.add_argument('--out',default='v9/results/independent_validation.json')
    args=ap.parse_args();print(json.dumps(evaluate(args.csv,args.out),indent=2))
