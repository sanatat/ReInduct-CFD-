"""Explicitly separate analytical verification, same-study holdout,
and unavailable independent industrial validation."""
from __future__ import annotations
import json
from pathlib import Path
from v8.lco_holdout import estimate
from v8.lco_kinetics import OBSERVED_T240_APPROX
from v8.physics import R
import numpy as np
from scipy.optimize import least_squares
from v8.lco_holdout import alpha

SRC='https://doi.org/10.1007/s11663-022-02712-1'
OTHER='https://pmc.ncbi.nlm.nih.gov/articles/PMC9548595/'

def leave_one_temperature_out():
    cases=[]
    temps=sorted(OBSERVED_T240_APPROX)
    for test in temps:
        train=tuple(x for x in temps if x!=test)
        # stage activation energies fixed at published estimates; 2 prefactors
        # fitted on 3 APPROXIMATE endpoints, evaluated on remaining one.
        tr=estimate(train_C=train,test_C=(test,))
        value=tr['held_out']['rows'][0]
        cases.append({'held_out_temperature_C':test,'training_temperature_C':list(train),
                      'observed_conversion_APPROX':value['observed_alpha_240min_approx'],
                      'predicted_conversion':value['predicted_alpha_240min'],
                      'abs_error':abs(value['predicted_alpha_240min']-value['observed_alpha_240min_approx'])})
    return {'source':SRC,'data_quality':'approximate 240-min endpoints interpreted from source; NOT digitized raw measurements',
            'validation_level':'within-study leave-one-temperature-out (4 folds), NOT independent lab/reactor validation',
            'summary_MAE_conversion':float(np.mean([x['abs_error'] for x in cases])),
            'cases':cases,'independent_reactor_validation':False,
            'independent_literature_context':{'source':OTHER,
                    'cathode_system':'spent lithium-ion mixed cathode materials reduced in H2',
                    'published_activation_energy_500_700_C_kJ_mol':84.86,
                    'caveat':'not numerically comparable with LiCoO2/C, which has different reactants and mechanisms'}}

def main(out='v9/results/kinetics_protocol.json'):
    d=leave_one_temperature_out();p=Path(out);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(d,indent=2)+'\n');return d

if __name__=='__main__':print(json.dumps(main(),indent=2))
