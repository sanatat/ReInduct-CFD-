"""Independent-study *qualitative* applicability audit for LiCoO2+C.

No artificial conversion time-series is generated. Li recovery/CoO recovery
are not the same observable as extent of LiCoO2 reduction alpha.
"""
from __future__ import annotations
import json
from pathlib import Path
from v10.lco_published_rates import EXTERNAL_QUALITATIVE
from v9.lco_process import scenario

def audit(out='v10/results/external_phase_audit.json'):
    source=EXTERNAL_QUALITATIVE
    prediction=scenario(temp_c=source['temperature_C'],time_min=source['hold_min'],
                        lco_mol=4,carbon_mol=5)
    predicted_solids={sp for sp, mol in prediction['inventories_mol'].items()
                      if mol>1e-10 and sp not in ('C','LiCoO2') and sp not in ('CO','CO2')}
    reported=set(source['observed_product_classes'])
    mismatch='Co' in predicted_solids and 'cobalt_oxide' in reported
    res={'independent_publication':source,
         'v9_model_predicted_solids':sorted(predicted_solids),
         'reported_product_classes':sorted(reported),
         'metallic_cobalt_vs_cobalt_oxide_mismatch':mismatch,
         'lithium_carbonate_phase_compatible':'Li2CO3' in predicted_solids and 'Li2CO3' in reported,
         'same_atmosphere':False,'same_reactor':False,
         'can_quantitatively_validate_conversion':False,
         'status':'QUALITATIVE PARTIAL CONSISTENCY AND LIMITATION; NOT A REACTOR VALIDATION',
         'explanation':'The 2022 experiment uses air and separately reports extracted product recoveries; the 2023 kinetics source uses inert gas. CoO recovery at 800 C must not be conflated with Co metal conversion or alpha.',
         'next_experiment':'Obtain inert-atmosphere LiCoO2+C conversion vs time at controlled carbon ratio, sampling phase-specific Co/CoO and Li2CO3/ Li2O by XRD and offgas analysis.'}
    p=Path(out);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(res,indent=2)+'\n')
    return res
if __name__=='__main__':print(json.dumps(audit(),indent=2))
