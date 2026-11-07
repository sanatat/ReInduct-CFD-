"""Exploratory lithium cobalt oxide / carbon conversion curves.

The stage models and activation energies are taken from LiCoO2 carbothermal
kinetics literature (2023); all prefactors/reference rate constants and transition
choice in this module are *assumed*, not experimentally estimated or validated.
Measured 240-minute conversions are used only as separate comparison anchors.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq
from v8.physics import R

EA_GB_J_MOL=121e3
EA_A3_J_MOL=95e3
TREF_K=1173.15
# Assumed reference constants (not fitted to experimental curves).
KREF_GB_PER_S=2.0e-5
KREF_A3_PER_S=1.1e-4
TRANSITION_MIN=25.
OBSERVED_T240_APPROX={700:.19,800:.54,1000:.92,1100:.94}

def ginstling_brounshtein(x):
    if not 0<=x<1:raise ValueError('X outside [0,1)')
    return 1-2*x/3-(1-x)**(2/3)

def gb_inverse(value):
    if value<0 or value>=1/3:raise ValueError('GB integral outside [0,1/3)')
    upper=1-1e-12
    if value>=ginstling_brounshtein(upper): return upper
    return brentq(lambda x:ginstling_brounshtein(x)-value,0.,upper)

def arrhenius(k_ref,Ea_J_mol,T):
    if T<=0:raise ValueError('T <=0')
    return k_ref*np.exp(-Ea_J_mol/R*(1./T-1./TREF_K))

def lco_conversion(t_s,T_K,transition_s=TRANSITION_MIN*60):
    if not 973.15<=T_K<=1373.15 or t_s<0:
        raise ValueError('Only 700–1100C and nonnegative time: no extrapolation')
    k1=arrhenius(KREF_GB_PER_S,EA_GB_J_MOL,T_K)
    k2=arrhenius(KREF_A3_PER_S,EA_A3_J_MOL,T_K)
    if t_s<=transition_s:
        return gb_inverse(min(k1*t_s,1/3-1e-12))
    x0=gb_inverse(min(k1*transition_s,1/3-1e-12))
    a0=(-np.log1p(-x0))**(1/3)
    # Continuity of conversion at transition. Stage rates need not be continuous.
    return float(1-np.exp(-(a0+k2*(t_s-transition_s))**3))

def benchmark():
    rows=[]
    for tc,observed in OBSERVED_T240_APPROX.items():
        pred=lco_conversion(240*60,tc+273.15)
        rows.append({'temperature_C':tc,'observed_alpha_240min_approx':observed,
                     'unfitted_model_prediction':pred,'error':pred-observed})
    return {'source':'https://pmc.ncbi.nlm.nih.gov/articles/PMC9821372/',
            'source_activation_energies_kJ_mol':[121,95],
            'reference_prefactors_calibrated':False,
            'qualifications':'Observed endpoints approximate, extracted from article prose; not external independent validation.',
            'rows':rows,'MAE_unfitted':float(np.mean([abs(x['error']) for x in rows]))}

if __name__=='__main__':
    import json
    print(json.dumps(benchmark(),indent=2))
