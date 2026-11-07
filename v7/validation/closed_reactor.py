"""Independent zero-D *numerical conservation* benchmark for the OpenFOAM equations.

This benchmark is not an experiment and not an industrial process prediction.
Mass sources: -R/epsilon H2, +R/epsilon H2O, +R/site NiO conversion.
Energy source: q_em - deltaH*R, interphase hA(Ts-Tg) cancels in total.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
R=8.314462618

def solve(dt=.002,end_s=2.,q_ref_W_m3=18000.,feed_back=True):
    eps=.45;site=450.;Cg=eps*800.;Cs=1.05e6;hA=600.;dH=-2200.
    Tg=825.;Ts=850.;h2=11.;h2o=0.;X=0.
    M0=eps*h2+site*X;Mwo0=eps*h2o-site*X;E0=Cg*Tg+Cs*Ts
    cumulative=0.;hist=[]
    n=round(end_s/dt)
    assert abs(n*dt-end_s)<1e-12
    for i in range(n):
        q=q_ref_W_m3*(1+0.002*(Ts-850.))**(-1) if feed_back else q_ref_W_m3
        k=.0015*np.exp(-85000/R*(1/Ts-1/850.))
        rate=min(k*max(h2,0.)*max(1-X,0.),eps*h2/dt,site*(1-X)/dt)
        delta=rate*dt
        # Conservative paired source increments.
        h2-=delta/eps;h2o+=delta/eps;X+=delta/site
        dQg=hA*(Ts-Tg)*dt
        Tg+=dQg/Cg
        Ts+=((q-dH*rate)*dt-dQg)/Cs
        cumulative+=(q-dH*rate)*dt
        hist.append((i*dt+dt,h2,h2o,X,Ts,Tg,q,rate))
    mass_err=eps*h2+site*X-M0
    water_err=eps*h2o-site*X-Mwo0
    energy_err=(Cg*Tg+Cs*Ts)-(E0+cumulative)
    return {'n_steps':n,'max_Ts_K':Ts,'final_X':X,'final_H2_mol_m3_pore':h2,
            'H2_inventory_error_mol_m3_bulk':mass_err,
            'H2O_inventory_error_mol_m3_bulk':water_err,
            'energy_error_J_m3':energy_err,
            'initial_energy_J_m3':E0,'injected_energy_J_m3':cumulative,
            'qualification':'0D conservation benchmark; not a CFD/experimental validation'},np.array(hist)

def main(out='v7/results/conservation'):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    a,h=solve();b,_=solve(q_ref_W_m3=0)
    results={'heating_on':a,'heating_off':b,'temperature_feedback_effect_K':a['max_Ts_K']-b['max_Ts_K']}
    (out/'zeroD_conservation.json').write_text(json.dumps(results,indent=2)+'\n')
    np.savetxt(out/'zeroD_history.csv',h,delimiter=',',header='time_s,H2_mol_m3,H2O_mol_m3,X,Ts_K,Tg_K,q_W_m3,R_mol_m3_s',comments='')
    return results
if __name__=='__main__':print(json.dumps(main(),indent=2))
