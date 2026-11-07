#!/usr/bin/env python3
"""0-D H2 -> 2 H -> H2 radical bookkeeping coupled to NiO reduction.

Reactions (molar bulk volumes):
 H2 -> 2H      (effective electron-impact dissociation, prescribed k_diss)
 2H -> H2      (effective recombination, k_rec)
 NiO + 2H -> Ni + H2O  (illustrative apparent surface reaction)

This is a stoichiometrically conservative SCREENING MODEL ONLY: no electron
energy distribution, discharge, plasma sheath, EM, gas hydrodynamics, or
quantitatively calibrated radical kinetics. Rates are illustrative.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp


def solve(end=12.,n0=150.,h2_initial=15.,k_diss=.06,k_recomb=.20,k_reduct=.03):
    if min(end,n0,h2_initial)<=0 or min(k_diss,k_recomb,k_reduct)<0:
        raise ValueError('invalid plasma-surrogate parameters')
    def rhs(t,y):
        c2,ch,cw,x=y
        rd=k_diss*max(c2,0.)
        rc=k_recomb*max(ch,0.)**2
        rr=k_reduct*max(ch,0.)**2*max(0,1.-x)
        return [-rd+rc,2*rd-2*rc-2*rr,rr,rr/n0]
    times=np.linspace(0,end,201)
    result=solve_ivp(rhs,(0,end),[h2_initial,0,0,0],t_eval=times,
                     method='BDF',rtol=1e-9,atol=1e-11)
    if not result.success: raise RuntimeError(result.message)
    C2,CH,CW,X=result.y
    hatoms=2*C2+CH+2*CW
    residual_h=float(np.max(abs(hatoms-2*h2_initial)))
    residual_o=float(np.max(abs(CW-n0*X)))
    if min(C2.min(),CH.min(),CW.min(),X.min())<-1e-10:
        raise RuntimeError('unphysical negative species')
    summary={'duration_s':end,'k_diss_s-1':k_diss,'k_rec_m3_mol_s':k_recomb,
             'k_red_m3_mol_s':k_reduct,'radical_max_mol_m3':float(CH.max()),
             'NiO_conversion_final':float(X[-1]),'H_atoms_balance_abs_mol_m3':residual_h,
             'solid_O_to_H2O_balance_abs_mol_m3':residual_o,
             'limitations':'prescribed effective kinetics; not plasma discharge or electron-energy simulation'}
    return result,summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--end',type=float,default=12.)
    p.add_argument('--kd',type=float,default=.06)
    p.add_argument('--out',default='results/v3/radicals')
    a=p.parse_args();r,s=solve(end=a.end,k_diss=a.kd)
    dest=Path(a.out);dest.mkdir(parents=True,exist_ok=True)
    np.savetxt(dest/'species.csv',np.column_stack((r.t,r.y.T)),delimiter=',',comments='',
               header='time_s,H2_mol_m3,H_mol_m3,H2O_mol_m3,NiO_conversion')
    (dest/'summary.json').write_text(json.dumps(s,indent=2)+'\n')
    print(json.dumps(s,indent=2))
if __name__=='__main__':main()
