#!/usr/bin/env python3
"""Conservative 1-D conjugate heat transfer across *separate* material meshes.

Gas layer x in [0,Lg] with Ng cells, solid layer [Lg,Lg+Ls] with Ns cells.
Harmonic face conductance at interface enforces equal-and-opposite heat flux.
Boundary temperatures prescribed; no advection/chemical reaction, benchmark only.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.linalg import solve_banded


def solve(Ng=24,Ns=32,Lg=.01,Ls=.02,kg=.05,ks=12.,TL=900.,TR=800.,Rcontact=0.):
    if min(Ng,Ns)<2 or min(Lg,Ls,kg,ks)<=0 or Rcontact<0 or TL<=0 or TR<=0:
        raise ValueError('bad multiregion thermal input')
    dxg=Lg/Ng;dxs=Ls/Ns
    dx=np.concatenate([np.full(Ng,dxg),np.full(Ns,dxs)])
    kval=np.concatenate([np.full(Ng,kg),np.full(Ns,ks)])
    N=Ng+Ns
    conduct=1./(dx[:-1]/(2*kval[:-1])+dx[1:]/(2*kval[1:])+np.where(np.arange(N-1)==Ng-1,Rcontact,0.))
    left=2*kg/dxg;right=2*ks/dxs
    diag=np.zeros(N);low=np.zeros(N);upper=np.zeros(N);rhs=np.zeros(N)
    diag[0]+=left;rhs[0]+=left*TL
    diag[-1]+=right;rhs[-1]+=right*TR
    for i,c in enumerate(conduct):
        diag[i]+=c;diag[i+1]+=c
        upper[i]=-c;low[i+1]=-c
    ab=np.zeros((3,N));ab[0,1:]=upper[:-1];ab[1,:]=diag;ab[2,:-1]=low[1:]
    T=solve_banded((1,1),ab,rhs)
    q_all=conduct*(T[:-1]-T[1:])
    qleft=left*(TL-T[0]);qright=right*(T[-1]-TR)
    qanalytic=(TL-TR)/(Lg/kg+Ls/ks+Rcontact)
    x=np.cumsum(dx)-dx/2
    summary={
        'Ng':Ng,'Ns':Ns,'interface_T_gas_K':float(T[Ng-1]),'interface_T_solid_K':float(T[Ng]),
        'heat_flux_left_W_m2':float(qleft),'heat_flux_right_W_m2':float(qright),
        'interface_heat_flux_W_m2':float(q_all[Ng-1]),
        'max_face_heat_flux_discrepancy_W_m2':float(np.max(abs(q_all-qanalytic))),
        'analytic_heat_flux_W_m2':float(qanalytic),
        'relative_flux_error':float(abs(qleft-qanalytic)/abs(qanalytic)),
        'limitation':'steady conduction-only separate-region interface; no OpenFOAM, 2D geometry, radiation or reacting gas'}
    return x,T,summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--Ng',type=int,default=24);p.add_argument('--Ns',type=int,default=32)
    p.add_argument('--out',default='results/v3/multiregion')
    a=p.parse_args();x,T,s=solve(Ng=a.Ng,Ns=a.Ns)
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    np.savetxt(out/'temperature.csv',np.column_stack([x,T]),delimiter=',',comments='',header='x_m,T_K')
    (out/'summary.json').write_text(json.dumps(s,indent=2)+'\n')
    print(json.dumps(s,indent=2))
if __name__=='__main__':main()
