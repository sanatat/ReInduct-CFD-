#!/usr/bin/env python3
"""Sequential EM-to-thermal coupling for a conducting slab, implicit FV heat equation.

rhoCp dT/dt = k d²T/dx² + q_EM(x) - hvol*(T-Tgas).
Zero-heat-flux at both ends; gas is a prescribed constant-temperature reservoir.
This is a ONE-WAY coupling (sigma constant), not CFD, MHD, or plasma.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
try:
    from .induction1d import solve
except ImportError:
    from induction1d import solve

def thermal(em, duration=20., dt=0.25, rho_cp=3.0e6, k_thermal=18.,
            hvol=7500., tgas=850.):
    x=em['x']; n=len(x); dx=x[1]-x[0]
    alpha=k_thermal*dt/(rho_cp*dx*dx)
    beta=dt*hvol/rho_cp
    T=np.full(n,tgas,dtype=float)
    # flux conservative finite-volume with insulated ends
    a=np.full(n,-alpha);b=np.full(n,1+2*alpha+beta);c=np.full(n,-alpha)
    a[0]=0;c[0]=-alpha;b[0]=1+alpha+beta
    c[-1]=0;a[-1]=-alpha;b[-1]=1+alpha+beta
    for step in range(int(round(duration/dt))):
        rhs=T+dt*(em['q']+hvol*tgas)/rho_cp
        bb=b.copy();rr=rhs.copy()
        for j in range(1,n):
            w=a[j]/bb[j-1];bb[j]-=w*c[j-1];rr[j]-=w*rr[j-1]
        T[-1]=rr[-1]/bb[-1]
        for j in range(n-2,-1,-1):T[j]=(rr[j]-c[j]*T[j+1])/bb[j]
    return T

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',default='results/induction_thermal')
    p.add_argument('--h0',type=float,default=10000.)
    p.add_argument('--duration',type=float,default=20.)
    args=p.parse_args()
    e=solve(h0=args.h0)
    T=thermal(e,duration=args.duration)
    path=Path(args.out);path.mkdir(parents=True,exist_ok=True)
    with (path/'fields.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['x_m','q_W_m3','T_K'])
        w.writerows(zip(e['x'],e['q'],T))
    d={'Tmax_K':float(T.max()),'Tmean_K':float(T.mean()),
       'Tmin_K':float(T.min()),'skin_depth_m':e['delta'],
       'EM_input_power_per_area_W_m2':e['power_per_area'],
       'caveat':'one-way EM/thermal, prescribed gas sink, illustrative thermal properties'}
    (path/'summary.json').write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps(d,indent=2))
if __name__=='__main__':main()
