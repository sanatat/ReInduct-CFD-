#!/usr/bin/env python3
"""1-D frequency-domain magnetic diffusion in a conducting slab.

Solve d²H/dx² = i*omega*mu*sigma*H on x in [0,L],
H(0)=H0 (imposed tangential magnetic-field amplitude), dH/dx(L)=0.
J = -dH/dx (A/m²), q = |J|²/(2*sigma) (W/m³).
This is a planar, linear, magnetoquasistatic benchmark, not a coil or plasma model.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np

MU0 = 4 * math.pi * 1e-7

def solve(n=301, length=0.02, frequency=100_000., conductivity=1e6,
          mu_r=1.0, h0=10_000.):
    if n < 5 or min(length,frequency,conductivity,mu_r,h0) <= 0:
        raise ValueError('positive physical inputs and n>=5 required')
    dx=length/(n-1)
    x=np.linspace(0,length,n)
    mu=MU0*mu_r
    omega=2*math.pi*frequency
    a=np.zeros(n,dtype=complex)
    b=np.zeros(n,dtype=complex)
    c=np.zeros(n,dtype=complex)
    d=np.zeros(n,dtype=complex)
    b[0]=1.; d[0]=h0
    k2=1j*omega*mu*conductivity
    for j in range(1,n-1):
        a[j]=1.; b[j]=-2.-k2*dx**2; c[j]=1.
    a[-1]=-1.; b[-1]=1. # zero field gradient at insulated symmetry boundary
    # Thomas tridiagonal algorithm with complex coefficients
    for j in range(1,n):
        f=a[j]/b[j-1]; b[j]-=f*c[j-1]; d[j]-=f*d[j-1]
    H=np.empty(n,dtype=complex)
    H[-1]=d[-1]/b[-1]
    for j in range(n-2,-1,-1): H[j]=(d[j]-c[j]*H[j+1])/b[j]
    J=-np.gradient(H,dx,edge_order=2)
    J[-1]=0.0
    q=np.abs(J)**2/(2*conductivity)
    delta=math.sqrt(2/(omega*mu*conductivity))
    return {'x':x,'H':H,'J':J,'q':q,'delta':delta,
            'power_per_area':float(np.trapezoid(q,x)),
            'params':{'frequency_Hz':frequency,'sigma_S_m':conductivity,
                      'mu_r':mu_r,'H0_A_m':h0,'length_m':length,'n':n}}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--nx',type=int,default=301)
    p.add_argument('--frequency',type=float,default=100000)
    p.add_argument('--conductivity',type=float,default=1e6)
    p.add_argument('--h0',type=float,default=10000)
    p.add_argument('--length',type=float,default=0.02)
    p.add_argument('--out',default='results/induction')
    args=p.parse_args()
    r=solve(args.nx,args.length,args.frequency,args.conductivity,h0=args.h0)
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    with (out/'field.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['x_m','H_real_A_m','H_imag_A_m','J_magnitude_A_m2','q_W_m3'])
        for x,h,j,q in zip(r['x'],r['H'],r['J'],r['q']):
            w.writerow([x,h.real,h.imag,abs(j),q])
    summary={**r['params'],'skin_depth_m':r['delta'],'power_per_area_W_m2':r['power_per_area'],
             'model':'linear planar magnetoquasistatic diffusion; no full coil or MHD'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
