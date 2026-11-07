"""Harmonic 1-D magnetic diffusion for temperature-dependent conductivity.

    d/dy[(1/sigma(T)) dH/dy] = i omega mu H

H(0)=H0, H'(L)=0. Linear MQS slab; no actual coil geometry,
Lorentz-force feedback or plasma chemistry. Conductivity is mapped
from an axially averaged solid-temperature profile by the reactor.
"""
from __future__ import annotations
import numpy as np
from .induction1d import MU0


def solve_variable(sigma, *, length=0.02, frequency=100_000.0,
                   h0=10_000.0, mu_r=1.0):
    sigma=np.asarray(sigma,dtype=float)
    n=sigma.size
    if n<5 or np.any(~np.isfinite(sigma)) or np.any(sigma<=0) or min(length,frequency,h0,mu_r)<=0:
        raise ValueError('need n>=5 and positive finite electromagnetic parameters')
    dy=length/(n-1)
    y=np.linspace(0,length,n)
    invsigma=1.0/sigma
    invface=2*invsigma[1:]*invsigma[:-1]/(invsigma[1:]+invsigma[:-1])
    low=np.zeros(n,complex); diag=np.zeros(n,complex); up=np.zeros(n,complex)
    b=np.zeros(n,complex)
    diag[0]=1.;b[0]=h0
    for j in range(1,n-1):
        low[j]=invface[j-1]/dy**2
        up[j]=invface[j]/dy**2
        diag[j]=-low[j]-up[j]-1j*(2*np.pi*frequency)*MU0*mu_r
    low[-1]=-1.; diag[-1]=1.
    for j in range(1,n):
        w=low[j]/diag[j-1]
        diag[j]-=w*up[j-1]
        b[j]-=w*b[j-1]
    h=np.zeros(n,complex)
    h[-1]=b[-1]/diag[-1]
    for j in range(n-2,-1,-1):
        h[j]=(b[j]-up[j]*h[j+1])/diag[j]
    current=-np.gradient(h,dy,edge_order=2)
    current[-1]=0.
    q=np.abs(current)**2/(2*sigma)
    return {'y':y,'H':h,'J':current,'q':q,
            'power_per_area_W_m2':float(np.trapezoid(q,y))}


def conductivity_from_temperature(T, sigma_ref=1e6, alpha=0.003,
                                  reference_K=850.):
    if sigma_ref<=0 or alpha<0:
        raise ValueError('conductivity reference must be positive; alpha nonnegative')
    t=np.asarray(T,dtype=float)
    denom=1+alpha*(t-reference_K)
    if np.any(denom<=0):
        raise ValueError('conductivity model invalid outside temperature range')
    return sigma_ref/denom
