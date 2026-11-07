"""Illustrative NiO + H2 -> Ni + H2O model (NOT literature-calibrated).

K(T) is ideal-gas reaction equilibrium proxy for unit-activity solids.
For gas quotient Q = pH2O/pH2, forward flux is suppressed as Q approaches K.
The apparent Arrhenius rate is not a validated intrinsic kinetic law.
"""
import numpy as np
R=8.314462618

def equilibrium_constant(T, dh_ref=-2200.0, ds_ref=50.0, tref=850.0, delta_cp=0.):
    t=np.asarray(T,dtype=float)
    if np.any(t<=0): raise ValueError('temperature must be positive')
    dh=dh_ref+delta_cp*(t-tref)
    ds=ds_ref+delta_cp*np.log(t/tref)
    dg=dh-t*ds
    return np.exp(np.clip(-dg/(R*t),-500,500))

def reaction_rate(T, h2, h2o, conversion, *, n0=150., preexp=1500., Ea=85000.,
                  c_ref=15., dh_ref=-2200., ds_ref=50., delta_cp=0.):
    T=np.asarray(T,dtype=float)
    h2=np.maximum(0,np.asarray(h2,dtype=float))
    h2o=np.maximum(0,np.asarray(h2o,dtype=float))
    conversion=np.asarray(conversion,dtype=float)
    if np.any(T<=0) or min(n0,c_ref)<=0 or preexp<0 or Ea<0:
        raise ValueError('invalid kinetic parameters')
    eq=equilibrium_constant(T,dh_ref,ds_ref,delta_cp=delta_cp)
    driving=np.maximum(h2-h2o/eq,0.)/c_ref
    return n0*preexp*np.exp(-Ea/(R*T))*driving*np.maximum(1.-conversion,0.)
