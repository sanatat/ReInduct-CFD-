"""Reference ideal-gas *sensible* enthalpy and Cp for CO, CO2, N2, H2, H.

CO/CO2/N2/H2: NIST Chemistry WebBook Shomate coefficients (Chase/JANAF).
H atom: Cp = 5/2 R, monatomic ideal gas assumption (sensible heat only).
T ref = 298.15 K. Chemical formation enthalpy is deliberately EXCLUDED;
H2 bond dissociation is separately accounted for in the radical energy ledger.
Not a calibrated high-temperature plasma thermodynamic model.
"""
from __future__ import annotations
import numpy as np
from v13.compressible_transport import R
SPECIES=('CO','CO2','N2','H2','H')
# t=T/1000; cp=A+B*t+C*t^2+D*t^3+E/t^2 ; h sensible (kJ/mol)
# Each: (Tswitch, low=[A,B,C,D,E,F,H], high=[...])
DATA={
 'CO':(1300.,[25.56759,6.096130,4.054656,-2.671301,.131021,-118.0089,-110.5271],
                [35.15070,1.300095,-.205921,.013550,-3.282780,-127.8375,-110.5271]),
 'CO2':(1200.,[24.99735,55.18696,-33.69137,7.948387,-.136638,-403.6075,-393.5224],
                 [58.16639,2.720074,-.492289,.038844,-6.447293,-425.9186,-393.5224]),
 'N2':(500.,[28.98641,1.853978,-9.647459,16.63537,.000117,-8.671914,0.],
                [19.50583,19.88705,-8.598535,1.369784,.527601,-4.935202,0.]),
 'H2':(1000.,[33.066178,-11.363417,11.432816,-2.772874,-.158558,-9.980797,0.],
                [18.563083,12.257357,-2.859786,.268238,1.977990,-1.147438,0.])}


def _evaluate(T, name, enthalpy=False):
    t=np.asarray(T,float)
    if np.any(t<298.15) or np.any(t>2000):
        raise ValueError('reference Cp/h evaluation requires 298.15–2000 K')
    if name=='H':
        return (2.5*R*(t-298.15)) if enthalpy else (np.zeros_like(t)+2.5*R)
    threshold,lo,hi=DATA[name]
    coeff=np.where((t[...,None]<threshold),np.asarray(lo),np.asarray(hi))
    a,b,c,d,e,f,h=np.moveaxis(coeff,-1,0)
    x=t/1000.
    if enthalpy:
        raw=1000*(a*x+b*x*x/2+c*x**3/3+d*x**4/4-e/x+f-h)
        # Tabulated Shomate fits are rounded independently on either side
        # of each breakpoint. Remove the (few J/mol) spurious enthalpy jump.
        z=threshold/1000.
        def integral(coef):
            A,B,C,D,E,F,H=coef
            return 1000*(A*z+B*z*z/2+C*z**3/3+D*z**4/4-E/z+F-H)
        offset=integral(lo)-integral(hi)
        return raw+np.where(t>=threshold,offset,0.)
    return a+b*x+c*x*x+d*x**3+e/x**2


def cp_molar(T):
    """Five-species array of Cp(T) [J/mol/K]; trailing dimensions match T."""
    return np.stack([_evaluate(T,s,False) for s in SPECIES])


def enthalpy_molar(T):
    """Five-species array of sensible enthalpies [J/mol] relative to 298.15 K."""
    return np.stack([_evaluate(T,s,True) for s in SPECIES])


def gas_energy(gas, Tgas, porosity):
    return porosity*np.sum(gas*enthalpy_molar(Tgas),axis=0)


def invert_energy(energy,gas,porosity,T_guess=None,tol=1e-9):
    energy=np.asarray(energy,float)
    if np.min(gas.sum(axis=0))<=0:raise ValueError('zero gas content')
    t=np.asarray(T_guess,float).copy() if T_guess is not None else np.full_like(energy,1050.)
    for _ in range(30):
        err=gas_energy(gas,t,porosity)-energy
        cp_tot=porosity*np.sum(gas*cp_molar(t),axis=0)
        t_next=t-err/np.maximum(cp_tot,1.e-16)
        if np.any((t_next<298.15)|(t_next>2000)):
            raise FloatingPointError('inverted temperature outside thermochemical table range')
        if np.max(np.abs(t_next-t))<tol:return t_next
        t=t_next
    raise RuntimeError('gas enthalpy inversion did not converge')
