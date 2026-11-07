"""Pressure/temperature Picard correction for the *reduced-order* porous gas.

The implicit total-molar pressure equation and the conservative species/EOS
field are iterated with sensible-energy advection to a consistent end-step
T and p. Properties are still frozen within an iteration and do NOT represent
compressible Navier--Stokes or thermo-consistent reaction energetics.
"""
from __future__ import annotations
import numpy as np
from v13.compressible_transport import pressure_corrected_step, R
from v13.gas_properties import enthalpy_molar, invert_energy


def advance(gas, Tgas, Eg, source_bulk, extra_energy_bulk, *, dt, length,
            porosity, pin, pout, inlet_fraction, T_in,
            max_iterations=60, relax=0.7, temperature_tolerance_K=2.e-8,
            pressure_tolerance_Pa=1.e-5, pressure_relative_tolerance=1.e-9):
    n=gas.shape[1]
    dz=length/n
    guess=Tgas.copy()
    final=None
    for count in range(1,max_iterations+1):
        trial, flux, sflux, info=pressure_corrected_step(
            gas,guess,source_bulk,dt,length,porosity,pin,pout,
            inlet_fraction=inlet_fraction)
        heat_flux=np.empty(n+1)
        for face, J in enumerate(flux):
            if face==0: up=T_in if J>=0 else Tgas[0]
            elif face==n:up=Tgas[-1] if J>=0 else T_in
            else:up=Tgas[face-1] if J>=0 else Tgas[face]
            heat_flux[face]=float(np.dot(sflux[:,face],enthalpy_molar(up).reshape(-1)))
        e_next=Eg+dt*(-np.diff(heat_flux)/dz+extra_energy_bulk)
        if np.any(e_next<=0) or np.any(trial.sum(axis=0)<=0):
            raise FloatingPointError('nonphysical coupled gas energy/concentration')
        t_next=invert_energy(e_next,trial,porosity,T_guess=guess)
        delta=float(np.max(abs(t_next-guess)))
        eos=float(np.max(np.abs(trial.sum(axis=0)*R*t_next-info['pressure_Pa'])))
        final=(trial,flux,sflux,info,e_next,t_next,heat_flux,
               {'iterations':count,'temperature_iteration_error_K':delta,
                'pressure_EOS_residual_Pa':eos})
        effective_p_tol=max(pressure_tolerance_Pa,
                             pressure_relative_tolerance*max(1.,float(np.max(np.abs(info['pressure_Pa'])))))
        if delta<=temperature_tolerance_K and eos<=effective_p_tol:
            return final
        guess=relax*t_next+(1.-relax)*guess
    raise RuntimeError(f'Pressure/temperature iteration did not converge after {max_iterations} iterations; '
                       f'last dT={final[-1]["temperature_iteration_error_K"]:.4g} K, '
                       f'EOS residual={final[-1]["pressure_EOS_residual_Pa"]:.4g} Pa')
