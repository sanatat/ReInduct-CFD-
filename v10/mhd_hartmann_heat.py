"""Physically dimensioned plane Hartmann flow, induced Joule heat & viscous heat.

Fixed transverse magnetic field B, incompressible weakly conductive Newtonian
liquid, no imposed cross-channel electric field: mu*u''-sigma*B^2*u+G=0.
Temperature solves -k*T'' = mu*(u')^2 + sigma*B^2*u^2.
No plasma ionization and no reactive bed in this benchmark.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.linalg import solve_banded

def solve(n=101, a=.01, mu=.01, sigma=80., B=3., pressure_gradient=100., kappa=2.0, wall_kelvin=300., out=None):
    if n<10 or min(a,mu,sigma,kappa)<=0 or B<0:raise ValueError('Bad physical parameters')
    y=np.linspace(-a,a,n);dy=y[1]-y[0]
    band=np.zeros((3,n-2));band[1]=2*mu/dy**2+sigma*B**2
    band[0,1:]=-mu/dy**2;band[2,:-1]=-mu/dy**2
    u=np.zeros(n);u[1:-1]=solve_banded((1,1),band,np.full(n-2,pressure_gradient))
    Ha=a*B*np.sqrt(sigma/mu)
    if Ha<1e-8:exact=pressure_gradient/(2*mu)*(a*a-y*y)
    else:exact=pressure_gradient/(sigma*B*B)*(1-np.cosh(Ha*y/a)/np.cosh(Ha))
    # Face gradients evaluate exact discrete viscous wall dissipation budget.
    grad=(u[1:]-u[:-1])/dy
    visc_integrated=mu*np.sum(grad**2)*dy
    joule_local=sigma*B*B*u*u
    joule_integrated=np.trapezoid(joule_local,y)
    power_in=pressure_gradient*np.trapezoid(u,y)
    mechanical_error=power_in-visc_integrated-joule_integrated
    # Approximate cell-centered heat sources with centered velocity gradient.
    face_heat=mu*np.gradient(u,dy,edge_order=2)**2+joule_local
    heatband=np.zeros((3,n-2));heatband[1]=2*kappa/dy**2
    heatband[0,1:]=-kappa/dy**2;heatband[2,:-1]=-kappa/dy**2
    temperature=np.full(n,wall_kelvin)
    temperature[1:-1]+=solve_banded((1,1),heatband,face_heat[1:-1])
    error=np.linalg.norm(u-exact)/max(np.linalg.norm(exact),1e-12)
    report={'n':n,'Hartmann_number':float(Ha),'relative_L2_flow_error':float(error),
            'pressure_input_W_m2':float(power_in),'viscous_dissipation_W_m2':float(visc_integrated),
            'Joule_heating_W_m2':float(joule_integrated),
            'mechanical_budget_relative_error':float(mechanical_error/max(abs(power_in),1e-12)),
            'max_temperature_K':float(np.max(temperature)),
            'assumptions':'1D steady laminar Hartmann flow, uniform imposed B, constant conductivity; not a plasma or industrial reactor'}
    if out:
        p=Path(out);p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(json.dumps(report,indent=2)+'\n')
        np.savez(p.with_suffix('.npz'),y=y,u=u,exact=exact,temperature=temperature,joule_W_m3=joule_local)
    return report

if __name__=='__main__':
    print(json.dumps(solve(out='v10/results/hartmann.json'),indent=2))
