"""Low-magnetic-Reynolds-number *motional* Lorentz proxy for conducting gas.

For an out-of-plane electric-current closure (E_z=0) and applied harmonic
in-plane B, J_z = sigma_g (u_x B_y - u_y B_x). The time-averaged force is
F = 0.5 Re(J_z x conj(B)). This is a limited MHD source term, NOT
self-consistent H2 plasma (no plasma conductivity prediction, charged
species or electron energy distribution). No electrically insulated endcaps.
"""
from __future__ import annotations
import numpy as np

def motional_lorentz(ux,uy,Bx,By,sigma_g=3.,rho_g=.10):
    ux=np.asarray(ux,dtype=float);uy=np.asarray(uy,dtype=float)
    Bx=np.asarray(Bx,complex);By=np.asarray(By,complex)
    if sigma_g<0 or rho_g<=0:raise ValueError('invalid gas properties')
    jz=sigma_g*(ux*By-uy*Bx)
    ax=-.5*np.real(jz*np.conj(By))/rho_g
    ay=.5*np.real(jz*np.conj(Bx))/rho_g
    dissipated_W_m3=rho_g*(ax*ux+ay*uy)
    return ax,ay,dissipated_W_m3
