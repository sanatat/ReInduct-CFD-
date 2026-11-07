"""NiO + H2 ⇌ Ni + H2O thermodynamic feasibility calculator.

Sources:
 - OECD NEA / CRC: NiO standard formation enthalpy -239.7 kJ mol-1,
   entropy 37.99 J mol-1 K-1.
 - NIST JANAF: Ni(cr) S°298=29.870 J mol-1 K-1.
 - NIST WebBook Shomate: H2(g), H2O(g) Cp,H,S.
 - JAEA NiO and Ni(fcc) thermodynamic tables: sensible H and S increments
   relative to 298.15 K, re-anchored to consistent standard values above.
This combines experimental compilations and DOES NOT replace a validated
thermodynamic software or gas/solid equilibrium phase database.
Temperature interpolation only 500–1000 K in this example.
"""
from pathlib import Path
import numpy as np
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=8.314462618
# NIST gas phase Shomate valid H2O 500-1700 K, H2 298-1000 K.
H2={'A':33.066178,'B':-11.363417,'C':11.432816,'D':-2.772874,'E':-.158558,
    'F':-9.980797,'G':172.707974,'H':0.}
H2O={'A':30.09200,'B':6.832514,'C':6.793435,'D':-2.534480,'E':.082139,
     'F':-250.8810,'G':223.3967,'H':-241.8264}
# JAEA tabulated condensed-phase H (J/mol), S (J/mol/K)
T_TAB=np.array([298.15,300.,400.,500.,600.,700.,800.,900.,1000.])
NIO_H=np.array([-236915.,-236837.,-232406.,-227639.,-222662.,-217538.,-212288.,-206943.,-201509.])
NIO_S=np.array([43.81885,44.07750,56.79430,67.42147,76.49167,84.39335,91.39650,97.69134,103.41524])
NI_H=np.array([-.0055,47.65638,2723.35551,5648.01813,8970.88449,12333.08952,15409.56575,18512.77408,21681.94551])
NI_S=np.array([29.79545,29.95482,37.63816,44.15197,50.19563,55.39565,59.50388,63.15832,66.49667])

def shomate_h_s(T,coef):
    t=np.asarray(T)/1000.
    A,B,C,D,E,F,G,H=[coef[k] for k in 'ABCDEFGH']
    # Enthalpy relative 298.15 K (kJ/mol); entropy J/mol/K
    return (A*t+B*t*t/2+C*t**3/3+D*t**4/4-E/t+F-H)*1000, \
           A*np.log(t)+B*t+C*t*t/2+D*t**3/3-E/(2*t*t)+G

def thermo(T):
    T=np.asarray(T,dtype=float)
    if np.any((T<500)|(T>1000)):raise ValueError('this curated range supports 500 <= T <= 1000 K only')
    # Anchor solid enthalpy to consistent 298-K reaction heat; use sensible increments.
    nio_h=-239700.+np.interp(T,T_TAB,NIO_H)-NIO_H[0]
    ni_h=np.interp(T,T_TAB,NI_H)-NI_H[0]
    nio_s=37.99+np.interp(T,T_TAB,NIO_S)-NIO_S[0]
    ni_s=29.87+np.interp(T,T_TAB,NI_S)-NI_S[0]
    h2_h,h2_s=shomate_h_s(T,H2)
    h2o_h,h2o_s=shomate_h_s(T,H2O)
    h2o_h+=-241826.0
    dH=ni_h+h2o_h-nio_h-h2_h
    dS=ni_s+h2o_s-nio_s-h2_s
    dG=dH-T*dS
    logK=-dG/(R*T)
    return {'delta_H_Jmol':dH,'delta_S_JmolK':dS,'delta_G_Jmol':dG,'Kp':np.exp(logK)}

def affinity_factor(T,p_h2,p_h2o):
    """Dimensionless kinetic driving force. Partial pressure ratio Q=pH2O/pH2."""
    K=thermo(T)['Kp']
    return 1.-(np.asarray(p_h2o)/np.asarray(p_h2))/K

def report(out='v5/results/thermo'):
    folder=Path(out);folder.mkdir(parents=True,exist_ok=True)
    Ts=np.linspace(500,1000,101);props=thermo(Ts)
    data={'reaction':'NiO(s) + H2(g) <=> Ni(s) + H2O(g)',
          'temperature_K':[600,800,1000],
          'values':[{k:float(thermo(float(t))[k]) for k in props} for t in [600,800,1000]],
          'limitations':['Ideal gas/solid activities; no Ni/NiO solid solution or particle pore model',
                         'Gas Shomate parameter ranges and rounded solid tabular interpolation',
                         'NiO reaction parameters not transferred to LiCoO2 chemistry',
                         'Not fitted to laboratory thermodynamic equilibrium data']}
    (folder/'equilibrium.json').write_text(json.dumps(data,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(7,4))
    ax.plot(Ts,props['delta_G_Jmol']/1000,label='ΔrG° [kJ/mol]')
    ax.plot(Ts,props['delta_H_Jmol']/1000,label='ΔrH° [kJ/mol]')
    ax.set_xlabel('Temperature [K]');ax.set_ylabel('Reaction property [kJ/mol]')
    ax.axhline(0,color='gray',lw=.7);ax.legend();ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(folder/'nio_thermodynamics.png',dpi=160);plt.close(fig)
    print(json.dumps(data,indent=2))
    return data

if __name__=='__main__':report()
