"""Spatial-map error and power convergence for v6 EM -> annular sector CFD cells.

Power scaling enforces exact integral, but local cell source may still be inaccurate.
A 3x3 Gauss-Legendre quadrature gives a higher-order cell-volume reference.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from v6.axisymmetric_em import run as solve_em,Grid
from v7.mapping import synthetic_annular_wedge,map_conservatively
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def analyse(out='v7/results/mapping_convergence'):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    _,em=solve_em(Grid(nr=48,nz=77))
    interp=RegularGridInterpolator((em['z'],em['r']),em['q'],bounds_error=False,fill_value=0)
    gauss,weights=np.polynomial.legendre.leggauss(3)
    records=[]
    for nr,nz in [(6,24),(12,48),(24,96),(36,144)]:
        C,V=synthetic_annular_wedge(nr=nr,nz=nz)
        rs=np.linspace(.005,.025,nr+1);zs=np.linspace(.04,.10,nz+1)
        refs=[]
        for j in range(nz):
            for i in range(nr):
                ra,rb=rs[i:i+2];za,zb=zs[j:j+2]
                ri=.5*(rb-ra)*gauss+.5*(ra+rb)
                zj=.5*(zb-za)*gauss+.5*(za+zb)
                zz,rr=np.meshgrid(zj,ri,indexing='ij')
                q=interp(np.column_stack((zz.ravel(),rr.ravel()))).reshape(3,3)
                num=np.sum(q*rr*(weights[:,None]*weights[None,:]))
                den=np.sum(rr*(weights[:,None]*weights[None,:]))
                refs.append(num/den)
        refs=np.asarray(refs)
        Pref=float(refs@V)
        qMapped,report=map_conservatively(em,C,V,Pref)
        weighted_l2=float(np.sqrt(np.sum(V*(qMapped-refs)**2)/np.sum(V*refs**2)))
        records.append({'nr':nr,'nz':nz,'power_W':Pref,'pre_scale_power_error_pct':100*report['raw_relative_power_error'],
                        'post_scale_power_error_pct':100*report['post_mapping_relative_power_error'],
                        'relative_local_L2_after_power_normalization':weighted_l2})
    (out/'mapping_convergence.json').write_text(json.dumps(records,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(7,4.4))
    ax.plot([x['nr'] for x in records],[x['relative_local_L2_after_power_normalization'] for x in records],'-o')
    ax.set(xlabel='Radial cells',ylabel='Volume-weighted relative L2 source error',title='EM-to-CFD local heat-source mapping refinement')
    ax.grid(alpha=.25);fig.tight_layout();fig.savefig(out/'mapping_convergence.png',dpi=160);plt.close(fig)
    return records
if __name__=='__main__':print(json.dumps(analyse(),indent=2))
