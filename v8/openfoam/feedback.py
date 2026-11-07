"""Offline/OpenFOAM-file EM<->temperature field exchange for annular wedge cells.

Reads actual cell geometry and a scalar Ts field from an ASCII OpenFOAM output.
Recomputes EM electrical conductivity from Ts, then power-conservatively maps qEM.
This tool is independently tested on synthetic OpenFOAM-like fields; the coupling
has NOT been exercised with a compiled OpenFOAM runtime in this environment.
"""
from __future__ import annotations
import argparse,json,re
from pathlib import Path
import numpy as np
from scipy.interpolate import griddata
from v6.axisymmetric_em import Grid, run as em_solve
from v7.mapping import load_actual_cells,map_conservatively,write_foam_scalar,em_ref_power


def read_foam_internal_scalar(file,expected_count):
    """Fail closed on unsupported OpenFOAM formats (binary/vector/invalid counts)."""
    src=Path(file).read_text()
    if re.search(r'\bformat\s+binary\s*;',src):raise ValueError('Only ASCII scalar fields supported')
    m=re.search(r'\binternalField\s+uniform\s+([+\-\d.eE]+)\s*;',src)
    if m:
        v=float(m.group(1))
        if not np.isfinite(v):raise ValueError('nonfinite field')
        return np.full(expected_count,v)
    m=re.search(r'\binternalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\((.*?)\)\s*;',src,flags=re.S)
    if not m:raise ValueError('Unable to parse ASCII scalar internalField')
    count=int(m.group(1))
    if count!=expected_count:raise ValueError(f'Cell count {count} != geometry count {expected_count}')
    vals=np.fromstring(m.group(2),sep=' ')
    if vals.size!=count or not np.isfinite(vals).all():raise ValueError('Invalid scalar values/count')
    return vals


def remap_from_temperature(cell_csv,temperature_foam_field,out, nr=24,nz=40,
                           default_temperature_K=850.,source_coil_current_A_m2=2e6):
    C,V=load_actual_cells(cell_csv)
    Ts_cells=read_foam_internal_scalar(temperature_foam_field,len(V))
    if np.any(Ts_cells<200) or np.any(Ts_cells>2500):raise ValueError('Unsafe/unsupported solid temperature range')
    grid=Grid(nr=nr,nz=nz)
    z,r=grid.coords()
    cr=np.hypot(C[:,0],C[:,1])
    # Nearest neighbor extrapolation only inside annular bed; outside remains reference.
    interp=griddata(np.column_stack([C[:,2],cr]),Ts_cells,np.column_stack([z.ravel(),r.ravel()]),method='nearest').reshape(nz,nr)
    bed=(r>=.005)&(r<=.025)&(z>=.040)&(z<=.100)
    em_temp=np.where(bed,interp,default_temperature_K)
    em_sum,em_fields=em_solve(grid,bed_temperature_K=em_temp,coil_current_A_m2=source_coil_current_A_m2)
    reference_power=em_ref_power(em_fields,.005,.025,.040,.100,nr=90,nz=120)
    q,report=map_conservatively(em_fields,C,V,reference_power)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    write_foam_scalar(out/'qEM',q)
    report.update({'T_cells_min_K':float(Ts_cells.min()),'T_cells_max_K':float(Ts_cells.max()),
                   'electric_conductivity_T_feedback':True,
                   'input_geometry_source':str(cell_csv),
                   'input_temperature_field':str(temperature_foam_field),
                   'em_nz_nr':[nz,nr],
                   'note':'Conservative file exchange; OpenFOAM feedback-loop runtime is not verified'})
    (out/'feedback_report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--cell-csv',required=True)
    p.add_argument('--ts',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--nr',type=int,default=24)
    p.add_argument('--nz',type=int,default=40)
    a=p.parse_args()
    print(json.dumps(remap_from_temperature(a.cell_csv,a.ts,a.out,a.nr,a.nz),indent=2))
