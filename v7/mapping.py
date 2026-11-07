"""Conservative q_EM mapping from v6 axisymmetric (r,z) to OpenFOAM cells.

The input cell geometry must come from the actual mesh utility (NOT assumed
blockMesh ordering). Energy conservation refers to the chosen CFD subdomain;
conserving power does not imply local interpolation error is small.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from v6.axisymmetric_em import Grid,run as solve_em
ROOT=Path(__file__).resolve().parent

def em_ref_power(em_fields,r0,r1,z0,z1,angle_rad=np.deg2rad(4.),nr=160,nz=220):
    rc=r0+(np.arange(nr)+.5)*(r1-r0)/nr
    zc=z0+(np.arange(nz)+.5)*(z1-z0)/nz
    interp=RegularGridInterpolator((em_fields['z'],em_fields['r']),em_fields['q'],bounds_error=False,fill_value=0.)
    zz,rr=np.meshgrid(zc,rc,indexing='ij')
    q=interp(np.stack([zz.ravel(),rr.ravel()],axis=1)).reshape(nz,nr)
    # azimuthal sector volume wedge: r dtheta dr dz
    return float(np.sum(q*rr)*angle_rad*((r1-r0)/nr)*((z1-z0)/nz))

def synthetic_annular_wedge(nr=12,nz=48,r0=.005,r1=.025,z0=.040,z1=.100,angle=np.deg2rad(4.)):
    """For OFFLINE V&V only. Does not claim OpenFOAM cell indexing."""
    re=np.linspace(r0,r1,nr+1);ze=np.linspace(z0,z1,nz+1)
    centres=[];vol=[]
    for j in range(nz):
        for i in range(nr):
            # exact area-weighted centroid radius of cylindrical annular sector
            a,b=re[i],re[i+1]
            rc=(2/3)*(b**3-a**3)/(b*b-a*a)
            centres.append([rc,0.,.5*(ze[j]+ze[j+1])])
            vol.append(.5*(b*b-a*a)*angle*(ze[j+1]-ze[j]))
    return np.asarray(centres),np.asarray(vol)

def load_actual_cells(csv_path):
    with open(csv_path) as f:rows=list(csv.DictReader(f))
    if not rows or not {'x_m','y_m','z_m','volume_m3'}.issubset(rows[0]):
        raise ValueError('Require x_m,y_m,z_m,volume_m3 columns from exportCellGeometry')
    xyz=np.array([[float(row[k]) for k in ('x_m','y_m','z_m')] for row in rows]);V=np.array([float(row['volume_m3']) for row in rows])
    if not np.all(np.isfinite(xyz)) or not np.all(np.isfinite(V)) or np.any(V<=0):
        raise ValueError('Invalid cell geometry')
    return xyz,V

def map_conservatively(em_fields,centres,V,reference_power,positive_only=True):
    interp=RegularGridInterpolator((em_fields['z'],em_fields['r']),em_fields['q'],bounds_error=False,fill_value=0.)
    r=np.hypot(centres[:,0],centres[:,1]);z=centres[:,2]
    q_raw=interp(np.column_stack([z,r]));q_raw=np.maximum(0,q_raw)
    raw_power=float(q_raw@V)
    if raw_power<=0 or reference_power<=0:raise ValueError('Power source must be positive')
    scale=reference_power/raw_power
    q_out=q_raw*scale
    result={'cell_count':len(V),'reference_power_W':float(reference_power),'raw_mapped_power_W':raw_power,
            'raw_relative_power_error':float(raw_power/reference_power-1),
            'scale':float(scale),'final_mapped_power_W':float(q_out@V),
            'post_mapping_relative_power_error':float((q_out@V)/reference_power-1),
            'min_q_W_m3':float(q_out.min()),'max_q_W_m3':float(q_out.max())}
    return q_out,result

def write_foam_scalar(path,values):
    values=np.asarray(values).reshape(-1)
    if not np.isfinite(values).all():raise ValueError('nonfinite q')
    boundary='\n'.join(f'{key} {{ type {val}; }}' for key,val in [('inlet','zeroGradient'),('outlet','zeroGradient'),('innerWall','zeroGradient'),('outerWall','zeroGradient'),('front','wedge'),('back','wedge')])
    header='''FoamFile
{
 version 2.0;
 format ascii;
 class volScalarField;
 location "0";
 object qEM;
}
dimensions [1 -1 -3 0 0 0 0];
'''
    body='\n'.join(f'{v:.15g}' for v in values)
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(header+f'internalField nonuniform List<scalar>\n{len(values)}\n(\n{body}\n);\n'+f'boundaryField {{\n{boundary}\n}}\n')

def demo(out=ROOT/'results',cell_csv=None,grid_n=40):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    report,fields=solve_em(Grid(nr=grid_n,nz=round(grid_n*1.6)))
    if cell_csv:
        C,V=load_actual_cells(cell_csv);src='OpenFOAM-generated actual mesh geometry'
    else:
        C,V=synthetic_annular_wedge();src='synthetic annular geometry (offline, NOT OpenFOAM)'
    pref=em_ref_power(fields,.005,.025,.040,.100)
    mapped,metrics=map_conservatively(fields,C,V,pref)
    metrics.update({'source_cells':src,'em_n_r':grid_n,'em_power_entire_domain_W':report['joule_power_W'],
                    'mapping_domain':'annular reactor: r=0.005..0.025 m; z=0.040..0.100 m; wedge=4 degrees',
                    'claim':'conservative offline field mapping; not OpenFOAM run'})
    write_foam_scalar(out/'qEM',mapped)
    (out/'power_mapping.json').write_text(json.dumps(metrics,indent=2)+'\n')
    return metrics

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cell-csv',default=None);p.add_argument('--out',default='v7/results/em_mapping');p.add_argument('--grid-n',type=int,default=40)
    args=p.parse_args();print(json.dumps(demo(args.out,args.cell_csv,args.grid_n),indent=2))
