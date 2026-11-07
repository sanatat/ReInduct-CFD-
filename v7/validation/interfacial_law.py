"""Cross-scale *Ni product-layer growth* speed -> Ni consumption [mol/(m3_bulk s)].

Uses Ni product-layer effective molar density rather than the NiO reactant density.
Dense-Ni layer assumed unless product-layer porosity is explicitly supplied.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from v7.validation.nio_holdout import read_data
ROOT=Path(__file__).resolve().parents[1]
R=8.314462618

def rate_table_from_all_data():
    data=read_data(); rows=[]
    for t in np.unique(data[:,0]):
        subset=data[data[:,0]==t]
        w=1/subset[:,3]**2
        k=np.sum(w*subset[:,1]*subset[:,2])/np.sum(w*subset[:,1]**2)
        rows.append((t,float(k)))
    return np.asarray(rows)

def ni_surface_speed_um_s(T_C,yH2,table=None):
    table=rate_table_from_all_data() if table is None else table
    if not (table[:,0].min() <= T_C <= table[:,0].max()):
        raise ValueError('Temperature outside published calibration range 500–1000 C')
    return float(np.interp(T_C,table[:,0],table[:,1])*yH2)

def interfacial_to_bulk_rate(T_K,cH2_mol_m3,X,area_m2_m3=2.,rhoNi_kg_m3=8908.,MNi_kg_mol=.0586934,pressure_Pa=101325.,Tgas_K=None,ni_product_porosity=0.):
    """R_bulk = Ni layer speed * Ni product-layer molar density * interface area * (1-X).

    cH2 [mol/m3 pore gas]; H2 mole fraction uses gas temperature Tg,
    not solid reaction temperature Ts. Product-layer porosity is unknown.
    Not an industrial intrinsic rate or validated conversion law.
    """
    if not (0<=X<=1 and 0<=ni_product_porosity<1 and area_m2_m3>=0 and cH2_mol_m3>=0):
        raise ValueError("invalid kinetic or geometrical parameters")
    if Tgas_K is None:Tgas_K=T_K
    from v8.physics import h2_mole_fraction
    y=h2_mole_fraction(cH2_mol_m3,Tgas_K,pressure_Pa)
    rate_m_s=ni_surface_speed_um_s(T_K-273.15,y)*1e-6
    return rate_m_s*(rhoNi_kg_m3/MNi_kg_mol)*(1-ni_product_porosity)*area_m2_m3*(1-X)

def export_foam_table(path=ROOT/'openfoam/case/constant/kineticTable',out=ROOT/'results/interfacial_law.json'):
    table=rate_table_from_all_data(); path=Path(path);path.parent.mkdir(exist_ok=True,parents=True)
    def numlist(v):return '('+' '.join(f'{float(x):.10g}' for x in v)+')'
    path.write_text('''FoamFile
{
 version 2.0;
 format ascii;
 class dictionary;
 object kineticTable;
}
// Table III Hidayat et al. 2009: interpolated pure-H2-equivalent initial Ni layer growth
// Rates [um/s]; temperatures [C]. Must not extrapolate to other temperatures.
'''+f'temperatureC {numlist(table[:,0])};\nrateUmPerSAtYH2One {numlist(table[:,1])};\n'
      + 'reactiveAreaPerBulkVolume 2.0; // [m2/m3bulk], illustrative\n'
      + 'NiProductLayerMolarDensity 151771.75; // [mol/m3 Ni metal layer]; correct for product porosity\n'
      + 'NiProductLayerPorosity 0.0; // illustrative dense Ni layer, not measured\n')
    record={'table_T_C':table[:,0].tolist(),'pure_H2_rate_um_s':table[:,1].tolist(),
        'a_interfacial_m2_m3_bulk_assumed':2.,'rho_Ni_product_kg_m3':8908.,'M_Ni_product_kg_mol':.0586934,
        'scope':'Empirical dense Ni layer-growth / assumed reactive area. Not reactor validation. Full 19-point data for deployed table; temperature holdouts are evaluated separately.'}
    Path(out).parent.mkdir(exist_ok=True,parents=True)
    Path(out).write_text(json.dumps(record,indent=2)+'\n')
    return record

if __name__=='__main__':print(json.dumps(export_foam_table(),indent=2))
