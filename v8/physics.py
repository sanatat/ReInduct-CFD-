"""Auditable physical correlations for the v8 *illustrative* reactor.

No correlation here is claimed to be calibrated for industrial battery black mass.
All concentrations are mol/m^3 of void gas; porosity is void fraction; q in W/m^3 bulk.
"""
from __future__ import annotations
from dataclasses import dataclass
import math

R=8.314462618
SIGMA_SB=5.670374419e-8
MW_NI=0.0586934    # kg / mol of metallic nickel
RHO_NI=8908.0      # kg / m^3 of dense metallic nickel
MW_NIO=0.0746928  # kg / mol NiO
RHO_NIO=6550.0     # kg / m^3 dense NiO (used for bed inventory only)


def layer_speed_to_bulk_rate(speed_um_s:float, area_m2_per_m3_bulk:float,
                             ni_layer_porosity:float=0.0, unconverted_fraction:float=1.0)->float:
    """Convert measured Ni product-LAYER growth speed to NiO molar conversion rate.

    Nickel-layer growth thickness v [m/s] * Ni atoms per m^3 of nickel product
    *(1-void fraction)* active interface area per m^3 of bulk bed = mol/(m^3_bulk s).
    The NiO/Ni reaction is 1:1. Geometry and product-layer porosity are unknown,
    so inferred bulk rate is illustrative (not validated by NiO slab kinetics).
    """
    if not (0<=ni_layer_porosity<1 and 0<=unconverted_fraction<=1):
        raise ValueError('porosity and remaining conversion outside physical bounds')
    if speed_um_s<0 or area_m2_per_m3_bulk<0 or not math.isfinite(speed_um_s+area_m2_per_m3_bulk):
        raise ValueError('nonnegative finite input required')
    return speed_um_s*1e-6*(RHO_NI/MW_NI)*(1-ni_layer_porosity)*area_m2_per_m3_bulk*unconverted_fraction


def molar_density_ratio()->float:
    return (RHO_NI/MW_NI)/(RHO_NIO/MW_NIO)


def h2_mole_fraction(c_h2_mol_m3:float, Tgas_K:float, pressure_pa:float)->float:
    if c_h2_mol_m3<0 or Tgas_K<=0 or pressure_pa<=0:
        raise ValueError('nonphysical gas inputs')
    y=c_h2_mol_m3*R*Tgas_K/pressure_pa
    if y>1+1e-8:
        raise ValueError('H2 concentration exceeds ideal-gas total concentration')
    return min(1.,y)


def ideal_gas_mixture(h2_mole_fraction_val:float, h2o_mole_fraction_val:float,
                      T_K:float, pressure_Pa:float)->dict:
    """Explicit H2/H2O/N2 mixture: ideal-gas rho and illustrative Cp(T).
    Cp polynomials are approximate engineering placeholders, not NASA-fit correlations.
    """
    yh, yw=h2_mole_fraction_val,h2o_mole_fraction_val
    if not (T_K>0 and pressure_Pa>0 and yh>=0 and yw>=0 and yh+yw<=1+1e-12):
        raise ValueError('invalid mole fractions / state')
    yn=max(0.,1-yh-yw)
    molar_mass=yh*.00201588+yw*.01801528+yn*.0280134
    # Baseline heat capacities 300 K and weak illustrative T-dependence; J/(mol K).
    cp_molar=yh*(28.8+.002*(T_K-300))+yw*(33.6+.010*(T_K-300))+yn*(29.1+.003*(T_K-300))
    rho=pressure_Pa*molar_mass/(R*T_K)
    return {'yH2':yh,'yH2O':yw,'yN2':yn,'rho_kg_m3':rho,'cp_molar_J_mol_K':cp_molar,
            'cp_mass_J_kg_K':cp_molar/molar_mass,
            'volumetric_heat_capacity_J_m3_K':pressure_Pa*cp_molar/(R*T_K)}


def volumetric_radiation_loss(Tsolid_K:float,Tamb_K:float,emissivity:float,
                              exposed_wall_area_per_bulk_volume_m_inv:float)->float:
    """Grey-wall diffuse emission loss [W/m^3 bulk]; >0 when solid is hotter.
    A/V must be the *EXPOSED wall* area, not the total particle surface area.
    """
    if min(Tsolid_K,Tamb_K)<0 or not 0<=emissivity<=1 or exposed_wall_area_per_bulk_volume_m_inv<0:
        raise ValueError('invalid radiation parameters')
    return emissivity*SIGMA_SB*exposed_wall_area_per_bulk_volume_m_inv*(Tsolid_K**4-Tamb_K**4)


def ergun_pressure_gradient(superficial_U_m_s:float,epsilon:float,particle_diameter_m:float,
                            mu_Pa_s:float,rho_kg_m3:float)->float:
    if not (0<epsilon<1 and particle_diameter_m>0 and mu_Pa_s>0 and rho_kg_m3>0):
        raise ValueError('invalid Ergun inputs')
    return (150*(1-epsilon)**2*mu_Pa_s*superficial_U_m_s/(epsilon**3*particle_diameter_m**2)
            +1.75*(1-epsilon)*rho_kg_m3*abs(superficial_U_m_s)*superficial_U_m_s
                /(epsilon**3*particle_diameter_m))
