// Independent C++17 SI unit/stoichiometry smoke test for Ni product-layer conversion.
// NOT a compiled OpenFOAM solver and not a black-mass experimental rate law.
#include <cmath>
#include <iostream>
int main() {
    constexpr double k_um_per_s=2.744874796;
    constexpr double R=8.314462618, P=101325., Tgas=850.;
    constexpr double NiMolarDensity=8908./0.0586934;
    constexpr double yH2=.53, area=2., X=.25, eps=.45, site=450.,dt=.001;
    const double cH2=yH2*P/(R*Tgas);
    const double yComputed=cH2*R*Tgas/P;
    const double molarRate=k_um_per_s*1e-6*yComputed*NiMolarDensity*area*(1-X);
    const double delta=molarRate*dt;
    const double h2=cH2-delta/eps,h2o=delta/eps,newX=X+delta/site;
    const double hBalance=eps*(h2+h2o-cH2);
    const double solidBalance=eps*h2+site*newX-(eps*cH2+site*X);
    if(!std::isfinite(molarRate)||molarRate<=0||std::abs(hBalance)>1e-11||std::abs(solidBalance)>1e-11) return 1;
    std::cout << "{\"rate_mol_m3_bulk_s\":"<<molarRate
              <<",\"molar_density_Ni_layer_mol_m3\":"<<NiMolarDensity
              <<",\"mole_fraction_consistency\":"<<yComputed
              <<",\"H_balance\":"<<hBalance
              <<",\"solid_balance\":"<<solidBalance<<"}\n";
    return 0;
}
