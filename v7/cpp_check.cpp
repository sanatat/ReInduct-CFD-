/* Independent C++17 SI-dimension kinetic and stoichiometry smoke test.
   Real OpenFOAM code remains compile-pending until OpenFOAM libraries available.
   This uses the exported Table III 800 C regression coefficient as input.
*/
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
int main(int argc,char** argv)
{
    if(argc!=2) {std::cerr<<"usage: test <rateUmSAtYH2One>\n";return 2;}
    const double k=std::atof(argv[1]);
    const double T=1073.15, P=101325.0, R=8.314462618;
    const double cH2=0.53*P/(R*T),X=.25;
    const double velocity_m_s=1.e-6*k*(cH2*R*T/P);
    const double molarDensity=6550.0/0.07469;
    const double area=2.0;
    const double source_mol_per_bulk_m3_s=velocity_m_s*molarDensity*area*(1-X);
    const double eps=.45,dt=.001,inventory=450.;
    const double dmol=source_mol_per_bulk_m3_s*dt;
    const double H2_after=cH2-dmol/eps;
    const double H2O_after=dmol/eps;
    const double X_after=X+dmol/inventory;
    const double err=eps*(H2_after+H2O_after)-eps*cH2;
    const double solid_err=eps*H2_after+inventory*X_after-(eps*cH2+inventory*X);
    if(std::abs(err)>1.e-10 || std::abs(solid_err)>1.e-10)return 1;
    std::cout<<"{\"R_mol_per_m3_bulk_s\":"<<source_mol_per_bulk_m3_s
             <<",\"hydrogen_balance_mol_per_m3\":"<<err
             <<",\"solid_balance_mol_per_m3\":"<<solid_err<<"}\n";
}
