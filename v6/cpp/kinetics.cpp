// Standalone C++17 volumetric surrogate for NiO + H2 -> Ni + H2O.
// Units: T [K], gas concentrations [mol/m3_gas], Kp dimensionless,
// output mol/(m3_bed s). This is NOT a packed-bed experimentally fitted law.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>

static double reaction_rate(const double T, const double cH2,
                           const double cH2O, const double X,
                           const double Kp)
{
    if (!(T >= 500. && T <= 1000. && cH2 >= 0. && cH2O >= 0.
        && X >= 0. && X <= 1. && Kp > 0.))
        throw std::invalid_argument("invalid thermodynamic range or concentrations");
    constexpr double R=8.314462618;
    constexpr double epsilon=0.45;
    constexpr double k0=750.0;
    constexpr double Ea=48000.0;
    const double driving=std::max(cH2-cH2O/Kp,0.0);
    return epsilon*k0*std::exp(-Ea/(R*T))*driving*(1.0-X);
}
int main(int argc,char**argv){
    if(argc!=6){std::cerr<<"usage: kinetics T_K cH2 cH2O X Kp\n";return 2;}
    try{
        const double v=reaction_rate(std::stod(argv[1]),std::stod(argv[2]),
                           std::stod(argv[3]),std::stod(argv[4]),std::stod(argv[5]));
        std::cout<<std::setprecision(17)<<v<<"\n";
    } catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 3;}
}
