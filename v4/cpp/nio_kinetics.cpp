// Literature rate-law demonstration, NOT validated volumetric reactor kinetics.
// Non-isothermal NiO reduction kinetics: Chemical Engineering Science 63 (2008)
// 567-575, DOI 10.1016/j.ces.2007.09.043.
// Ea=90.8 kJ/mol, ln(A[min^-1])=19.50; f(alpha)=alpha^0.63 (1-alpha)^1.39.
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
static double dadt(double T_K, double alpha) {
    const double R=8.314462618, Ea=90800.0;
    if (!(T_K>0 && alpha>=0 && alpha<=1)) throw std::runtime_error("invalid input");
    const double per_second=std::exp(19.50-Ea/(R*T_K))/60.0;
    return per_second*std::pow(alpha,.63)*std::pow(1-alpha,1.39);
}
int main(int argc,char**argv){
  try {
    const double temp=argc>1?std::stod(argv[1]):800.0;
    const double alpha=argc>2?std::stod(argv[2]):0.2;
    std::cout<<std::setprecision(16)<<dadt(temp,alpha)<<"\n";
    return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 2;}
}
