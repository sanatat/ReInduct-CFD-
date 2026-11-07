// Standalone C++17 cross-check of the NIST Shomate sensible Cp and h fits.
// NOT an OpenFOAM solver or plasma thermochemistry model.
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using A7 = std::array<double,7>;
struct Fit { double breakpoint; A7 lo,hi; };
const std::array<Fit,4> Fits{{
  {1300,{25.56759,6.096130,4.054656,-2.671301,.131021,-118.0089,-110.5271},{35.15070,1.300095,-.205921,.013550,-3.282780,-127.8375,-110.5271}},
  {1200,{24.99735,55.18696,-33.69137,7.948387,-.136638,-403.6075,-393.5224},{58.16639,2.720074,-.492289,.038844,-6.447293,-425.9186,-393.5224}},
  {500,{28.98641,1.853978,-9.647459,16.63537,.000117,-8.671914,0},{19.50583,19.88705,-8.598535,1.369784,.527601,-4.935202,0}},
  {1000,{33.066178,-11.363417,11.432816,-2.772874,-.158558,-9.980797,0},{18.563083,12.257357,-2.859786,.268238,1.977990,-1.147438,0}}
}};
double raw_h(double T,const A7 &v) {
  const double z=T/1000.0;
  return 1000.0*(v[0]*z+v[1]*z*z/2+v[2]*std::pow(z,3)/3+v[3]*std::pow(z,4)/4-v[4]/z+v[5]-v[6]);
}
double cp(int idx,double T) {
  if(T<298.15 || T>2000.0)throw std::invalid_argument("Temperature outside Shomate range");
  if(idx==4)return 2.5*8.31446261815324;
  auto f=Fits.at(idx);auto a=(T<f.breakpoint?f.lo:f.hi);double z=T/1000.0;
  return a[0]+a[1]*z+a[2]*z*z+a[3]*z*z*z+a[4]/(z*z);
}
double h(int idx,double T) {
  if(T<298.15 || T>2000.0)throw std::invalid_argument("Temperature outside Shomate range");
  if(idx==4)return 2.5*8.31446261815324*(T-298.15);
  auto f=Fits.at(idx);
  return raw_h(T,T<f.breakpoint?f.lo:f.hi)
       +(T>=f.breakpoint?raw_h(f.breakpoint,f.lo)-raw_h(f.breakpoint,f.hi):0.0);
}
int main() {
  const double T=1080;
  std::cout<<std::setprecision(15);
  for(int i=0;i<5;++i) {
    const double derivative=(h(i,T+0.001)-h(i,T-0.001))/0.002;
    if(std::abs(derivative-cp(i,T))>0.01) return 1;
    if(cp(i,T)<19 || cp(i,T)>100)return 2;
    std::cout<<i<<","<<cp(i,T)<<","<<h(i,T)<<"\n";
  }
  return 0;
}
