// Standalone C++17 consistency test for four v12 reaction pathways.
// Does not compile or validate the OpenFOAM solver.
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
int main() {
 constexpr int ns=11, ne=6, nr=4;
 constexpr std::array<std::array<double,ne>,ns> a{{
  {{1,1,0,2,0,0}},{{0,0,1,0,0,0}},{{2,0,0,1,0,0}},
  {{0,1,0,1,0,0}},{{0,1,0,0,0,0}},{{2,0,1,3,0,0}},
  {{0,0,1,1,0,0}},{{0,0,1,2,0,0}},{{0,0,0,0,2,0}},
  {{0,0,0,0,0,2}},{{0,0,0,0,0,1}}
 }};
 constexpr std::array<std::array<double,nr>,ns> s{{
  {{-2,0,0,0}},{{-1,-2,0,0}},{{1,0,-1,0}},{{2,-2,0,0}},
  {{0,2,0,0}},{{0,0,1,0}},{{1,2,0,0}},{{0,0,-1,0}},
  {{0,0,0,0}},{{0,0,0,-1}},{{0,0,0,2}}
 }};
 double maxResidual=0;
 for(int e=0;e<ne;++e)for(int r=0;r<nr;++r){
  double v=0;for(int i=0;i<ns;++i)v+=a[i][e]*s[i][r];
  maxResidual=std::max(maxResidual,std::abs(v));
 }
 if(maxResidual>1e-12)throw std::runtime_error("element balance failed");
 std::cout<<"PASS C++17 11-species, 4-pathway elemental conservation: max="<<maxResidual<<"\n";
}
