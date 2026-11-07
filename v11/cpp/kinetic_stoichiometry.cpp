// Portable C++17 cross-check: published stage fits, local elemental closure.
// These k are effective model coefficients; not validated local reaction kinetics.
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <iomanip>
int main() {
    const double R = 8.31446261815324;
    const double T = 1080.;
    const double Ea1 = 121167.4518068571, Ea2=95278.5676391636;
    const double A1 =584.4837776336469, A2=19.794878001298635;
    if(T<973.15||T>1373.15)throw std::out_of_range("published kinetic domain");
    const double k1=A1/60.*std::exp(-Ea1/(R*T));
    const double k2=A2/60.*std::exp(-Ea2/(R*T));
    // Li Co C O N elemental counts for nine species: LCO,C,Li2O,CoO,Co,Li2CO3,CO,CO2,N2
    const std::array<std::array<int,5>,9> atoms{{
      {{1,1,0,2,0}},{{0,0,1,0,0}},{{2,0,0,1,0}},{{0,1,0,1,0}},
      {{0,1,0,0,0}},{{2,0,1,3,0}},{{0,0,1,1,0}},{{0,0,1,2,0}},{{0,0,0,0,2}}
    }};
    // R1 = 2LCO+C->Li2O+2CoO+CO;
    // R2 = 2CoO+2C->2Co+2CO;
    // R3 = Li2O+CO2->Li2CO3.
    const std::array<std::array<int,3>,9> stoich{{
      {{-2,0,0}},{{-1,-2,0}},{{1,0,-1}},{{2,-2,0}},{{0,2,0}},
      {{0,0,1}},{{1,2,0}},{{0,0,-1}},{{0,0,0}}
    }};
    for(int e=0;e<5;++e) for(int j=0;j<3;++j){
      int sum=0;for(int s=0;s<9;++s)sum+=atoms[s][e]*stoich[s][j];
      if(sum!=0)throw std::runtime_error("unbalanced reaction");
    }
    if(!(k1>0 && k2>0))throw std::runtime_error("rates nonpositive");
    std::cout<<std::setprecision(16)<<"{\"temperature_K\":"<<T
        <<",\"k1_per_second\":"<<k1<<",\"k2_per_second\":"<<k2
        <<",\"all_reactions_atom_balanced\":true}\n";
}
