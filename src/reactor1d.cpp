// ReInduct-CFD: reduced-order 1-D reactive porous-bed benchmark.
// Copyright (c) 2026. Educational prototype; parameters illustrative.
// C++17; no third-party dependencies.
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
constexpr double pi=3.141592653589793;
constexpr double R=8.314462618;

struct Config {
    int n=60;
    double dt=0.01, end=120.0, q=600000.0;
    double Ea=85000.0, A=1500.0;
    double H2in=15.0, U=0.016, Tin=850.0;
    std::string out="results/base";
};

Config arguments(int argc, char** argv) {
    Config c;
    for (int i=1;i<argc;i++) {
        const std::string s=argv[i];
        if (s=="--help") {
            std::cout << "reactor1d [--nx INT] [--dt SEC] [--tend SEC] "
                         "[--heat W_PER_M3] [--Ea J_PER_MOL] [--A PER_SEC] "
                         "[--h2-in MOL_PER_M3] [--velocity M_PER_S] "
                         "[--inlet-temp K] [--out DIR]\n";
            std::exit(0);
        }
        if (i+1>=argc) throw std::invalid_argument("Missing value after "+s);
        const std::string v=argv[++i];
        if(s=="--nx") c.n=std::stoi(v);
        else if(s=="--dt") c.dt=std::stod(v);
        else if(s=="--tend") c.end=std::stod(v);
        else if(s=="--heat") c.q=std::stod(v);
        else if(s=="--Ea") c.Ea=std::stod(v);
        else if(s=="--A") c.A=std::stod(v);
        else if(s=="--h2-in") c.H2in=std::stod(v);
        else if(s=="--velocity") c.U=std::stod(v);
        else if(s=="--inlet-temp") c.Tin=std::stod(v);
        else if(s=="--out") c.out=v;
        else throw std::invalid_argument("Unknown argument: "+s);
    }
    if(c.n<4 || c.n>100000 || !(c.dt>0) || !(c.end>0) || c.end>1e6 ||
       c.q<0 || c.H2in<0 || c.Ea<0 || c.A<0 || c.U<=0 ||
       c.Tin<=0)
        throw std::invalid_argument("Unphysical or out-of-range input");
    return c;
}

// One-dimensional surrogate for the species and energy balances of a fixed solid bed.
// Gas: superficial advection + interphase reaction and heat exchange.
// Solid: thermal diffusion + assumed volumetric induction-like heat input.
// This is NOT an OpenFOAM CFD calculation or a Maxwell/MHD/plasma solution.
int main(int argc, char** argv) {
 try {
    const Config p=arguments(argc,argv);
    const int N=p.n;
    const double L=0.12;               // reactor length [m]
    const double eps=0.45;             // gas volume fraction [-]
    const double U=p.U;                // superficial gas velocity [m/s]
    const double nNiO=150.0;           // initial active NiO [mol/m3 bulk]
    const double CpG=450.0;           // gas volumetric heat capacity [J/m3 gas/K]
    const double CpS=1.188e6;         // solid volumetric heat capacity [J/m3 bulk/K]
    const double hVol=800.0;          // interphase heat transfer [W/m3 bulk/K]
    const double kw=1.0;              // effective axial bed conductivity [W/m/K]
    const double hLoss=35.0;          // distributed wall heat loss [W/m3 bulk/K]
    const double Tin=p.Tin, Tamb=850.0;
    const double dH=-2500.0;          // assumed reaction enthalpy [J/mol reaction]
    const double cRef=15.0;           // reference gas concentration [mol/m3 gas]
    const double dz=L/N;
    const double cfl=U*p.dt/(eps*dz);
    const double gasExchange=p.dt*hVol/(eps*CpG);
    const double solidDiff=p.dt*kw/(CpS*dz*dz);
    if(cfl>0.9 || gasExchange>0.4 || 2*solidDiff>0.4)
        throw std::invalid_argument("dt too large: decrease dt (CFL/heat stability)");
    const int steps=static_cast<int>(std::ceil(p.end/p.dt));
    const double dt=p.end/steps;
    fs::create_directories(p.out);
    std::ofstream field(p.out+"/fields.csv"), series(p.out+"/history.csv");
    if(!field || !series) throw std::runtime_error("Cannot open output files");
    field << "z_m,Tgas_K,Tsolid_K,H2_mol_m3gas,H2O_mol_m3gas,solid_conversion\n";
    series << "time_s,Tsolid_mean_K,Tsolid_max_K,H2_out_mol_m3gas,"
              "H2O_out_mol_m3gas,conversion_mean,H2_balance_abs_mol_m2,"
              "H2_balance_relative\n";

    std::vector<double> h2(N,0.0), h2o(N,0.0), Tg(N,Tin), Ts(N,Tin), X(N,0.0);
    std::vector<double> nh2(N),nh2o(N),nTg(N),nTs(N),nX(N);
    double totalIn=0.0,totalOut=0.0, lastResidual=0.0;
    double maxTs=0.0,meanX=0.0,meanTs=0.0;
    auto report=[&](double time){
        meanTs=0.0;meanX=0.0;maxTs=0.0;
        double inventory=0.0,reduced=0.0;
        for(int i=0;i<N;++i){
            meanTs+=Ts[i]/N;
            meanX+=X[i]/N;
            maxTs=std::max(maxTs,Ts[i]);
            inventory+=eps*h2[i]*dz;
            reduced+=nNiO*X[i]*dz;
        }
        // Conservative accounting: inlet H2 = outlet H2 + gas inventory + reduced NiO.
        lastResidual=totalIn-totalOut-inventory-reduced;
        double rel=std::abs(lastResidual)/std::max(totalIn,1.0);
        series << std::setprecision(12) << time << ',' << meanTs << ',' << maxTs
               << ',' << h2.back() << ',' << h2o.back() << ',' << meanX
               << ',' << std::abs(lastResidual) << ',' << rel << '\n';
    };
    report(0.0);
    for(int k=0;k<steps;++k) {
        // Explicit, upwind face fluxes for a positive superficial gas velocity.
        totalIn+=dt*U*p.H2in;
        totalOut+=dt*U*h2.back();
        for(int i=0;i<N;++i) {
            const double h2up=(i==0?p.H2in:h2[i-1]);
            const double h2oup=(i==0?0.0:h2o[i-1]);
            const double Tup=(i==0?Tin:Tg[i-1]);
            const double h2Adv=std::max(0.0,h2[i]+cfl*(h2up-h2[i]));
            const double h2oAdv=std::max(0.0,h2o[i]+cfl*(h2oup-h2o[i]));
            double TgAdv=Tg[i]+cfl*(Tup-Tg[i]);
            // An illustrative irreversible apparent first-order conversion model.
            // It is deliberately NOT fitted to literature kinetics or radical chemistry.
            const double kv=p.A*std::exp(-p.Ea/(R*Ts[i]));
            const double rawRate=nNiO*kv*(h2Adv/cRef)*(1.0-X[i]);
            const double rate=std::max(0.0,std::min({rawRate,
                     eps*h2Adv/dt,nNiO*(1.0-X[i])/dt}));
            nh2[i]=h2Adv-dt*rate/eps;
            nh2o[i]=h2oAdv+dt*rate/eps;
            nX[i]=X[i]+dt*rate/nNiO;
            const double qInter=hVol*(Ts[i]-TgAdv);
            nTg[i]=TgAdv+dt*qInter/(eps*CpG);
            const double lTs=(i==0?Ts[i]:Ts[i-1]);
            const double rTs=(i==N-1?Ts[i]:Ts[i+1]);
            const double diffusion=kw*(lTs-2.0*Ts[i]+rTs)/(dz*dz);
            const double z=(i+0.5)*dz;
            const double localization=std::exp(-std::pow((z-0.065)/0.04,2.0));
            const double qInd=p.q*localization;   // assumed heating source, NOT EM-solved
            nTs[i]=Ts[i]+dt*(diffusion+qInd-qInter-hLoss*(Ts[i]-Tamb)
                             -dH*rate)/CpS;
        }
        h2.swap(nh2);h2o.swap(nh2o);Tg.swap(nTg);Ts.swap(nTs);X.swap(nX);
        if((k+1)%std::max(1,steps/120)==0 || k+1==steps) report((k+1)*dt);
    }
    for(int i=0;i<N;++i){
        field << std::setprecision(12) << (i+0.5)*dz << ',' << Tg[i] << ','
              << Ts[i] << ',' << h2[i] << ',' << h2o[i] << ',' << X[i] << '\n';
    }
    // Calculate an analytic electromagnetic skin depth as a separate reference.
    // Valid under the ideal non-magnetic, uniform-conductor approximation.
    const double mu0=4.0*pi*1e-7;
    const double frequency=100000.0;     // Hz
    const double sigma=1e6;             // S/m; illustrative conducting solid
    const double delta=std::sqrt(2.0/(2.0*pi*frequency*mu0*sigma));
    std::ofstream summary(p.out+"/summary.txt");
    summary << std::setprecision(10)
            << "ReInduct-CFD 1D educational surrogate (NOT OpenFOAM)\n"
            << "cells=" << N << " dt_s=" << dt << " t_end_s=" << p.end << '\n'
            << "U_m_s=" << U << " Tin_K=" << Tin << " heat_peak_W_m3=" << p.q << '\n'
            << "mean_conversion=" << meanX << '\n'
            << "mean_solid_temperature_K=" << meanTs << '\n'
            << "maximum_solid_temperature_K=" << maxTs << '\n'
            << "max_absolute_H2_balance_error_mol_m2=" << std::abs(lastResidual) << '\n'
            << "H2_balance_error_definition=inflow-outflow-gas_inventory-reduced_solid\n"
            << "analytical_skin_depth_m=" << delta << '\n'
            << "skin_depth_assumptions=mu_r=1,f=100kHz,sigma=1e6S/m\n"
            << "WARNING: Heating is prescribed; not computed from Maxwell equations.\n"
            << "WARNING: Rates and parameters are illustrative, unvalidated.\n";
    summary.flush();
    std::cout << std::setprecision(7) << "DONE, N=" << N << " conversion=" << meanX
              << " Tmax[K]=" << maxTs << " |H2 balance error|="
              << std::abs(lastResidual) << " mol/m2; output=" << p.out << '\n';
    return std::abs(lastResidual)>1e-8 ? 2 : 0;
 } catch(const std::exception& e){ std::cerr << "ERROR: " << e.what() << '\n'; return 1; }
}
