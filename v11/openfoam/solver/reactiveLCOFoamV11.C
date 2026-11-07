/*
 * Research candidate for OpenFOAM.com v2606; compilation and runtime NOT VERIFIED.
 * Incompressible gas, fixed porous solids, illustrative LiCoO2/C reaction network.
 * This is NOT a plasma MHD solver and NOT an industrial reactor digital twin.
 * All species mol/m3; gas concentrations per pore volume, solid per bulk volume.
 */
#include "fvCFD.H"
#include "pimpleControl.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    pimpleControl pimple(mesh);

    #define READ_FLD(F) volScalarField F(IOobject(#F,runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh)
    READ_FLD(p);
    READ_FLD(CO); READ_FLD(CO2); READ_FLD(N2);
    READ_FLD(LCO); READ_FLD(Carbon); READ_FLD(Li2O); READ_FLD(CoO);
    READ_FLD(CoMetal); READ_FLD(Li2CO3);
    READ_FLD(Tg); READ_FLD(Ts); READ_FLD(qEM);
    #undef READ_FLD
    volVectorField U(IOobject("U",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh);
    surfaceScalarField phi(IOobject("phi",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),fvc::flux(U));

    IOdictionary kin(IOobject("kineticsV11",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE));
    const scalar A1=readScalar(kin.lookup("A1_per_min"));
    const scalar A2=readScalar(kin.lookup("A2_per_min"));
    const scalar Ea1=readScalar(kin.lookup("Ea1_J_mol"));
    const scalar Ea2=readScalar(kin.lookup("Ea2_J_mol"));
    const scalar k3=readScalar(kin.lookup("k3_s"));
    const scalar porosity=0.45;
    const dimensionedScalar eps("eps",dimless,porosity);
    const dimensionedScalar invEps("invEps",dimless,1.0/porosity);
    const dimensionedScalar nu("nu",dimensionSet(0,2,-1,0,0,0,0),4e-5);
    const dimensionedScalar D("D",dimensionSet(0,2,-1,0,0,0,0),1.2e-5);
    const dimensionedScalar alphaGas("alphaGas",dimensionSet(0,2,-1,0,0,0,0),5e-5);
    const dimensionedScalar alphaSolid("alphaSolid",dimensionSet(0,2,-1,0,0,0,0),1e-6);
    // Constant volumetric capacities are reduced-model placeholders.
    const dimensionedScalar Cs("Cs",dimensionSet(1,-1,-2,-1,0,0,0),1.05e6);
    const dimensionedScalar CgPore("CgPore",dimensionSet(1,-1,-2,-1,0,0,0),340.);
    const dimensionedScalar hA("hA",dimensionSet(1,-1,-3,-1,0,0,0),48.);
    const dimensionedScalar r0("r0",dimensionSet(0,-3,-1,0,1,0,0),0.);
    const dimensionedScalar drag0("drag0",dimensionSet(0,0,-1,0,0,0,0),0.);
    const dimensionedScalar q0("q0",dimensionSet(1,-1,-3,0,0,0,0),0.);
    volScalarField ra(IOobject("ra",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),mesh,r0);
    volScalarField rb(IOobject("rb",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),mesh,r0);
    volScalarField rc(IOobject("rc",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),mesh,r0);
    volScalarField drag(IOobject("drag",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,drag0);
    volScalarField qRad(IOobject("qRad",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),mesh,q0);
    const dimensionedScalar gammaGas=hA/(eps*CgPore);
    const dimensionedScalar gammaSolid=hA/Cs;

    while(runTime.run())
    {
        ++runTime;
        Info<< "Time = " << runTime.timeName() << nl;
        while(pimple.loop())
        {
            const scalar delta=runTime.deltaTValue();
            forAll(ra,celli)
            {
                const scalar T=Ts[celli];
                if(T<973.15 || T>1373.15)
                {
                    FatalErrorInFunction << "Temperature out of 700-1100C kinetic domain: "
                        << T << exit(FatalError);
                }
                const scalar k1=A1/60.0*Foam::exp(-Ea1/(8.314462618*T));
                const scalar k2=A2/60.0*Foam::exp(-Ea2/(8.314462618*T));
                scalar a=Foam::min(0.5*k1*Foam::max(0.0,LCO[celli]),
                    Foam::min(Foam::max(0.0,LCO[celli])/2,Foam::max(0.0,Carbon[celli]))/delta);
                scalar b=Foam::min(0.5*k2*Foam::max(0.0,CoO[celli]),
                    Foam::min(Foam::max(0.0,CoO[celli])/2,Foam::max(0.0,Carbon[celli])/2)/delta);
                const scalar carbonRate=a+2*b;
                const scalar scale=Foam::min(1.0,Foam::max(0.0,Carbon[celli])/
                    (delta*Foam::max(carbonRate,SMALL)));
                a*=scale; b*=scale;
                const scalar yCO2=Foam::min(1.0,Foam::max(0.0,
                    CO2[celli]*8.314462618*Foam::max(Tg[celli],1.0)/101325.0));
                const scalar c=Foam::min(k3*Foam::max(0.0,Li2O[celli])*yCO2,
                    Foam::min(Foam::max(0.0,Li2O[celli]),
                              porosity*Foam::max(0.0,CO2[celli]))/delta);
                ra[celli]=a; rb[celli]=b; rc[celli]=c;
                const scalar dp=0.003;
                drag[celli]=150*Foam::pow(1-porosity,2)*nu.value()/
                    (Foam::pow(porosity,3)*dp*dp)
                    +1.75*(1-porosity)*Foam::mag(U[celli])/
                    (Foam::pow(porosity,3)*dp);
                const scalar wallAreaOverV=2*.025/(.025*.025-.005*.005)*.15;
                qRad[celli]=.70*5.670374419e-8*wallAreaOverV*
                    (Foam::pow(Foam::max(T,1.0),4)-Foam::pow(1010.,4));
            }
            ra.correctBoundaryConditions();rb.correctBoundaryConditions();rc.correctBoundaryConditions();
            drag.correctBoundaryConditions();qRad.correctBoundaryConditions();

            fvVectorMatrix UEqn
            (
                fvm::ddt(U)+fvm::div(phi,U)-fvm::laplacian(nu,U)+fvm::Sp(drag,U)
            );
            UEqn.relax();
            if(pimple.momentumPredictor())solve(UEqn == -fvc::grad(p));
            volScalarField rAU("rAU",1.0/UEqn.A());
            volVectorField HbyA("HbyA",rAU*UEqn.H());
            surfaceScalarField phiHbyA("phiHbyA",fvc::flux(HbyA));
            while(pimple.correct())
            {
                while(pimple.correctNonOrthogonal())
                {
                    fvScalarMatrix pEqn(fvm::laplacian(rAU,p)==fvc::div(phiHbyA));
                    pEqn.solve();
                    if(pimple.finalNonOrthogonalIter())phi=phiHbyA-pEqn.flux();
                }
                U=HbyA-rAU*fvc::grad(p);
                U.correctBoundaryConditions();
            }
            // A: 2LiCoO2 + C -> Li2O + 2CoO + CO
            // B: 2CoO + 2C -> 2Co + 2CO
            // C: Li2O + CO2 -> Li2CO3
            solve(fvm::ddt(CO)+invEps*fvm::div(phi,CO)-fvm::laplacian(D,CO)
                  ==(ra+2*rb)/eps);
            solve(fvm::ddt(CO2)+invEps*fvm::div(phi,CO2)-fvm::laplacian(D,CO2)
                  ==-rc/eps);
            solve(fvm::ddt(N2)+invEps*fvm::div(phi,N2)-fvm::laplacian(D,N2)
                  ==r0/eps);
            solve(fvm::ddt(LCO)==-2*ra);
            solve(fvm::ddt(Carbon)==-ra-2*rb);
            solve(fvm::ddt(Li2O)==ra-rc);
            solve(fvm::ddt(CoO)==2*ra-2*rb);
            solve(fvm::ddt(CoMetal)==2*rb);
            solve(fvm::ddt(Li2CO3)==rc);
            solve(fvm::ddt(Tg)+invEps*fvm::div(phi,Tg)-fvm::laplacian(alphaGas,Tg)
                  ==gammaGas*(Ts-Tg));
            // Reaction enthalpies are UNCALIBRATED illustrative values J/mol extent.
            const dimensionedScalar dH1("dH1",dimensionSet(1,2,-2,0,-1,0,0),45000.);
            const dimensionedScalar dH2("dH2",dimensionSet(1,2,-2,0,-1,0,0),-30000.);
            const dimensionedScalar dH3("dH3",dimensionSet(1,2,-2,0,-1,0,0),-25000.);
            solve(fvm::ddt(Ts)-fvm::laplacian(alphaSolid,Ts)
                  ==gammaSolid*(Tg-Ts)+(qEM-qRad-dH1*ra-dH2*rb-dH3*rc)/Cs);
        }
        const scalar maxDiv=gMax(mag(fvc::div(phi))).value();
        if(gMin(LCO).value()<-1e-8 || gMin(Carbon).value()<-1e-8 ||
           gMin(CO).value()<-1e-8 || gMin(CO2).value()<-1e-8)
        {
            FatalErrorInFunction<< "Negative inventory" << exit(FatalError);
        }
        Info<< "V11 meanLCO=" << LCO.weightedAverage(mesh.V()).value()
            << " maxTs=" << gMax(Ts).value()
            << " maxDivPhi=" << maxDiv
            << " EMpower=" << fvc::domainIntegrate(qEM).value()
            << " reactionA=" << fvc::domainIntegrate(ra).value() << nl;
        runTime.write();
    }
    Info<< "End" << endl;
    return 0;
}
