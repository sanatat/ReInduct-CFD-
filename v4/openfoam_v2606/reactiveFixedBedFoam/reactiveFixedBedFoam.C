/*
 * PENDING-COMPILE OpenFOAM.com v2606 proof-of-concept port for a stationary gas/solid
 * continuum. 2-D incompressible PISO flow + H2/H2O + two temperatures + X.
 * The fixed-bed chemistry is representative NiO + H2 -> Ni + H2O.
 * qEM is an offline FVM cell field exported by v5/em2d.py, NOT a
 * fully coupled EM-PISO solver, plasma, or moving particles.
 * Validate units, numerical stability and OpenFOAM API after first compilation.
 */
#include "fvCFD.H"
#include "pisoControl.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading U p H2 H2O Tg Ts X qEM" << nl;
    volVectorField U(IOobject("U",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh);
    volScalarField p(IOobject("p",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh);
    volScalarField H2(IOobject("H2",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh);
    volScalarField H2O(IOobject("H2O",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh);
    volScalarField Tg(IOobject("Tg",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh);
    volScalarField Ts(IOobject("Ts",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh);
    volScalarField X(IOobject("X",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh);
    volScalarField qEM(IOobject("qEM",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
    surfaceScalarField phi(IOobject("phi",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),fvc::flux(U));

    const dimensionedScalar nu("nu",dimensionSet(0,2,-1,0,0,0,0),5e-5);
    const dimensionedScalar D("D",dimensionSet(0,2,-1,0,0,0,0),1e-5);
    const dimensionedScalar aGas("aGas",dimensionSet(0,2,-1,0,0,0,0),0.008/(0.45*450.));
    const dimensionedScalar aSolid("aSolid",dimensionSet(0,2,-1,0,0,0,0),1./1.188e6);
    const dimensionedScalar gammaGas("gammaGas",dimensionSet(0,0,-1,0,0,0,0),600./(0.45*450.));
    const dimensionedScalar gammaSolid("gammaSolid",dimensionSet(0,0,-1,0,0,0,0),600./1.188e6);
    // qEM [W/m3] from Python EM field; rhoCpS [J/(m3 K)] for physical K/s.
    const dimensionedScalar rhoCpS("rhoCpS",dimensionSet(1,-1,-2,-1,0,0,0),1.188e6);
    const dimensionedScalar rxnHeat("rxnHeat",dimensionSet(0,0,0,1,0,0,0),2200.*150./1.188e6); // K per unit-X/s
    const scalar gasConsumption = 150./(0.45*15.);
    volScalarField porousDrag(IOobject("porousDrag",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,
        dimensionedScalar("zeroDrag",dimensionSet(0,0,-1,0,0,0,0),0.));
    volScalarField Rvol(IOobject("Rvol",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,
        dimensionedScalar("zeroRate",dimensionSet(0,0,-1,0,0,0,0),0.));

    Info<< "qEM max: " << gMax(qEM).value() << " W/m3" << nl;
    pisoControl piso(mesh);
    while (runTime.run())
    {
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl;
        forAll(porousDrag, celli)
        {
            const scalar x=mesh.C()[celli].x();
            porousDrag[celli] = (x>=0.036 && x<=0.084)
                ? 0.8*max(0.2,1.0+0.008*(Ts[celli]-850.))*(1.0-0.3*X[celli]) : 0.;
            const scalar T=max(Ts[celli],300.);
            const scalar k=1500.0*exp(-85000.0/(8.314462618*T));
            const scalar possible=k*max(H2[celli],0.)*max(1.-X[celli],0.);
            const scalar available=max(H2[celli],0.)/max(runTime.deltaTValue()*gasConsumption,SMALL);
            const scalar solidAvailable=max(1.-X[celli],0.)/max(runTime.deltaTValue(),SMALL);
            Rvol[celli]=min(possible,min(available,solidAvailable));
        }
        porousDrag.correctBoundaryConditions();
        Rvol.correctBoundaryConditions();

        fvVectorMatrix UEqn
        (
            fvm::ddt(U) + fvm::div(phi,U) - fvm::laplacian(nu,U)
            + fvm::Sp(porousDrag,U)
        );
        // Projection predictor/corrector; pressure is kinematic [m2/s2].
        solve(UEqn == -fvc::grad(p));
        volScalarField rAU("rAU",1.0/UEqn.A());
        U = rAU*UEqn.H();
        U.correctBoundaryConditions();
        phi=fvc::flux(U);
        while (piso.correct())
        {
            fvScalarMatrix pEqn(fvm::laplacian(rAU,p)==fvc::div(phi));
            pEqn.setReference(0,0.);
            pEqn.solve();
            phi-=pEqn.flux();
        }
        U-=rAU*fvc::grad(p);
        U.correctBoundaryConditions();
        solve(fvm::ddt(H2)+fvm::div(phi,H2)-fvm::laplacian(D,H2)==-gasConsumption*Rvol);
        solve(fvm::ddt(H2O)+fvm::div(phi,H2O)-fvm::laplacian(D,H2O)==gasConsumption*Rvol);
        solve(fvm::ddt(X)==Rvol);
        solve(fvm::ddt(Tg)+fvm::div(phi,Tg)-fvm::laplacian(aGas,Tg)
              ==gammaGas*(Ts-Tg));
        solve(fvm::ddt(Ts)-fvm::laplacian(aSolid,Ts)
              ==gammaSolid*(Tg-Ts)+qEM/rhoCpS+rxnHeat*Rvol);
        runTime.write();
        Info<< "mean(X): " << X.weightedAverage(mesh.V()).value()
            << " min(H2): " << min(H2).value() << endl;
    }
    Info<< "End" << endl;
    return 0;
}
