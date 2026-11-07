/*
 * OpenFOAM.com v2606 candidate. COMPILE/RUNTIME PENDING external installation.
 * Constant-density incompressible stationary bed, axisymmetric wedge geometry.
 * IMPORTANT: No claim of a validated battery-recycling reactor.
 *
 * H2,H2O [mol/m^3 of pore gas], X [-], qEM [W/m^3 of bulk bed].
 * eps*c balances, solid NiO capacity [mol/m^3 of bulk bed].
 * Constant gas volumetric Cp, two-temperature LTE departure, Ergun drag.
 */
#include "fvCFD.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    volVectorField U(IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    volScalarField p(IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    volScalarField H2(IOobject("H2", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    volScalarField H2O(IOobject("H2O", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    volScalarField Tg(IOobject("Tg", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    volScalarField Ts(IOobject("Ts", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    volScalarField X(IOobject("X", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    volScalarField qEM(IOobject("qEM", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    surfaceScalarField phi(IOobject("phi",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),fvc::flux(U));

    // SI dimensions: nu, diffusivity, thermal diffusivities [m2/s].
    const dimensionedScalar nu("nu", dimensionSet(0,2,-1,0,0,0,0),3.2e-5);
    const dimensionedScalar D("D", dimensionSet(0,2,-1,0,0,0,0),1.2e-5);
    const dimensionedScalar alphaGas("alphaGas", dimensionSet(0,2,-1,0,0,0,0),5.0e-5);
    const dimensionedScalar alphaSolid("alphaSolid", dimensionSet(0,2,-1,0,0,0,0),1.0e-6);
    const dimensionedScalar porosity("eps",dimless,0.45);
    const scalar eps=0.45;
    // The solid inventory is per bulk m3, not per m3 of solids.
    const dimensionedScalar solidNiO("solidNiO",dimensionSet(0,-3,0,0,1,0,0),450.0);
    // Effective Cp_solid per bulk volume [J/(m3*K)] and Cp_gas pore volume [J/(m3*K)].
    const dimensionedScalar Csolid("Csolid",dimensionSet(1,-1,-2,-1,0,0,0),1.05e6);
    const dimensionedScalar CgasPore("CgasPore",dimensionSet(1,-1,-2,-1,0,0,0),800.);
    const dimensionedScalar heatTransfer("heatTransfer",dimensionSet(1,-1,-3,-1,0,0,0),600.);
    const dimensionedScalar deltaH("deltaH",dimensionSet(1,2,-2,0,-1,0,0),-2200.);
    // Grey-body radiation is represented by a *lumped wall-area/bed-volume* sink.
    // Geometric approximation for annular cylindrical outer wall, not P1/DOM radiation.
    const scalar radiativeWallAreaPerBulkVolume=2*.025/(.025*.025-.005*.005); // 1/m
    const scalar emissivity=.70, wallAmbientK=700.;
    const dimensionedScalar radZero("radZero",dimensionSet(1,-1,-3,0,0,0,0),0.0);
    volScalarField qRad(IOobject("qRad",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),mesh,radZero);
    const dimensionedScalar rateZero("rateZero",dimensionSet(0,-3,-1,0,1,0,0),0.0);
    const dimensionedScalar dragZero("dragZero",dimensionSet(0,0,-1,0,0,0,0),0.0);
    volScalarField reactionRate(IOobject("reactionRate",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),mesh,rateZero);
    volScalarField porousDrag(IOobject("porousDrag",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dragZero);
    const dimensionedScalar cGas("Cg",dimless,eps);
    const dimensionedScalar invEPS("invEPS",dimless,1.0/eps);
    // hA/C [s-1]; heat of reaction -deltaH*R (exothermic if deltaH < 0).
    const dimensionedScalar gammaGas=heatTransfer/(porosity*CgasPore);
    const dimensionedScalar gammaSolid=heatTransfer/Csolid;

    IOdictionary kinetics
    (
        IOobject("kineticTable",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE)
    );
    scalarList tempC(kinetics.lookup("temperatureC"));
    scalarList rateUm(kinetics.lookup("rateUmPerSAtYH2One"));
    const scalar aSurface=readScalar(kinetics.lookup("reactiveAreaPerBulkVolume"));
    const scalar rhoNi_mol_m3=readScalar(kinetics.lookup("NiProductLayerMolarDensity"));
    const scalar niLayerPorosity=readScalar(kinetics.lookup("NiProductLayerPorosity"));
    if(niLayerPorosity<0 || niLayerPorosity>=1) FatalErrorInFunction << "Invalid Ni layer porosity" << exit(FatalError);
    const scalar idealGasPressurePa=101325.0;
    if(tempC.size()!=rateUm.size() || tempC.size()<2)
    {
        FatalErrorInFunction << "Invalid kineticTable" << exit(FatalError);
    }
    while (runTime.run())
    {
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl;
        const label maxCoupling=5;
        for (label outer=0;outer<maxCoupling;++outer)
        {
            const scalarField prevH2(H2.primitiveField());
            const scalarField prevTs(Ts.primitiveField());
            const scalar dp=0.003; // representative particle diameter [m]
            forAll(reactionRate,i)
            {
                const scalar T=max(Ts[i],300.0);
                // Table-III *Ni-layer growth speed* [um/s] from dense NiO samples.
                // Translate to porous-bed volumetric rate using ASSUMED exposed area.
                // This is an illustrative cross-scale mapping, not a validated kinetic law.
                const scalar tc=T-273.15;
                // Never silently extrapolate the data-based kinetic law.
                if(tc<tempC[0] || tc>tempC[tempC.size()-1])
                {
                    FatalErrorInFunction << "T out of 500-1000C kinetic range: " << tc
                        << exit(FatalError);
                }
                label j=0;
                while(j+2<tempC.size() && tc>tempC[j+1]) ++j;
                const scalar weight=(tc-tempC[j])/(tempC[j+1]-tempC[j]);
                const scalar kUmS=(1-weight)*rateUm[j]+weight*rateUm[j+1];
                const scalar yH2=min(1.,max(0.,H2[i]*8.314462618*max(Tg[i],1.0)/idealGasPressurePa));
                const scalar speed_m_s=1.e-6*kUmS*yH2;
                scalar rate=speed_m_s*rhoNi_mol_m3*(1-niLayerPorosity)*aSurface*max(0.,1.-X[i]);
                // Explicit bound via local inventories per time step.
                const scalar availableGas=eps*max(0.0,H2[i])/max(runTime.deltaTValue(),SMALL);
                const scalar availableSolid=450.0*max(0.0,1.0-X[i])/max(runTime.deltaTValue(),SMALL);
                reactionRate[i]=min(rate,min(availableGas,availableSolid));
                qRad[i]=emissivity*5.670374419e-8*radiativeWallAreaPerBulkVolume
                        *(pow(max(Ts[i],1.0),4)-pow(wallAmbientK,4));
                // Ergun resistance coefficient [1/s], U is superficial m/s.
                const scalar linear=150*pow(1-eps,2)*nu.value()/(pow(eps,3)*dp*dp);
                const scalar quadratic=1.75*(1-eps)*mag(U[i])/(pow(eps,3)*dp);
                porousDrag[i]=linear+quadratic;
            }
            reactionRate.correctBoundaryConditions();
            porousDrag.correctBoundaryConditions();
            qRad.correctBoundaryConditions();
            fvVectorMatrix UEqn
            (
                fvm::ddt(U) + fvm::div(phi,U) - fvm::laplacian(nu,U)
                + fvm::Sp(porousDrag,U)
            );
            // PISO-like tentative velocity and single non-orthogonal pressure projection.
            // Candidate API: requires compilation against OpenFOAM.com v2606.
            solve(UEqn==-fvc::grad(p));
            volScalarField rAU("rAU",1.0/UEqn.A());
            volVectorField HbyA("HbyA",rAU*UEqn.H());
            surfaceScalarField phiHbyA("phiHbyA",fvc::flux(HbyA));
            fvScalarMatrix pEqn(fvm::laplacian(rAU,p)==fvc::div(phiHbyA));
            pEqn.solve();
            phi=phiHbyA-pEqn.flux();
            U=HbyA-rAU*fvc::grad(p);
            U.correctBoundaryConditions();
            // 0th-order stationary solid. Gas species concentrations refer to pore gas.
            solve(fvm::ddt(H2)+invEPS*fvm::div(phi,H2)
                -fvm::laplacian(D,H2)==-reactionRate/porosity);
            solve(fvm::ddt(H2O)+invEPS*fvm::div(phi,H2O)
                -fvm::laplacian(D,H2O)==reactionRate/porosity);
            solve(fvm::ddt(X)==reactionRate/solidNiO);
            solve(fvm::ddt(Tg)+invEPS*fvm::div(phi,Tg)
                -fvm::laplacian(alphaGas,Tg)==gammaGas*(Ts-Tg));
            solve(fvm::ddt(Ts)-fvm::laplacian(alphaSolid,Ts)
                ==gammaSolid*(Tg-Ts)+(qEM-qRad-deltaH*reactionRate)/Csolid);
            scalar maxChange=0.;
            forAll(H2,i)
            {
                maxChange=max(maxChange,mag(H2[i]-prevH2[i])/max(1.0,mag(H2[i])));
                maxChange=max(maxChange,mag(Ts[i]-prevTs[i])/max(1.0,mag(Ts[i])));
            }
            Info<< "coupling " << outer+1 << " relativeChange " << maxChange << nl;
            if(maxChange<1e-7)break;
        }
        // Reject out-of-range states instead of quietly clipping concentrations.
        if(gMin(H2).value() < -1e-7 || gMin(H2O).value() < -1e-7
           || gMin(X).value() < -1e-7 || gMax(X).value()>1.+1e-7)
        {
            FatalErrorInFunction << "Negative species or out-of-range X" << exit(FatalError);
        }
        Info<< "meanX=" << X.weightedAverage(mesh.V()).value()
            << " minH2=" << gMin(H2).value()
            << " maxTs=" << gMax(Ts).value()
            << " totalEMPower=" << fvc::domainIntegrate(qEM).value()
            << " totalRadiativeLoss=" << fvc::domainIntegrate(qRad).value()
            << " maxDivergencePhi=" << gMax(mag(fvc::div(phi))).value() << endl;
        runTime.write();
    }
    Info<< "End" << endl;
    return 0;
}
