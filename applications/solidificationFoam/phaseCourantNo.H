scalar maxPhaseCo(runTime.controlDict().get<scalar>("maxPhaseCo"));

// --- Advective Co imported from CourantNo.H
volScalarField meshV = volScalarField::New
(
    "cellVolume",
    mesh,
    dimVolume,
    mesh.V(),
    zeroGradientFvPatchField<scalar>::typeName
);
meshV.correctBoundaryConditions();
Co = 0.5*fvc::surfaceSum(mag(phi))/meshV*runTime.deltaT();

// --- Estimate the interface velocity via rhs
volScalarField rhsPhase = mag(sqr(interfaceWidth)*fvc::laplacian(phase) - fPrime(phase)
    + a1*interfaceWidth/alloy.interfaceEnergy()
        *gPrime(phase)*alloy.chemicalDrivingForce(phase, T))/tau;
scalar phaseCoNum = gMax(rhsPhase)*runTime.deltaTValue();
Info<< "Interface Courant number max = " << phaseCoNum << endl;
Co = max(Co, rhsPhase*runTime.deltaT());

// --- Estimate evolution of concentration via rhs
forAllIters(alloy.components(), iter)
{
    alloyComponent& C = iter();
    volScalarField h = alloy.partition(phase);
    volVectorField normal = calcNormal(fvc::grad(phase));
    volScalarField rhsC = mag
    (
        fvc::laplacian(C.diffusion(phase), (C - C.equilibrium(phase, T))/h)
      + fvc::div(interfaceWidth/Foam::sqrt(2.)*C.deltaA()*normal*rhsPhase)
    )/C.equilibrium(phase, T);
    scalar CoNumC = gMax(rhsC)*runTime.deltaTValue();
    Info<< C.name() << " Courant number max = " << CoNumC << endl;
    phaseCoNum = max(phaseCoNum, CoNumC);
    Co = max(Co, rhsC*runTime.deltaT());
}

