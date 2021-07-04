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
const volScalarField rhsPhase = mag(sqr(interfaceWidth)*fvc::laplacian(phase) - fPrime(phase)
    + a1*interfaceWidth/alloy.interfaceEnergy()
        *gPrime(phase)*alloy.chemicalDrivingForce(phase, T))/tau;

scalar phaseCoNum = gMax(rhsPhase)*runTime.deltaTValue();
Info<< "Interface Courant number max = " << phaseCoNum << endl;
Co = max(Co, rhsPhase*runTime.deltaT());

// --- Estimate evolution of concentration via rhs
const volScalarField k = alloy.partition(phase, T);
const volVectorField nHat = calcNormal(fvc::grad(phase));

for (const auto& C : alloy.components())
{
    const volScalarField Ceq = C.equilibrium(phase, T);
    const volScalarField rhsC =
        mag(fvc::div(aAT*interfaceWidth*C.delta(T)*nHat*rhsPhase))/Ceq*k;

    const scalar CoNumC = gMax(rhsC)*runTime.deltaTValue();
    Info<< C.name() << " Courant number max = " << CoNumC << endl;
    phaseCoNum = max(phaseCoNum, CoNumC);
    Co = max(Co, rhsC*runTime.deltaT());
}

