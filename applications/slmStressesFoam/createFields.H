Info<< "Reading transport and thermal properties\n" << endl;
quiescentGasMetalMixture mixture(mesh);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

// Melt indicator (the maximum liquid fraction in time)
volScalarField wasMelted
(
    IOobject
    (
        "wasMelted",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar()
);

// Need to store rho for ddt(rho, h)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);

Info<< "Reading laser properties\n" << endl;
autoPtr<volumetricLaserHeatSource> laserHeatSource
(
    volumetricLaserHeatSource::New(mesh)
);
Info<< endl;
