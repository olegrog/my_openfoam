bool writeAllFields(runTime.controlDict().get<bool>("writeAllFields"));

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"
volScalarField divU("divU", fvc::div(phi)); // Used in correctPhi.H

Info<< "Reading transport and thermal properties\n" << endl;
incompressibleGasMetalMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());
const volScalarField& T = mixture.T();

const volScalarField& rho1 = mixture.rhoM();
const dimensionedScalar& rho2 = mixture.rho2();

// Uncomment this to analyze the difference in implementations of |\grad\alpha|
/*
{
    Info<< "Write magGradAlpha fields\n" << endl;
    word name = "magGradAlpha";
    volScalarField& alpha = mixture.alpha1();
    std::vector<volScalarField> magGradAlphaFields
    {
        volScalarField(name + "1", mag(fvc::grad(alpha))),
        volScalarField(name + "2", mag(fvc::reconstruct(fvc::snGrad(alpha)*mesh.magSf()))),
        volScalarField(name + "3", 3*mag(fvc::average(fvc::snGrad(alpha)))),
        volScalarField(name + "4", mag(fvc::reconstructMag(fvc::snGrad(alpha)*mesh.magSf()))),
    };
    for (const auto& f : magGradAlphaFields) f.write();
}
*/

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

volScalarField radiativeCooling
(
    IOobject
    (
        "radiativeCooling",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeAllFields ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimPower/dimArea, Zero)
);

volScalarField evaporativeCooling
(
    IOobject
    (
        "evaporativeCooling",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeAllFields ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimPower/dimArea, Zero)
);

volScalarField heatConvection
(
    IOobject
    (
        "heatConvection",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeAllFields ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimPower/dimVolume, Zero)
);

volScalarField heatConduction
(
    IOobject
    (
        "heatConduction",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeAllFields ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimPower/dimVolume, Zero)
);

volScalarField heatConduction2
(
    IOobject
    (
        "heatConduction2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeAllFields ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimPower/dimVolume, Zero)
);

volScalarField energySource
(
    IOobject
    (
        "energySource",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        writeAllFields ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimPower/dimMass, Zero)
);


Info<< "Reading problemProperties\n" << endl;
IOdictionary problemProperties
(
    IOobject
    (
        "problemProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
dimensionedScalar ambientPressure("ambientPressure", problemProperties);
dimensionedScalar ambientTemperature("ambientTemperature", problemProperties);


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        writeAllFields ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);

// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    mixture
);

// Density for buoyancy force
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
    rho*(1.0 - beta*(T - ambientTemperature))
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rhok*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(mixture.phase1Name());

isoAdvection advector(mixture.alpha1(), phi, U);

// NB: MRF should be registered before creating laserHeatSource
autoPtr<movingReferenceFrame> MRF(movingReferenceFrame::New(mesh));

Info<< "Reading laser properties\n" << endl;
autoPtr<surfaceLaserHeatSource> laserHeatSource
(
    surfaceLaserHeatSource::New(mixture, advector)
);
Info<< endl;
