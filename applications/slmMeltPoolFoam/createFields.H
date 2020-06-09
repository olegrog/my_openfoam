#include "Read/laserProperties.H"

bool writeProperties(runTime.controlDict().get<bool>("writeProperties"));

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

Info<< "Reading transport and thermal properties\n" << endl;
gasMetalMixture mixture(U, phi);

volScalarField& rho = mixture.rho();
const volScalarField& T = mixture.T();

const gasMetalThermo& thermo = mixture.thermo();
const dimensionedScalar molarMass = thermo.metalM();
const dimensionedScalar Hfusion = thermo.Hfusion();
const dimensionedScalar Hvapour = thermo.Hvapour();
const dimensionedScalar Tmelting = thermo.Tmelting();
const dimensionedScalar Tboiling = thermo.Tboiling();

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

volVectorField gradT
(
    IOobject
    (
        "gradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeProperties ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    fvc::grad(T)
);

volScalarField laserHeatSource
(
    IOobject
    (
        "laserHeatSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeProperties ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimPower/dimArea
);

volScalarField radiativeCooling
(
    IOobject
    (
        "radiativeCooling",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeProperties ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimPower/dimArea
);

volScalarField evaporativeCooling
(
    IOobject
    (
        "evaporativeCooling",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeProperties ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimPower/dimArea
);

volScalarField heatConvection
(
    IOobject
    (
        "heatConvection",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeProperties ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimPower/dimVolume
);

volScalarField heatConduction
(
    IOobject
    (
        "heatConduction",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeProperties ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimPower/dimVolume
);

volScalarField heatConduction2
(
    IOobject
    (
        "heatConduction2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeProperties ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimPower/dimVolume
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

{
    scalar factor = (molarMass*Hvapour/constant::physicoChemical::R).value();
    Info<< "Vapor pressure = " << ambientPressure.value()
        << "*exp(" << factor/Tboiling.value()
        << " - " << factor << "/T)\n" << endl;
}

volScalarField vaporPressure
(
    IOobject
    (
        "vaporPressure",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimPressure)
);

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

