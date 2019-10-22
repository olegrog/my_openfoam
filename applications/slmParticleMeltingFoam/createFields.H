#include "Read/thermalProperties.H"
#include "Read/laserProperties.H"
dimensionedScalar he_solidus = (0.5*T_solidus*(2*Cp_sol + dCp_sol*T_solidus));
dimensionedScalar he_liquidus = enthalpyFusion + 0.5*T_solidus*(2*Cp_sol + dCp_sol*T_solidus);

Info << "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField liquidFraction
(
    IOobject
    (
        "liquidFraction",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    tanhSmooth(T, (T_solidus + T_liquidus)/2, (T_liquidus - T_solidus)/2)
);

Info << "Reading field he\n" << endl;
volScalarField he
(
    IOobject
    (
        "he",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    enthalpyCalc(T, liquidFraction, Cp_sol, Cp_liq, dCp_sol, dCp_liq, T_solidus, T_liquidus, enthalpyFusion)
);

volScalarField diffusivity
(
    IOobject
    (
        "diffusivity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1, -1, -1, 0, 0, 0, 0)
);

#include "createRDeltaT.H"

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

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
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
    p_rgh + rho*gh
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
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"

Info<< "Reading problemProperties\n" << endl;
IOdictionary problemProperties(
    IOobject(
        "problemProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
dimensionedScalar permeability(problemProperties.lookup("permeability"));
dimensionedScalar marangoni(problemProperties.lookup("marangoni"));
dimensionedScalar ballRadius(problemProperties.lookup("ballRadius"));
dimensionedScalar molarMass(problemProperties.lookup("molarMass"));
dimensionedScalar ambientPressure(problemProperties.lookup("ambientPressure"));

label nBalls(problemProperties.lookupOrDefault<label>("nBalls", 0));
scalar anisimov(problemProperties.lookupOrDefault<scalar>("anisimov", 0.54));

