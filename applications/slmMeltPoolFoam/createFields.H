#include "Read/thermalProperties.H"
#include "Read/laserProperties.H"

bool printProperties(runTime.controlDict().get<bool>("printProperties"));

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
    mesh,
    dimless
);
calcLiquidFraction(liquidFraction, T, T_solidus, T_liquidus);

volScalarField wasMelted
(
    IOobject
    (
        "wasMelted",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    liquidFraction
);

Info << "Reading field he\n" << endl;
volScalarField he
(
    IOobject
    (
        "he",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimEnergy/dimMass,
    T.primitiveField(),
    T.boundaryField()
);
calcEnthalpy(he, T, liquidFraction, Cp_sol, Cp_liq, dCp_sol, dCp_liq, T_solidus, T_liquidus, enthalpyFusion);

volScalarField dLiquidFraction
(
    IOobject
    (
        "dLiquidFraction",
        runTime.timeName(),
        mesh
    ),
    mesh,
    inv(he.dimensions())
);
calcLiquidFractionDer(dLiquidFraction, T, T_solidus, T_liquidus);

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
dimensionedScalar permeability("permeability", mixture);

autoPtr<Function1<scalar>> pSigma
(
    Function1<scalar>::New("sigma", mixture.subDict("sigma"))
);
dimensionedScalar dsigma_dT
(
    "dsigma_dT",
    mixture.sigmaK()().dimensions() / dimTemperature * dimLength,
    pSigma->value(1) - pSigma->value(0)
);
Info<< "Surface tension = " << pSigma->value(0) << " + (" << dsigma_dT.value() << ")*T" << endl;

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

volScalarField heatCapacity
(
    IOobject
    (
        "heatCapacity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        printProperties ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimSpecificHeatCapacity
);

volScalarField thermalConductivity
(
    IOobject
    (
        "thermalConductivity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        printProperties ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimPower/dimLength/dimTemperature
);

volScalarField laserHeatSource
(
    IOobject
    (
        "laserHeatSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        printProperties ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimPower/dimVolume
);

// Used only for writing
volScalarField kinematicViscosity
(
    IOobject
    (
        "kinematicViscosity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        printProperties ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mixture.nu()
);

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
        IOobject::READ_IF_PRESENT,
        printProperties ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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

isoAdvection advector(alpha1, phi, U);

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
dimensionedScalar molarMass("molarMass", problemProperties);
dimensionedScalar ambientPressure("ambientPressure", problemProperties);

Switch segregation(
    pimple.dict().lookupOrDefault<Switch>("segregation", false)
);
autoPtr<multicomponentAlloy> pAlloy;
if (segregation) {
    // This field is necessary for creating multicomponentAlloy
    volScalarField phase
    (
        IOobject
        (
            "phase",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh, dimless, zeroGradientFvPatchField<scalar>::typeName
    );

    Info<< "Reading alloyProperties\n" << endl;
    pAlloy = autoPtr<multicomponentAlloy>::New(mesh);
}

