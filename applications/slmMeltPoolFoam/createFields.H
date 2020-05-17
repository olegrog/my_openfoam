#include "Read/thermalProperties.H"
#include "Read/laserProperties.H"

bool writeProperties(runTime.controlDict().get<bool>("writeProperties"));

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
dimensionedScalar mushyCoeff("mushyCoeff", mixture);

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
Info<< "Surface tension = " << pSigma->value(0) << " + (" << dsigma_dT.value() << ")*T\n" << endl;

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

// Uncomment this to analyze the difference in implementations of |\grad\alpha|
/*{
    Info<< "Write magGradAlpha fields\n" << endl;
    word name = "magGradAlpha";
    std::vector<volScalarField> magGradAlphaFields
    {
        volScalarField(name + "1", mag(fvc::grad(alpha1))),
        volScalarField(name + "2", mag(fvc::reconstruct(fvc::snGrad(alpha1) * mesh.magSf()))),
        volScalarField(name + "3", 3*mag(fvc::average(fvc::snGrad(alpha1)))),
        volScalarField(name + "4", mag(fvc::reconstructMag(fvc::snGrad(alpha1) * mesh.magSf()))),
    };
    for (const auto& f : magGradAlphaFields) f.write();
}*/

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
    T.boundaryField()   // copy BC from the temperature field
);

volScalarField heAtFusion(he_melting(alpha2));
volScalarField hePrimeAlpha(he);

LiquidFraction liquidFraction
(
    alpha1, he, heAtFusion, enthalpyFusion, heAtFusionPrime, thermalProperties
);

{
    scalar error, maxError = pimple.dict().get<scalar>("liquidFractionInitError");
    // TODO(olegrog): Use Newton's iterations for boosting
    Info<< "Fixed-point iterations for finding the enthalpy:" << endl;
    do {
        calcEnthalpy(he, T, liquidFraction, alpha2, Cp_sol, dCp_sol, Cp_liq, dCp_liq, Cp_gas,
            T_melting, enthalpyFusion);
        liquidFraction().storePrevIter();
        liquidFraction.correct();
        error = Foam::max(mag(liquidFraction() - liquidFraction().prevIter())).value();
        Info<< "liquidFraction error = " << error << endl;
    } while (error > maxError);
    Info<< endl;
}
// To debug calcTemperature & calcEnthalpy
// he.write();
// calcTemperature(T, he, liquidFraction, alpha2,
//     Cp_sol, dCp_sol, Cp_liq, dCp_liq, Cp_gas, T_melting, enthalpyFusion);
// T.write();

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

volScalarField heatCapacity
(
    IOobject
    (
        "heatCapacity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeProperties ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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
        writeProperties ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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

// Used only for writing
volScalarField kinematicViscosity
(
    IOobject
    (
        "kinematicViscosity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeProperties ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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
        writeProperties ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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
    rho*(1.0 - beta*(T - T_ambient))
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
mesh.setFluxRequired(alpha1.name());

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

{
    scalar factor = (molarMass * enthalpyBoiling / constant::physicoChemical::R).value();
    Info<< "Vapor pressure = " << ambientPressure.value() << " * exp(" << factor/T_boiling.value()
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

