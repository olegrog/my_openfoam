Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
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

Info<< "Reading transport properties\n" << endl;
singlePhaseTransportModel laminarTransport(U, phi);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

if (p.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

Info<< "Reading boundary conditions\n" << endl;
IOdictionary boundaryConditions
(
    IOobject
    (
        "boundaryConditions",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    boundaryConditions
);
T.dimensions().reset(dimTemperature);

volScalarField phase
(
    IOobject
    (
        "phase",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    boundaryConditions
);

volScalarField nGrain
(
    IOobject
    (
        "nGrain",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    boundaryConditions
);

volScalarField Co
(
    IOobject
    (
        "Co",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    boundaryConditions
);

boundaryConditions.set("internalField", "uniform (0 0 0)");

volVectorField grain
(
    IOobject
    (
        "grain",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    boundaryConditions
);

Info<< "Reading phase-field properties\n" << endl;
IOdictionary phaseFieldProperties
(
    IOobject
    (
        "phaseFieldProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar interfaceWidth("interfaceWidth", phaseFieldProperties);

dictionary crystallographicAngles = phaseFieldProperties.subDict("crystallographicAngles");

Info<< "Reading problem-specific properties\n" << endl;
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

dimensionedScalar undercooling("undercooling", problemProperties);
dimensionedScalar tempGradient("tempGradient", problemProperties);
dimensionedScalar coolingRate("coolingRate", problemProperties);

Info<< "Reading alloy properties\n" << endl;
multicomponentAlloy alloy(mesh);
