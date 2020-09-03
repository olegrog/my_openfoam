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

Info<< "Reading field Ul\n" << endl;
volVectorField Ul
(
    IOobject
    (
        "Ul",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Ul
);

#include "createPhi.H"

Info<< "Reading transport properties\n" << endl;
singlePhaseTransportModel laminarTransport(Ul, phi);

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
volScalarField boundaryConditions
(
    IOobject
    (
        "boundaryConditions",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
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
    dimTemperature,
    boundaryConditions.primitiveField(),
    boundaryConditions.boundaryField()
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
    dimless,
    boundaryConditions.primitiveField(),
    boundaryConditions.boundaryField()
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
    dimless,
    boundaryConditions.primitiveField(),
    boundaryConditions.boundaryField()
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
    dimless,
    boundaryConditions.primitiveField(),
    boundaryConditions.boundaryField()
);

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
    dimensionedVector(),
    boundaryConditions.boundaryField().types()
);

Info<< "Reading alloy properties\n" << endl;
multicomponentAlloy alloy(mesh);

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

dimensionedScalar undercooling("undercooling", dimTemperature, problemProperties);
dimensionedScalar tempGradient("tempGradient", dimTemperature/dimLength, problemProperties);
dimensionedScalar coolingRate("coolingRate", dimTemperature/dimTime, problemProperties);
dimensionedScalar interfaceWidth("interfaceWidth", dimLength, problemProperties);
dictionary crystallographicAngles = problemProperties.subDict("crystallographicAngles");

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

