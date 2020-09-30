Info<< "Reading field T\n" << endl;
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

Info<< "Reading field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading polymer properties\n" << endl;
multicomponentPolymer polymer(mesh);

volScalarField D1
(
    IOobject
    (
        "D1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimViscosity
);

volScalarField D2
(
    IOobject
    (
        "D2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimViscosity,
    zeroGradientFvPatchField<scalar>::typeName
);

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
    mesh,
    dimPressure
);

Info<< "Reading transport properties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar mu("mu", dimDynamicViscosity, transportProperties);
dimensionedScalar meanParticleSize("meanParticleSize", dimLength, transportProperties);
dimensionedScalar permeability0 = sqr(meanParticleSize)/180;

volScalarField permeability
(
    IOobject
    (
        "permeability",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimArea
);

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

dimensionedScalar heatingRate("heatingRate", dimTemperature/dimTime, problemProperties);

