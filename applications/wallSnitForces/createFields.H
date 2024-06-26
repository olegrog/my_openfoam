volScalarField T0
(
    IOobject
    (
        "T0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volVectorField U1
(
    IOobject
    (
        "U1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField p0
(
    IOobject
    (
        "p0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField p2
(
    IOobject
    (
        "p2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar g1("gamma1", transportProperties);
dimensionedScalar g2("gamma2", transportProperties);
dimensionedScalar g3("gamma3", transportProperties);
dimensionedScalar g7("gamma7", transportProperties);

IOdictionary nondimensionalProperties
(
    IOobject
    (
        "nondimensionalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

scalar s = nondimensionalProperties.lookupOrDefault("s", 0.5);

