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

//Info<< "Reading transport properties\n" << endl;
// read viscosity and permeability

Info<< "Reading polymer properties\n" << endl;
multicomponentPolymer polymer(mesh);

volScalarField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    polymer.diffusion(rho, T)
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
    polymer.pressure(rho, T)
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

