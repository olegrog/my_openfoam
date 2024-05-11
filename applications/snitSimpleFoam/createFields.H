    Info<< "Reading field T0\n" << endl;
    volScalarField T0
    (
        IOobject
        (
            "T0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U1\n" << endl;
    volVectorField U1
    (
        IOobject
        (
            "U1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field p0\n" << endl;
    volScalarField p0
    (
        IOobject
        (
            "p0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field p2\n" << endl;
    volScalarField p2
    (
        IOobject
        (
            "p2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U1 / T0) & mesh.Sf()
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p2, simple.dict(), pRefCell, pRefValue);
    dimensionedScalar initialPressure = fvc::domainIntegrate(p2);
    dimensionedScalar totalVolume = sum(mesh.V());

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

    Info<< "Reading transportProperties\n" << endl;

    Info<< "Reading gamma1\n" << endl;
    dimensionedScalar gamma1("gamma1", transportProperties);

    Info<< "Reading gamma2\n" << endl;
    dimensionedScalar gamma2("gamma2", transportProperties);

    Info<< "Reading gamma7\n" << endl;
    dimensionedScalar gamma7("gamma7", transportProperties);

    // the default value correspond to the hard-sphere model s = 0.5
    // Gamma1,2 \propto T^s, \Gamma_3 \propto T^{2s}, \Gamma_7 \propto T^{2s-1}
    scalar s = nondimensionalProperties.lookupOrDefault("s", 0.5);
