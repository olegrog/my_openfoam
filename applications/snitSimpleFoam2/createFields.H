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
    T = T0;

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
    U = U1;

    volScalarField curvature
    (
        IOobject
        (
            "curvature",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
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

    IOdictionary knudsenLayerProperties
    (
        IOobject
        (
            "knudsenLayerProperties",
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
    dimensionedScalar gamma1
    (
        transportProperties.lookup("gamma1")
    );

    Info<< "Reading gamma2\n" << endl;
    dimensionedScalar gamma2
    (
        transportProperties.lookup("gamma2")
    );

    Info<< "Reading gamma7\n" << endl;
    dimensionedScalar gamma7
    (
        transportProperties.lookup("gamma7")
    );
