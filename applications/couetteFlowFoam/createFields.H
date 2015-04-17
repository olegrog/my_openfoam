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

    Info<< "Reading field U0\n" << endl;
    volScalarField U0
    (
        IOobject
        (
            "U0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field T1\n" << endl;
    volScalarField T1
    (
        IOobject
        (
            "T1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U1\n" << endl;
    volScalarField U1
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

    Info<< "Reading transportProperties\n" << endl;
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

    Info<< "Reading s\n" << endl;
    dimensionedScalar s
    (
        transportProperties.lookup("s")
    );

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

    Info<< "Reading gamma3\n" << endl;
    dimensionedScalar gamma3
    (
        transportProperties.lookup("gamma3")
    );

    Info<< "Reading gamma7\n" << endl;
    dimensionedScalar gamma7
    (
        transportProperties.lookup("gamma7")
    );

    Info<< "Reading gamma8\n" << endl;
    dimensionedScalar gamma8
    (
        transportProperties.lookup("gamma8")
    );

    Info<< "Reading gamma9\n" << endl;
    dimensionedScalar gamma9
    (
        transportProperties.lookup("gamma9")
    );

    Info<< "Reading gamma10\n" << endl;
    dimensionedScalar gamma10
    (
        transportProperties.lookup("gamma10")
    );

    Info<< "Reading k0\n" << endl;
    dimensionedScalar k0
    (
        transportProperties.lookup("k0")
    );

    Info<< "Reading d1\n" << endl;
    dimensionedScalar d1
    (
        transportProperties.lookup("d1")
    );

    Info<< "Reading kn\n" << endl;
    dimensionedScalar kn
    (
        transportProperties.lookup("kn")
    );

    volSymmTensorField p1
    (
        IOobject
        (
            "p1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField q1
    (
        IOobject
        (
            "q1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volSymmTensorField p2
    (
        IOobject
        (
            "p2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField q2
    (
        IOobject
        (
            "q2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

