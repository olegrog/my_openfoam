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

    Info<< "Reading field phase\n" << endl;
    volScalarField phase
    (
        IOobject
        (
            "phase",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
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
        dimensionedScalar("1", dimless, 0)
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

    dimensionedScalar interfaceWidth(phaseFieldProperties.lookup("interfaceWidth"));

    scalar a1(phaseFieldProperties.lookupOrDefault<scalar>("a1", Foam::sqrt(2.)/6));
    scalar a2(phaseFieldProperties.lookupOrDefault<scalar>("a2", 2.35));
    scalar epsilon4(phaseFieldProperties.lookupOrDefault<scalar>("epsilon4", 0));

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

    dimensionedScalar undercooling(problemProperties.lookup("undercooling"));
    dimensionedScalar tempGradient(problemProperties.lookup("tempGradient"));
    dimensionedScalar coolingRate(problemProperties.lookup("coolingRate"));

    scalar frontPositionRel(problemProperties.lookupOrDefault<scalar>("frontPosition", 0));
    scalar interfaceNarrowing(problemProperties.lookupOrDefault<scalar>("interfaceNarrowing", 1));
    scalar seedNarrowing(problemProperties.lookupOrDefault<scalar>("seedNarrowing", 1));
    label nSeeds(problemProperties.lookupOrDefault<label>("nSeeds", 0));

