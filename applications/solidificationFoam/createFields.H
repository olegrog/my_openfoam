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

    volVectorField smallVector
    (
        IOobject("smallVector", runTime.constant(), mesh),
        mesh,
        1e-20 * dimensionedVector("1", dimless / dimLength, vector(0, 1, 0))
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

    scalar a1(phaseFieldProperties.lookupOrDefault<scalar>("a1", Foam::sqrt(2.)/6));
    scalar a2(phaseFieldProperties.lookupOrDefault<scalar>("a2", 2.35));
    scalar epsilon4(phaseFieldProperties.lookupOrDefault<scalar>("epsilon4", 0));
    scalar frontPosition(phaseFieldProperties.lookupOrDefault<scalar>("frontPosition", 0));

    dimensionedScalar interfaceWidth(phaseFieldProperties.lookup("interfaceWidth"));
    dimensionedScalar undercooling(phaseFieldProperties.lookup("undercooling"));
    dimensionedScalar tempGradient(phaseFieldProperties.lookup("tempGradient"));
    dimensionedScalar coolingRate(phaseFieldProperties.lookup("coolingRate"));

    multicomponentAlloy alloy(mesh);
    dimensionedScalar tau = a1 * a2 * pow3(interfaceWidth) * alloy.relaxationTime();

    volScalarField coordX = mesh.C().component(vector::X);
    volScalarField coordY = mesh.C().component(vector::Y);
