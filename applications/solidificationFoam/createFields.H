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

    dimensionedScalar sigma(phaseFieldProperties.lookup("sigma"));
    dimensionedScalar interfaceWidth(phaseFieldProperties.lookup("interfaceWidth"));
    dimensionedScalar density(phaseFieldProperties.lookup("density"));
    dimensionedScalar molarMass(phaseFieldProperties.lookup("molarMass"));
    dimensionedScalar molarVolume(molarMass / density);

    multicomponentAlloy alloy(mesh);
    dimensionedScalar tau = a1 * a2 * pow3(interfaceWidth) / sigma * alloy.kineticParameter()
        / molarVolume;

    volScalarField coordX = mesh.C().component(vector::X);
    volScalarField coordY = mesh.C().component(vector::Y);
