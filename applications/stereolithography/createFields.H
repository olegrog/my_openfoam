volScalarField totalExposure
(
    IOobject
    (
        "totalExposure",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimEnergy/dimArea
);

volScalarField totalSqrtExposure
(
    IOobject
    (
        "totalSqrtExposure",
        runTime.timeName(),
        mesh
    ),
    mesh,
    sqrt(dimEnergy/dimArea)
);

volScalarField polymerization
(
    IOobject
    (
        "polymerization",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimless
);

Info<< "Reading polymerProperties" << endl;
IOdictionary polymerProperties
(
    IOobject
    (
        "polymerProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


const scalar aConst(polymerProperties.get<scalar>("aConst"));
const scalar bConst(polymerProperties.get<scalar>("bConst"));
const scalar pGel(polymerProperties.get<scalar>("gelPoint"));

const dimensionedScalar Ec("criticalExposure", totalExposure.dimensions(), polymerProperties);
const dimensionedScalar Dp("penetrationDepth", dimLength, polymerProperties);

const laserScanner laser(mesh);
