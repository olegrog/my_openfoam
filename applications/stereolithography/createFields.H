volScalarField totalExposure = volScalarField::New
(
    "totalExposure",
    mesh,
    dimensionedScalar(dimEnergy/dimArea)
);

volScalarField totalSqrtExposure = volScalarField::New
(
    "totalSqrtExposure",
    mesh,
    dimensionedScalar(sqrt(dimEnergy/dimArea))
);

volScalarField polymerization = volScalarField::New
(
    "polymerization",
    mesh,
    dimensionedScalar()
);

volScalarField nLayer = volScalarField::New
(
    "nLayer",
    mesh,
    dimensionedScalar(dimless, 1)
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
const scalar pGel(polymerProperties.get<scalar>("pGel"));

const dimensionedScalar Ec("criticalExposure", totalExposure.dimensions(), polymerProperties);
const dimensionedScalar Dp("penetrationDepth", dimLength, polymerProperties);

const laserScanner laser(mesh);
