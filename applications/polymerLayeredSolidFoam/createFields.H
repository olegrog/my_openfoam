Info<< "Reading field D\n" << endl;
volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading polymer properties\n" << endl;
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

// Polymerization properties
const scalar aConst(polymerProperties.get<scalar>("aConst"));
const scalar bConst(polymerProperties.get<scalar>("bConst"));
const scalar pGel(polymerProperties.get<scalar>("pGel"));
const dimensionedScalar Ec("criticalExposure", dimEnergy/dimArea, polymerProperties);
const dimensionedScalar Dp("penetrationDepth", dimLength, polymerProperties);

// Mechanical properties
const dimensionedScalar Emax("Emax", dimPressure, polymerProperties);
const dimensionedScalar epsilonChemicalMax("epsilonChemicalMax", dimless, polymerProperties);
const dimensionedScalar E0 = Emax*polymerProperties.get<scalar>("E0ByEmax");

Info<< "Calculating polymerization properties" << endl;

const volVectorField& coord = mesh.C();
const volScalarField Z = coord.component(2);

const boundBox& bounds = mesh.bounds();
const dimensionedScalar ymin("ymin", dimLength, bounds.min().y());
const dimensionedScalar ymax("ymax", dimLength, bounds.max().y());
const dimensionedScalar zmin("zmin", dimLength, bounds.min().z());
const dimensionedScalar zmax("zmax", dimLength, bounds.max().z());
const dimensionedScalar L = ymax - ymin;
const dimensionedScalar small = SMALL*L;

if (mag(zmin) > small)
{
    FatalError
        << "Coordinate zmin is not equal to zero."
        << exit(FatalError);
}

using constant::mathematical::pi;
const laserScanner laser(mesh);
const scalar Fscatt = aConst*log(laser.E()/Ec).value() + bConst;
const dimensionedScalar PConst = -pow025(2*pi)*Foam::log(1 - pGel)/sqrt(Ec*laser.R()/laser.V(L));
const dimensionedScalar radius = laser.R()*Fscatt;

label nLayer = 0;
for (; laser.height(nLayer + 1) < zmax + small; nLayer++);

if (mag(zmax - laser.height(nLayer)) > small)
{
    FatalError
        << "Coordinate zmax does not correspond to integer number of layers:" << nl
        << "  zmax = " << zmax.value() << nl
        << "  laser.height(" << nLayer << ") = " << laser.height(nLayer).value() << nl
        << exit(FatalError);
}

Info<< " -- Scattering factor = " << Fscatt << endl
    << " -- Polymerization constant = " << PConst.value() << endl
    << " -- Number of layers = " << nLayer << nl << endl;

volScalarField totalSqrtExposure
(
    IOobject
    (
        "totalSqrtExposure",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(sqrt(dimEnergy/dimArea))
);

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    neg0(Z - laser.height(1) - small)
);

Info<< "Calculating mechanical properties\n" << endl;

volScalarField E("E", Emax*p);
const volScalarField nu = volScalarField::New
(
    "nu",
    mesh,
    dimensionedScalar("nu", dimless, polymerProperties)
);

volScalarField mu("mu", E/(2.0*(1.0 + nu)));
volScalarField lambda("lambda", nu*E/((1.0 + nu)*(1.0 - 2.0*nu)));
volScalarField threeK("threeK", E/(1.0 - 2.0*nu));
volScalarField active("active", p);

volTensorField gradD("gradD", fvc::grad(D));
volSymmTensorField epsilon("epsilon", symm(gradD));

volSymmTensorField epsilonRes
(
    IOobject
    (
        "epsilonRes",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    0*epsilon
);

volSymmTensorField sigma
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mu*(twoSymm(gradD) + 2*epsilonRes)
  + lambda*(tr(gradD) + tr(epsilonRes))*I
  - threeK*epsilonChemicalMax*p*I
);

volVectorField divSigmaExp
(
    "divSigmaExp",
    fvc::div(sigma) - fvc::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
);

volScalarField sigmaEq
(
    IOobject
    (
        "sigmaEq",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    sigma.dimensions()
);

// Flag to write auxiliary fields
bool writeAllFields(runTime.controlDict().get<bool>("writeAllFields"));

if (writeAllFields)
{
    E.writeOpt() = IOobject::AUTO_WRITE;
    active.writeOpt() = IOobject::AUTO_WRITE;
    sigma.writeOpt() = IOobject::AUTO_WRITE;
    epsilon.writeOpt() = IOobject::AUTO_WRITE;
}
