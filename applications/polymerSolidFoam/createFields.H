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

Info<< "Reading field polymerization\n" << endl;
volScalarField polymerization
(
    IOobject
    (
        "polymerization",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field nLayer\n" << endl;
volScalarField nLayer
(
    IOobject
    (
        "nLayer",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
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

const dimensionedScalar Emax("Emax", dimPressure, polymerProperties);
const dimensionedScalar pGel("pGel", dimless, polymerProperties);
const dimensionedScalar epsilonChemicalMax("epsilonChemicalMax", dimless, polymerProperties);
const dimensionedScalar interlayerCohesion("interlayerCohesion", dimless, polymerProperties);

Info<< "Calculating mechanical properties\n" << endl;

volScalarField E("E", Emax*(polymerization - pGel)/(1 - pGel));

if (gMin(E) < SMALL)
{
    FatalError
        << "Somewhere the Young's modulus is equal to zero."
        << exit(FatalError);
}

volScalarField nu = volScalarField::New
(
    "nu",
    mesh,
    dimensionedScalar("nu", dimless, polymerProperties)
);

Info<< "Calculating Lame's coefficients\n" << endl;

volScalarField mu("mu", E/(2.0*(1.0 + nu)));
volScalarField lambda("lambda", nu*E/((1.0 + nu)*(1.0 - 2.0*nu)));
volScalarField threeK("threeK", E/(1.0 - 2.0*nu));

volTensorField gradD("gradD", fvc::grad(D));
volSymmTensorField epsilon("epsilon", symm(gradD));

symmTensor Ixy = symmTensor::I; Ixy.zz() = 0;

volSymmTensorField sigma
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mu*twoSymm(gradD) + lambda*I*tr(gradD)
  - threeK*epsilonChemicalMax*(polymerization*I + interlayerCohesion*(nLayer - 1)*Ixy)
);

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
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

mesh.setFluxRequired(D.name());

// Flag to write auxiliary fields
bool writeAllFields(runTime.controlDict().get<bool>("writeAllFields"));

if (writeAllFields)
{
    polymerization.writeOpt() = IOobject::AUTO_WRITE;
    nLayer.writeOpt() = IOobject::AUTO_WRITE;
    E.writeOpt() = IOobject::AUTO_WRITE;
    sigma.writeOpt() = IOobject::AUTO_WRITE;
    epsilon.writeOpt() = IOobject::AUTO_WRITE;
}
