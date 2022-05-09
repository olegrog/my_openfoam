Info<< "Reading field lambda\n" << endl;
volScalarField lambda
(
    IOobject
    (
        "lambda",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
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

Info<< "Reading ZNDdict\n" << endl;
IOdictionary ZNDdict
(
    IOobject
    (
        "ZNDdict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const vector nHat(ZNDdict.get<vector>("normal").normalise());
const dimensionedVector origin("origin", dimLength, ZNDdict);

Info<< "Reading chemical properties\n" << endl;
IOdictionary chemicalProperties
(
    IOobject
    (
        "chemicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

autoPtr<psiThermo> pThermo(psiThermo::New(mesh));
volScalarField gamma = pThermo->Cp()/pThermo->Cv();
chemicalProperties.add("g", gMax(gamma));
scalar WbyR = 1e-3*gMax(pThermo->W()())/constant::physicoChemical::R.value();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pThermo->rho()
);

const ZNDsystem ode(chemicalProperties);

autoPtr<ODESolver> odeSolver = ODESolver::New(ode, ZNDdict);

label timeIndex = 1;
autoPtr<volScalarField> normalisedGradRhoPtr; // Indicator for mesh refinement

if (isA<dynamicRefineFvMesh>(mesh))
{
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).optionalSubDict(mesh.typeName + "Coeffs")
    );
    timeIndex = refineDict.get<label>("refineInterval");

    normalisedGradRhoPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "normalisedGradRho",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimDensity/dimLength
        )
    );
}
