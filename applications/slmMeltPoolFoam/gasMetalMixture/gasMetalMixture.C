/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020-2021 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of slmMeltPoolFoam.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gasMetalMixture.H"

#include "fvcGrad.H"
#include "fvcReconstruct.H"
#include "constants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gasMetalMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gasMetalMixture::gasMetalMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    immiscibleIncompressibleTwoPhaseMixture(U, phi),
    mesh_(U.mesh()),
    thermo_(phase1Name(), phase2Name(), mesh_),
    writeProperties_(mesh_.time().controlDict().get<bool>("writeProperties")),
    T_
    (
        IOobject
        (
            "T",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    h_
    (
        IOobject
        (
            "h",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimEnergy/dimMass,
        T_.primitiveField(),
        T_.boundaryField()   // copy BC from the temperature field
    ),
    hAtMelting_(thermo_.hAtMelting(alpha2())),
    hAtBoiling_(thermo_.hAtBoiling(alpha2())),
    liquidFraction_
    (
        IOobject
        (
            "liquidFraction",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar()
    ),
    vapourFraction_
    (
        IOobject
        (
            "vapourFraction",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar()
    ),
    Cp_(volScalarField::New("Cp", mesh_, dimSpecificHeatCapacity)),
    k_(volScalarField::New("k", mesh_, dimViscosity*dimDensity*Cp_.dimensions())),
    rho_("rho", alpha1()*rho1_ + alpha2()*rho2_),
    TPrimeEnthalpy_
    (
        volScalarField::New("TPrimeEnthalpy", mesh_, dimTemperature*dimMass/dimEnergy)
    ),
    TPrimeMetalFraction_
    (
        volScalarField::New("TPrimeMetalFraction", mesh_, dimTemperature)
    ),
    pSigma_(Function1<scalar>::New("sigma", this->subDict("sigma"))),
    dSigmaDT_
    (
        "dSigmaDT",
        this->sigmaK()().dimensions()/dimTemperature*dimLength,
        pSigma_->value(1) - pSigma_->value(0)
    ),
    mushyCoeff_("mushyCoeff", *this)
{
    // --- Diagnostic info

    Info<< endl
        << "Surface tension = " << pSigma_->value(0) << " + ("
        << dSigmaDT_.value() << ")*T\n" << endl;

    // --- Diagnostic fields

    if (writeProperties_)
    {
        // nu_ is defined in incompressibleTwoPhaseMixture
        const std::vector<std::reference_wrapper<volScalarField>> scalarFieldsForDiag
        {
            std::ref(Cp_), std::ref(k_), std::ref(rho_), std::ref(nu_)
        };

        for (auto& field : scalarFieldsForDiag)
        {
            field.get().writeOpt() = IOobject::AUTO_WRITE;
        }
    }

    // --- Checks

    if (rho1() < rho2())
    {
        FatalError
            << "Metal is lighter than ambient gas."
            << exit(FatalError);
    }

    // --- Cycle for initial conditions

    const dictionary& h0Dict = mesh_.solverDict(h_.name() + "corr");
    const scalar tolerance = h0Dict.get<scalar>("tolerance");
    const scalar maxIter = h0Dict.getOrDefault<label>("maxIter", 1000);
    label nIter = 0;
    scalar residual;

    Info<< "Fixed-point iterations for correcting enthalpy:" << endl;
    // TODO(olegrog): Use Newton's iterations for boosting
    do
    {
        h_.storePrevIter();
        h_ = thermo_.h(T_, liquidFraction_, vapourFraction_, alpha2());
        liquidFraction_ = thermo_.liquidFraction(h_, hAtMelting_, alpha1());
        vapourFraction_ = thermo_.vapourFraction(h_, hAtBoiling_, alpha1());

        const volScalarField residualField = mag(h_ - h_.prevIter());
        residual = gMax(residualField);

        if (debug)
        {
            Info<< " -- residual = " << residual << ", iteration = " << nIter + 1 << endl;
        }
    }
    while (++nIter < maxIter && residual > tolerance);

    if (residual > tolerance)
    {
        Info<< " -- not converged within " << nIter << " iterations, final residual = "
            << residual << " > " << tolerance << endl;
    }
    else
    {
        Info<< " -- converged in " << nIter << " iterations, final residual = "
            << residual << " < " << tolerance << endl;
    }

    if (!h_.typeHeaderOk<volScalarField>())
    {
        Info<< "Updating the BC for " << h_.name() << endl;
        // operator== is used to force the assignment of the boundary field
        h_ == thermo_.h(T_, liquidFraction_, vapourFraction_, alpha2());
    }

    Info<< endl;

    // --- Correct all fields

    correctThermo();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::gasMetalMixture::marangoniForce() const
{
    const volVectorField normal(fvc::reconstruct(nHatf()));
    const volTensorField I_nn(tensor::I - sqr(normal));

    // gradT = fvc::grad(T_) is a less smooth alternative
    const volVectorField gradT =
        TPrimeEnthalpy_*fvc::grad(h_)
      + TPrimeMetalFraction_*fvc::grad(alpha1());

    return dSigmaDT_*(gradT & I_nn)*mag(fvc::grad(alpha2()));
}


Foam::tmp<Foam::volScalarField> Foam::gasMetalMixture::solidPhaseDamping() const
{
    return mushyCoeff_*alpha1()
        *sqr(alpha1() - liquidFraction_)/(sqr(liquidFraction_) + SMALL);
}


Foam::tmp<Foam::volScalarField> Foam::gasMetalMixture::vapourPressure
(
    const dimensionedScalar& p0
) const
{
    using constant::physicoChemical::R;
    return p0*exp(thermo_.metalM()*thermo_.Hvapour()/R*(1/thermo_.Tboiling() - 1/T_));
}


void Foam::gasMetalMixture::correct()
{
    immiscibleIncompressibleTwoPhaseMixture::correct();

    hAtMelting_ = thermo_.hAtMelting(alpha2());
    hAtBoiling_ = thermo_.hAtBoiling(alpha2());
    // rho_ = alpha1()*rho1_ + alpha2()*rho2_; // already evaluated in alphaEqnSubCycle.H
}


void Foam::gasMetalMixture::calcDerivatives()
{
    const dimensionedScalar& Hfus = thermo_.Hfusion();
    const dimensionedScalar& Hvap = thermo_.Hvapour();

    TPrimeEnthalpy_ =
        (
            1
          - Hfus*thermo_.liquidFractionPrimeEnthalpy(h_, hAtMelting_, alpha1())
          - Hvap*thermo_.vapourFractionPrimeEnthalpy(h_, hAtBoiling_, alpha1())
        )/Cp_;

    TPrimeMetalFraction_ =
        (
            thermo_.sensibleEnthalpyPrimeGasFraction(T_)
          + Hfus*thermo_.liquidFractionPrimeGasFraction(h_, hAtMelting_, alpha1())
          + Hvap*thermo_.vapourFractionPrimeGasFraction(h_, hAtBoiling_, alpha1())
        )/Cp_;
}


void Foam::gasMetalMixture::correctThermo()
{
    liquidFraction_ = thermo_.liquidFraction(h_, hAtMelting_, alpha1());
    vapourFraction_ = thermo_.vapourFraction(h_, hAtBoiling_, alpha1());

    T_ = thermo_.T(h_, hAtMelting_, liquidFraction_, vapourFraction_, alpha2());
    Cp_ = thermo_.Cp(T_, liquidFraction_, alpha2());
    k_ = thermo_.k(T_, liquidFraction_, alpha2());

    calcDerivatives();
}


// ************************************************************************* //
