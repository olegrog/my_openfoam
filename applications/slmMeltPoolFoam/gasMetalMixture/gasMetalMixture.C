/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020 Oleg Rogozin
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
    liquidFraction_(thermo_.liquidFraction(h_, hAtMelting_, alpha1())),
    hPrimeGasFraction_(thermo_.hPrimeGasFraction(T_, liquidFraction_.dGasFraction())),
    Cp_(thermo_.Cp(T_, liquidFraction_(), alpha2())),
    k_(thermo_.k(T_, liquidFraction_(), alpha2())),
    rho_("rho", alpha1()*rho1() + alpha2()*rho2()),
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

    const scalar maxError = mesh_.solverDict(liquidFraction_().name()).get<scalar>("tolerance");
    scalar error;

    // TODO(olegrog): Use Newton's iterations for boosting
    Info<< "Fixed-point iterations for finding the enthalpy:" << endl;
    do {
        // operator== is used to force the assignment of the boundary field
        h_ == thermo_.h(T_, liquidFraction_(), alpha2());
        liquidFraction_().storePrevIter();
        liquidFraction_.correct();
        const volScalarField errField = mag(liquidFraction_() - liquidFraction_().prevIter());
        error = gMax(errField);
        Info<< liquidFraction_().name() << " error = " << error << endl;
    } while (error > maxError);
    Info<< endl;

    correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::gasMetalMixture::marangoniForce
(
    const volVectorField& gradT
) const
{
    const volVectorField normal(fvc::reconstruct(nHatf()));
    const volTensorField I_nn(tensor::I - sqr(normal));

    return dSigmaDT_*(gradT & I_nn)*mag(fvc::grad(alpha2()));
}


Foam::tmp<Foam::volScalarField> Foam::gasMetalMixture::solidPhaseDamping() const
{
    return mushyCoeff_*alpha1()
        *sqr(1 - liquidFraction_.inMetal())/(sqr(liquidFraction_.inMetal()) + SMALL);
}


void Foam::gasMetalMixture::correct()
{
    immiscibleIncompressibleTwoPhaseMixture::correct();
    liquidFraction_.correct();

    T_ = thermo_.T(h_, hAtMelting_, liquidFraction_(), alpha2());
    hAtMelting_ = thermo_.hAtMelting(alpha2());
    hPrimeGasFraction_ = thermo_.hPrimeGasFraction(T_, liquidFraction_.dGasFraction());
    Cp_ = thermo_.Cp(T_, liquidFraction_(), alpha2());
    k_ = thermo_.k(T_, liquidFraction_(), alpha2());
    rho_ = alpha1()*rho1() + alpha2()*rho2();
}


void Foam::gasMetalMixture::finalCorrect()
{
    liquidFraction_.finalCorrect();
}


// ************************************************************************* //
