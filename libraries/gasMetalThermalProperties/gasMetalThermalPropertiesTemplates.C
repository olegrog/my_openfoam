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
    This file is part of gasMetalThermalProperties.

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

#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Mixture>
Foam::gasMetalThermalProperties<Mixture>::gasMetalThermalProperties
(
    const fvMesh& mesh,
    const Mixture& mixture
)
:
    thermo_(mesh),
    writeProperties_(mesh.time().controlDict().get<bool>("writeProperties")),
    mixture_(mixture),
    alphaM_(mixture.alpha1()),
    alphaG_(mixture.alpha2()),
    gradAlphaM_(fvc::grad(alphaM_, "nHat")),
    T_
    (
        IOobject
        (
            "T",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    h_
    (
        IOobject
        (
            "h",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        T_.primitiveField(),
        T_.boundaryField()   // copy BC from the temperature field
    ),
    hAtMelting_(thermo_.hAtMelting(alphaG_)),
    liquidFraction_
    (
        IOobject
        (
            "liquidFraction",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar()
    ),
    liquidFractionPrimeEnthalpy_
    (
        volScalarField::New("liquidFractionPrimeEnthalpy", mesh, dimMass/dimEnergy)
    ),
    Cp_(thermo_.Cp(T_, liquidFraction_, alphaG_)),
    k_(thermo_.k(T_, liquidFraction_, alphaG_)),
    HsPrimeAlphaG_(thermo_.HsPrimeAlphaG(T_)),
    updatedRedistr_(false),
    tredistribution_()
{
    // --- Activate auto-writing of additional fields

    if (writeProperties_)
    {
        const std::vector<std::reference_wrapper<volScalarField>> scalarFieldsForWriting
        {
            std::ref(Cp_), std::ref(k_)
        };

        for (auto& field : scalarFieldsForWriting)
        {
            field.get().writeOpt() = IOobject::AUTO_WRITE;
        }
    }

    // --- Checks

    if (mixture.phase1Name() != "metal")
    {
        FatalError << "The first phase of the mixture is not metal." << exit(FatalError);
    }

    if (mixture.phase2Name() != "gas")
    {
        FatalError << "The second phase of the mixture is not gas." << exit(FatalError);
    }

    if (mixture.rho1() < mixture.rho2())
    {
        FatalError << "Metal is lighter than ambient gas." << exit(FatalError);
    }

    // --- Correct fields

    if (!h_.typeHeaderOk<volScalarField>())
    {
        const dictionary& h0Dict = mesh.solverDict(h_.name() + "corr");
        const scalar tolerance = h0Dict.get<scalar>("tolerance");
        const scalar maxIter = h0Dict.getOrDefault<label>("maxIter", 1000);
        label nIter = 0;
        scalar residual;

        Info<< "Fixed-point iterations for correcting enthalpy:" << endl;
        do
        {
            h_.storePrevIter();
            h_ = thermo_.h(T_, liquidFraction_, alphaG_);
            calcMetalFractions();

            const volScalarField residualField = mag(h_ - h_.prevIter());
            residual = gMax(residualField);

            DebugInfo << " -- residual = " << residual << ", iteration = " << nIter + 1 << endl;
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

        Info<< "Updating the BC for " << h_.name() << endl;
        // operator== is used to force the assignment of the boundary field
        h_ == thermo_.h(T_, liquidFraction_, alphaG_);
    }

    correctThermo();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Mixture>
void Foam::gasMetalThermalProperties<Mixture>::calcMetalFractions()
{
    const volScalarField x = (h_ - hAtMelting_)/thermo_.Hfusion()/(alphaM_ + SMALL);

    liquidFraction_ = alphaM_*thermo_.sigmoid().value(x);
    liquidFractionPrimeEnthalpy_ = thermo_.sigmoid().derivative(x)/thermo_.Hfusion();
}


template<class Mixture>
void Foam::gasMetalThermalProperties<Mixture>::calcRedistribution() const
{
    const volScalarField CpM = thermo_.Cp(T_, liquidFraction_, geometricUniformField<scalar>(0));
    const volScalarField CpG = thermo_.Cp(T_, liquidFraction_, geometricUniformField<scalar>(1));
    const dimensionedScalar rhoM = mixture_.rho1();
    const dimensionedScalar rhoG = mixture_.rho2();
    const volScalarField rho = alphaM_*rhoM + alphaG_*rhoG;

    tredistribution_ = 2*rho*Cp_/(rhoM*CpM + rhoG*CpG);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Mixture>
void Foam::gasMetalThermalProperties<Mixture>::correct()
{
    hAtMelting_ = thermo_.hAtMelting(alphaG_);

    // "nHat" is used according to the interfaceProperties class
    gradAlphaM_ = fvc::grad(alphaM_, "nHat");
}


template<class Mixture>
void Foam::gasMetalThermalProperties<Mixture>::correctThermo()
{
    calcMetalFractions();

    T_ = thermo_.T(h_, hAtMelting_, liquidFraction_, alphaG_);
    Cp_ = thermo_.Cp(T_, liquidFraction_, alphaG_);
    k_ = thermo_.k(T_, liquidFraction_, alphaG_);
    HsPrimeAlphaG_ = thermo_.HsPrimeAlphaG(T_);

    updatedRedistr_ = false;
}

template<class Mixture>
const Foam::volScalarField&
Foam::gasMetalThermalProperties<Mixture>::surfaceHeatSourceRedistribution() const
{
    if (!updatedRedistr_)
    {
        calcRedistribution();
        updatedRedistr_ = true;
    }

    return tredistribution_();
}


// ************************************************************************* //
