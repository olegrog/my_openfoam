/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020-2023 Oleg Rogozin
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
#include "fixedGradientFvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Mixture>
Foam::gasMetalThermalProperties<Mixture>::gasMetalThermalProperties
(
    const fvMesh& mesh,
    const Mixture& mixture
)
:
    thermo_(mesh),
    writeAllFields_(mesh.time().controlDict().getOrDefault("writeAllFields", false)),
    solidificationFields_(thermo_.subDict("solidificationFields")),
    enabledSolidification_(solidificationFields_.getOrDefault("enabled", false)),
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
    liquidFractionInMetal_
    (
        volScalarField::New("liquidFractionInMetal", mesh, dimless)
    ),
    liquidFractionPrimeEnthalpy_
    (
        volScalarField::New("liquidFractionPrimeEnthalpy", mesh, dimMass/dimEnergy)
    ),
    Cp_(thermo_.Cp(T_, liquidFraction_, alphaG_)),
    kappa_(thermo_.kappa(T_, liquidFraction_, alphaG_)),
    HsPrimeAlphaG_(thermo_.HsPrimeAlphaG(T_)),
    tsolidificationGradient_(),
    tsolidificationSpeed_(),
    tredistribution_(),
    tgradT_(),
    updatedRedistribution_(false),
    updatedGradT_(false)
{
    // --- Field initialization

    if (enabledSolidification_)
    {
        tsolidificationGradient_ = tmp<volScalarField>::New
        (
            IOobject
            (
                "solidificationGradient",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTemperature/dimLength)
        );

        tsolidificationSpeed_ = tmp<volScalarField>::New
        (
            IOobject
            (
                "solidificationSpeed",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimVelocity)
        );
    }

    // --- Activate auto-writing of additional fields

    if (writeAllFields_)
    {
        const std::vector<std::reference_wrapper<volScalarField>> scalarFieldsForWriting
        {
            std::ref(Cp_), std::ref(kappa_)
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

        Info<< "Updating the BC for " << h_.name() << nl << endl;
        // operator== is used to force the assignment of the boundary field
        h_ == thermo_.h(T_, liquidFraction_, alphaG_);

        auto& hBf = h_.boundaryFieldRef();
        forAll(hBf, patchi)
        {
            if (isA<fixedGradientFvPatchScalarField>(hBf[patchi]))
            {
                refCast<fixedGradientFvPatchScalarField>(hBf[patchi]).gradient() =
                    hBf[patchi].fvPatchField::snGrad();
            }
            else if (isA<mixedFvPatchScalarField>(hBf[patchi]))
            {
                refCast<mixedFvPatchScalarField>(hBf[patchi]).refValue() = hBf[patchi];
            }
        }
    }

    correctThermo();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Mixture>
void Foam::gasMetalThermalProperties<Mixture>::calcMetalFractions()
{
    const volScalarField x = (h_ - hAtMelting_)/thermo_.Hfusion()/(alphaM_ + SMALL);

    liquidFractionInMetal_ = thermo_.sigmoid().value(x);
    liquidFraction_ = alphaM_*liquidFractionInMetal_;
    //liquidFractionPrimeEnthalpy_ = thermo_.sigmoid().derivative(x)/thermo_.Hfusion();
}


template<class Mixture>
void Foam::gasMetalThermalProperties<Mixture>::calcRedistribution() const
{
    const volScalarField CpM = thermo_.Cp(T_, liquidFraction_, geometricUniformField<scalar>(0));
    const volScalarField CpG = thermo_.Cp(T_, liquidFraction_, geometricUniformField<scalar>(1));
    const volScalarField& rhoM = mixture_.rhoM();
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
    kappa_ = thermo_.kappa(T_, liquidFraction_, alphaG_);
    HsPrimeAlphaG_ = thermo_.HsPrimeAlphaG(T_);

    updatedRedistribution_ = false;
    updatedGradT_ = false;
}


template<class Mixture>
const Foam::volScalarField&
Foam::gasMetalThermalProperties<Mixture>::surfaceHeatSourceRedistribution() const
{
    if (!updatedRedistribution_)
    {
        calcRedistribution();
        updatedRedistribution_ = true;
    }

    return tredistribution_();
}


template<class Mixture>
const Foam::volVectorField& Foam::gasMetalThermalProperties<Mixture>::gradT() const
{
    if (!updatedGradT_)
    {
        tgradT_ = fvc::grad(T_);
        updatedGradT_ = true;
    }

    return tgradT_();
}


template<class Mixture>
void Foam::gasMetalThermalProperties<Mixture>::correctPassiveFields()
{
    if (enabledSolidification_)
    {
        volScalarField& solidificationGradient = tsolidificationGradient_.ref();
        volScalarField& solidificationSpeed = tsolidificationSpeed_.ref();

        const scalar alphaTol = solidificationFields_.getOrDefault("alphaTol", 1.0);
        const volScalarField magGradT = mag(fvc::grad(T_));
        const volScalarField magGradTOld = mag(fvc::grad(T_.oldTime()));

        forAll(T_, cellI)
        {
            if (alphaM_[cellI] < alphaTol) continue;
            if (alphaM_.oldTime()[cellI] < alphaTol) continue;

            const scalar phi = liquidFraction_[cellI]/alphaM_[cellI];
            const scalar phiOld = liquidFraction_.oldTime()[cellI]/alphaM_.oldTime()[cellI];

            if (phiOld > 0.5 && phi < 0.5)
            {
                const scalar deltaT = T_.mesh().time().deltaTValue();
                const scalar coolingRate = (T_.oldTime()[cellI] - T_[cellI])/deltaT;
                const scalar coolingRateTol = SMALL*thermo_.Tmelting().value()/deltaT;

                if (coolingRate < -coolingRateTol)
                {
                    FatalErrorInFunction
                        << " coolingRate = " << coolingRate
                        << ", phi = " << phiOld << " -> " << phi
                        << ", alphaM = " << alphaM_.oldTime()[cellI] << " -> " << alphaM_[cellI]
                        << ", T = " << T_.oldTime()[cellI] << " -> " << T_[cellI]
                        << ", h = " << h_.oldTime()[cellI] << " -> " << h_[cellI]
                        << exit(FatalError);
                }

                const scalar A = (0.5 - phiOld)/(phi - phiOld);
                const scalar B = (phi - 0.5)/(phi - phiOld);
                const scalar oldSpeed = solidificationSpeed[cellI];

                // Linear weighted interpolation
                solidificationGradient[cellI] = A*magGradT[cellI] + B*magGradTOld[cellI];
                solidificationSpeed[cellI] = coolingRate/solidificationGradient[cellI];

                if (oldSpeed > 0)
                {
                    DebugInfo<< "Resolidification in cell #" << cellI
                        << ", alphaM = " << alphaM_[cellI]
                        << ", V = " << oldSpeed
                        << " -> " << solidificationSpeed[cellI] << endl;
                }
            }
        }
    }
}


// ************************************************************************* //
