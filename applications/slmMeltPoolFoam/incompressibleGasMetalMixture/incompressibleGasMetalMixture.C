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

#include "incompressibleGasMetalMixture.H"

#include "fvc.H"
#include "constants.H"

#include "generateGeometricField.H"
#include "updateGeometricField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleGasMetalMixture, 0);

    defineTemplateTypeNameAndDebugWithName
    (
        gasMetalThermalProperties<incompressibleGasMetalMixture>,
        "incompressibleGasMetalThermalProperties",
        0
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleGasMetalMixture::incompressibleGasMetalMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    immiscibleIncompressibleTwoPhaseMixture(U, phi),
    gasMetalThermalProperties(U.mesh(), *this),
    sigmaPtr_(Function1<scalar>::New("sigma", this->subDict("sigma"))),
    mushyCoeff_("mushyCoeff", dimDensity/dimTime, *this),
    Tcritical_("Tcritical", dimTemperature, thermo().subDict("metal")),
    metalDict_(subDict("metal")),
    quasiIncompressible_(metalDict_.getOrDefault("quasiIncompressible", false)),
    rhoJump_(rho1_.value() - metalDict_.get<scalar>("rhoSolid")),
    dRhoMDTSolid_(metalDict_.lookup(IOobject::groupName("dRhoDT", "solid"))),
    dRhoMDTLiquid_(metalDict_.lookup(IOobject::groupName("dRhoDT", "liquid"))),
    rhoM_(volScalarField::New("rhoM", U.mesh(), rho1_)),
    divPhi_(volScalarField::New("divPhi", U.mesh(), dimensionedScalar(inv(dimTime))))
{
    const scalar Tmelting = thermo().Tmelting().value();

    Info<< "Transport properties:" << endl
        << " -- Surface tension at Tmelting = " << sigmaPtr_->value(Tmelting) << endl
        << " -- Marangoni coefficient at Tmelting = " << dSigmaDT(Tmelting) << endl
        << " -- Gas density = " << rho2_.value() << endl;

    if (quasiIncompressible_)
    {
        Info<< " -- Liquid metal density at " << Tmelting + 1000 << " = "
            << rhoM(Tmelting, Tmelting + 1000, 1) << endl
            << " -- Liquid metal density at Tmelting = " << rhoM(Tmelting, Tmelting, 1) << endl
            << " -- Solid metal density at Tmelting = " << rhoM(Tmelting, Tmelting, 0) << endl
            << " -- Solid metal density at " << Tmelting - 1000 << " = "
            << rhoM(Tmelting, Tmelting - 1000, 0) << endl;
    }
    else
    {
        Info<< " -- Metal density = " << rho1_.value() << endl;
    }
    Info<< endl;

    // Activate auto-writing of additional fields
    if (writeAllFields_)
    {
        // nu_ is defined in incompressibleTwoPhaseMixture
        nu_.writeOpt() = IOobject::AUTO_WRITE;
        divPhi_.writeOpt() = IOobject::AUTO_WRITE;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::incompressibleGasMetalMixture::dSigmaDT(scalar T) const
{
    // Function1 does not support derivatives; therefore, we compute them manually
    const scalar dT = ROOTSMALL;

    return (sigmaPtr_->value(T+dT) - sigmaPtr_->value(T-dT))/2/dT;
}


Foam::tmp<Foam::volScalarField> Foam::incompressibleGasMetalMixture::dSigmaDT() const
{
    return generateGeometricField<volScalarField>
    (
        "dSigmaDT",
        T().mesh(),
        dimEnergy/dimArea/dimTemperature,
        [this](scalar T)
        {
            return dSigmaDT(T);
        },
        T()
    );
}


Foam::tmp<Foam::volScalarField> Foam::incompressibleGasMetalMixture::dRhoMDT() const
{
    const scalar Tm = thermo().Tmelting().value();

    return generateGeometricField<volScalarField>
    (
        "dRhoMDT",
        T().mesh(),
        dimDensity/dimTemperature,
        [this, Tm](scalar T, scalar phi)
        {
            return T <= Tm ? dRhoMDTSolid_.value(T) : dRhoMDTLiquid_.value(T);
        },
        T(), liquidFraction()
    );
}


Foam::scalar Foam::incompressibleGasMetalMixture::rhoM(scalar Tm, scalar T, scalar phi) const
{
    scalar piecewise = T <= Tm ? dRhoMDTSolid_.integral(Tm, T) : dRhoMDTLiquid_.integral(Tm, T);
    return rho1_.value() + piecewise - (1 - phi)*rhoJump_;
}


void Foam::incompressibleGasMetalMixture::updateRhoM()
{
    const scalar Tm = thermo().Tmelting().value();

    updateGeometricField<volScalarField>
    (
        rhoM_, [this, Tm](scalar& f, scalar T, scalar phi)
        {
            f = rhoM(Tm, T, phi);
        },
        T(), liquidFractionInMetal()
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::incompressibleGasMetalMixture::marangoniForce() const
{
    const volVectorField normal(fvc::reconstruct(nHatf()));
    const volTensorField I_nn(tensor::I - sqr(normal));

    const volVectorField& gradAlphaM = gasMetalThermalProperties::gradAlphaM();
    const volVectorField& gradT = gasMetalThermalProperties::gradT();

    // This is an alternative formula:
    //  gradT = TPrimeEnthalpy()*fvc::grad(h_) + TPrimeMetalFraction()*gradAlphaM;

    return dSigmaDT()*(gradT & I_nn)*mag(gradAlphaM);
}


Foam::tmp<Foam::volScalarField> Foam::incompressibleGasMetalMixture::solidPhaseDamping() const
{
    return mushyCoeff_*alphaM_
        *sqr(1 - liquidFractionInMetal())/(sqr(liquidFractionInMetal()) + SMALL);
}


Foam::tmp<Foam::volScalarField> Foam::incompressibleGasMetalMixture::vapourPressure
(
    const dimensionedScalar& p0
) const
{
    using constant::physicoChemical::R;
    return
        p0*exp(thermo_.metalM()*thermo_.Hvapour()/R
       *(1/thermo_.Tboiling() - 1/min(T(), Tcritical_)));
}


const Foam::volScalarField& Foam::incompressibleGasMetalMixture::divPhi()
{
    if (quasiIncompressible_)
    {
        const dimensionedScalar rhoJump("rhoJump", dimDensity, rhoJump_);
        divPhi_ =
        (
          - rhoJump*fvc::DDt(phi_, liquidFractionInMetal())
          - dRhoMDT()*fvc::DDt(phi_, T())
        )*alphaM_/rhoM_;
    }

    return divPhi_;
}


void Foam::incompressibleGasMetalMixture::correct()
{
    // incompressibleTwoPhaseMixture -> nu
    // interfaceProperties -> K = curvature, nHatf
    immiscibleIncompressibleTwoPhaseMixture::correct();
    // gasMetalThermalProperties -> hAtMelting, gradAlphaM
    gasMetalThermalProperties::correct();
}


void Foam::incompressibleGasMetalMixture::correctThermo()
{
    gasMetalThermalProperties::correctThermo();

    if (quasiIncompressible_)
    {
        updateRhoM();
    }
}


// ************************************************************************* //
