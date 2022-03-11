/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020-2022 Oleg Rogozin
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

#include "fvcGrad.H"
#include "fvcReconstruct.H"
#include "constants.H"

#include "generateGeometricField.H"

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
    mushyCoeff_("mushyCoeff", *this)
{
    const scalar Tmelting = thermo().Tmelting().value();

    Info<< endl
        << " -- Surface tension at Tmelting = " << sigmaPtr_->value(Tmelting) << endl
        << " -- Marangoni coefficient at Tmelting = " << dSigmaDT(Tmelting) << endl
        << endl;

    // Activate auto-writing of additional fields
    if (writeAllFields_)
    {
        // nu_ is defined in incompressibleTwoPhaseMixture
        nu_.writeOpt() = IOobject::AUTO_WRITE;
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
    volScalarField liquidFractionInMetal =
        thermo_.sigmoid().value((h() - hAtMelting())/thermo_.Hfusion()/(alphaM_ + SMALL));
    return mushyCoeff_*alphaM_
        *sqr(1 - liquidFractionInMetal)/(sqr(liquidFractionInMetal) + SMALL);
}


Foam::tmp<Foam::volScalarField> Foam::incompressibleGasMetalMixture::vapourPressure
(
    const dimensionedScalar& p0
) const
{
    using constant::physicoChemical::R;
    return p0*exp(thermo_.metalM()*thermo_.Hvapour()/R*(1/thermo_.Tboiling() - 1/T()));
}


void Foam::incompressibleGasMetalMixture::correct()
{
    immiscibleIncompressibleTwoPhaseMixture::correct();
    gasMetalThermalProperties::correct();
}


// ************************************************************************* //
