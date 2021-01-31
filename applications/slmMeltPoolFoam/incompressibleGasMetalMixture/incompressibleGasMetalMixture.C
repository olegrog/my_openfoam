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

#include "incompressibleGasMetalMixture.H"

#include "fvcGrad.H"
#include "fvcReconstruct.H"
#include "constants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleGasMetalMixture, 0);
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

    // --- Activate auto-writing of additional fields
    if (writeProperties_)
    {
        // nu_ is defined in incompressibleTwoPhaseMixture
        nu_.writeOpt() = IOobject::AUTO_WRITE;
    }

    // --- Checks

    if (rho1() < rho2())
    {
        FatalError << "Metal is lighter than ambient gas." << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::incompressibleGasMetalMixture::marangoniForce() const
{
    const volVectorField normal(fvc::reconstruct(nHatf()));
    const volTensorField I_nn(tensor::I - sqr(normal));
    const volVectorField gradAlphaM = fvc::grad(alphaM_);

    const volVectorField gradT = fvc::grad(T_);

    // This is an alternative formula:
    //  gradT = TPrimeEnthalpy()*fvc::grad(h_) + TPrimeMetalFraction()*gradAlphaM;

    return dSigmaDT_*(gradT & I_nn)*mag(gradAlphaM);
}


Foam::tmp<Foam::volScalarField> Foam::incompressibleGasMetalMixture::solidPhaseDamping() const
{
    volScalarField liquidFractionInMetal =
        thermo_.sigmoid().value((h_ - hAtMelting_)/thermo_.Hfusion()/(alphaM_ + SMALL));
    return mushyCoeff_*alphaM_
        *sqr(1 - liquidFractionInMetal)/(sqr(liquidFractionInMetal) + SMALL);
}


Foam::tmp<Foam::volScalarField> Foam::incompressibleGasMetalMixture::vapourPressure
(
    const dimensionedScalar& p0
) const
{
    using constant::physicoChemical::R;
    return p0*exp(thermo_.metalM()*thermo_.Hvapour()/R*(1/thermo_.Tboiling() - 1/T_));
}


void Foam::incompressibleGasMetalMixture::correct()
{
    immiscibleIncompressibleTwoPhaseMixture::correct();
    gasMetalThermalProperties::correct();
}


// ************************************************************************* //
