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

#include "LiquidFraction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LiquidFraction::LiquidFraction
(
    const volScalarField& metalFraction,
    const volScalarField& h,
    const volScalarField& hAtMelting,
    const dimensionedScalar& Hfusion,
    const dimensionedScalar& hAtMeltingPrime,
    const dictionary& dict
)
:
    pSigmoid_(sigmoidFunction::New(dict, 0, 1)),
    metalFraction_(metalFraction),
    h_(h),
    hAtMelting_(hAtMelting),
    Hfusion_(Hfusion),
    hAtMeltingPrime_(hAtMeltingPrime),
    liquidFraction_
    (
        IOobject
        (
            "liquidFraction",
            h.mesh().time().timeName(),
            h.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        h.mesh(),
        dimensionedScalar()
    ),
    inMetal_
    (
        IOobject
        (
            "liquidFractionInMetal",
            h.mesh().time().timeName(),
            h.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        h.mesh(),
        dimless
    ),
    dGasFraction_
    (
        IOobject
        (
            "dLiquidFractionDGasFraction",
            h.mesh().time().timeName(),
            h.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        h.mesh(),
        dimless
    ),
    dEnthalpy_
    (
        IOobject
        (
            "dLiquidFractionDEnthalpy",
            h.mesh().time().timeName(),
            h.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        h.mesh(),
        inv(h.dimensions())
    ),
    wasMelted_
    (
        IOobject
        (
            "wasMelted",
            h.mesh().time().timeName(),
            h.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        h.mesh(),
        dimless
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LiquidFraction::correct()
{
    volScalarField x(2*(h_ - hAtMelting_)/Hfusion_);
    inMetal_ = pSigmoid_->value(x);
    liquidFraction_ = metalFraction_*inMetal_;
    dEnthalpy_ = metalFraction_/Hfusion_*pSigmoid_->derivative(x);
    dGasFraction_ = -pSigmoid_->value(x) - dEnthalpy_*hAtMeltingPrime_;
}


void Foam::LiquidFraction::finalCorrect()
{
    wasMelted_ = Foam::max(wasMelted_, liquidFraction_);
}


// ************************************************************************* //
