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

#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LiquidFraction::LiquidFraction
(
    const volScalarField& metalFraction,
    const volScalarField& enthalpy,
    const volScalarField& enthalpyAtFusion,
    const dimensionedScalar& enthalpyFusion,
    const dimensionedScalar& enthalpyAtFusionPrime,
    const dictionary& thermalProperties
)
:
    volScalarField
    (
        IOobject
        (
            "liquidFraction",
            enthalpy.mesh().time().timeName(),
            enthalpy.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        enthalpy.mesh(),
        dimensionedScalar(),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    metalFraction_(metalFraction),
    enthalpy_(enthalpy),
    enthalpyFusion_(enthalpyFusion),
    enthalpyAtFusion_(enthalpyAtFusion),
    enthalpyAtFusionPrime_(enthalpyAtFusionPrime),
    sigmoid_(sigmoidFunction::New(thermalProperties, 0, 1)),
    inMetal_
    (
        IOobject
        (
            "liquidFractionInMetal",
            enthalpy.mesh().time().timeName(),
            enthalpy.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        enthalpy.mesh(),
        dimless
    ),
    dAlpha_
    (
        IOobject
        (
            "liquidFractionDAlpha",
            enthalpy.mesh().time().timeName(),
            enthalpy.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        enthalpy.mesh(),
        dimless
    ),
    dEnthalpy_
    (
        IOobject
        (
            "liquidFractionDEnthalpy",
            enthalpy.mesh().time().timeName(),
            enthalpy.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        enthalpy.mesh(),
        inv(enthalpy.dimensions())
    ),
    wasMelted_
    (
        IOobject
        (
            "wasMelted",
            enthalpy.mesh().time().timeName(),
            enthalpy.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        enthalpy.mesh(),
        dimless
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LiquidFraction::update()
{
    volScalarField x(2*(enthalpy_ - enthalpyAtFusion_) / enthalpyFusion_);
    inMetal_ = (*sigmoid_)(x);
    volScalarField _ = sigmoid_->der(x);
    *this == metalFraction_ * inMetal_;
    dEnthalpy_ = metalFraction_ / enthalpyFusion_ * sigmoid_->der(x);
    dAlpha_ = -(*sigmoid_)(x) - dEnthalpy_ * enthalpyAtFusionPrime_;
}

void Foam::LiquidFraction::finalUpdate()
{
    wasMelted_ = Foam::max(wasMelted_, *this);
}

