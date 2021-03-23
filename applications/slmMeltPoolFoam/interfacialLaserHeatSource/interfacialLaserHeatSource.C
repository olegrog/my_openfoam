/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2021 Oleg Rogozin
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

#include "interfacialLaserHeatSource.H"

#include "fvcGrad.H"
#include "constants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfacialLaserHeatSource, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::interfacialLaserHeatSource::surfaceGaussian
(
    const volVectorField& x,
    dimensionedVector x0
) const
{
    using constant::mathematical::pi;
    volVectorField r = x - x0;
    r.replace(vector::Z, 0);
    return 2*exp(-2*magSqr(r/radius_))/pi/sqr(radius_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacialLaserHeatSource::interfacialLaserHeatSource
(
    const volScalarField& alphaM,
    const volVectorField& gradAlphaM
)
:
    IOdictionary
    (
        IOobject
        (
            "laserProperties",
            alphaM.mesh().time().constant(),
            alphaM.mesh().thisDb(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(alphaM.mesh()),
    time_(alphaM.mesh().time()),
    powerDistribution_
    (
        volScalarField::New("powerDistribution", mesh_, dimPower/dimArea)
    ),
    absorptivity_
    (
        volScalarField::New("absorptivity", mesh_, inv(dimLength))
    ),
    alphaM_(alphaM),
    gradAlphaM_(gradAlphaM),
    radius_("radius", dimLength, *this),
    power_("power", dimPower, *this),
    velocity_("velocity", dimVelocity, *this),
    coordStart_("coordStart", dimLength, *this),
    timeStop_("timeStop", dimTime, *this),
    absorptivityModelPtr_(absorptivityModel::New(subDict("absorptivity")))
{
    if (time_.controlDict().get<bool>("writeProperties"))
    {
        powerDistribution_.writeOpt() = IOobject::AUTO_WRITE;
        absorptivity_.writeOpt() = IOobject::AUTO_WRITE;
    }

    const boundBox& bounds = mesh_.bounds();
    Info<< "Dimensions of the computational domain = " << bounds.max() - bounds.min() << endl
        << "Length of the printed track = " << (velocity_*timeStop_).value() << endl
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::interfacialLaserHeatSource::update()
{
    powerDistribution_ = (time_ < timeStop_)*power_*surfaceGaussian(mesh_.C(), position());
}


void Foam::interfacialLaserHeatSource::correct()
{
    absorptivity_ = absorptivityModelPtr_->value(alphaM_, gradAlphaM_, *this);
}


// ************************************************************************* //
