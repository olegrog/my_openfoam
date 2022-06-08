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
    This file is part of laserHeatSource.

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

#include "laserProperties.H"

#include "Constant.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(laserProperties, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laserProperties::laserProperties(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "laserProperties",
            mesh.time().constant(),
            mesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    time_(mesh.time()),
    powerPtr_(Function1<scalar>::New("power", *this)),
    radiusPtr_(Function1<scalar>::New("radius", *this)),
    velocityPtr_(Function1<vector>::New("velocity", *this)),
    coordStart_("coordStart", dimLength, *this),
    timeStop_("timeStop", dimTime, *this),
    beamPtr_(laserBeam::New(subDict("beam"), mesh, *this)),
    laserDict_
    (
        IOobject
        (
            "laser",
            mesh.time().timeName(),
            "uniform",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    position_("position", dimLength, coordStart_.value(), laserDict_),
    MRF_(mesh.lookupClass<movingReferenceFrame>().begin().val())
{
    if (isA<Function1Types::Constant<vector>>(velocityPtr_()))
    {
        const boundBox& bounds = mesh.bounds();
        Info<< " -- Dimensions of the computational domain = "
            << bounds.max() - bounds.min() << endl
            << " -- Length of the printed track = "
            << velocityPtr_->value(time_.value())*timeStop_.value() << endl;
    }

    if (!MRF_)
    {
        WarningInFunction << "MRF has not found" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedVector Foam::laserProperties::velocity() const
{
    vector value = velocityPtr_->value(time_.value());
    if (MRF_)
    {
        value -= MRF_->Urel().value();
    }
    return dimensionedVector("velocity", dimVelocity, value);
}


void Foam::laserProperties::update()
{
    const dimensionedVector vel = velocity();
    position_ += vel*time_.deltaT();

    DebugInfo
        << "Laser: power = " << power().value()*switchedOn() << " radius = " << radius().value()
        << " velocity = " << vel.value() << " position = " << position_.value() << endl;

    laserDict_.set(power().name(), power().value());
    laserDict_.set(radius().name(), radius().value());
    laserDict_.set(vel.name(), vel.value());
    laserDict_.set(position_.name(), position_.value());
    laserDict_.set("switchedOn", switchedOn());
}


// ************************************************************************* //
