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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(laserProperties);
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
    radius_("radius", dimLength, *this),
    velocity_("velocity", dimVelocity, *this),
    coordStart_("coordStart", dimLength, *this),
    timeStop_("timeStop", dimTime, *this),
    beamPtr_(laserBeam::New(subDict("beam"), mesh, *this))
{
    const boundBox& bounds = mesh.bounds();
    Info<< " -- Dimensions of the computational domain = " << bounds.max() - bounds.min() << endl
        << " -- Length of the printed track = " << (velocity_*timeStop_).value() << endl;
}


// ************************************************************************* //
