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

#include "laserHeatSource.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(laserHeatSource);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laserHeatSource::laserHeatSource(const fvMesh& mesh)
:
    laserProperties(mesh),
    mesh_(mesh),
    writeAllFields_(mesh.time().controlDict().get<bool>("writeAllFields")),
    source_
    (
        IOobject
        (
            "laserHeatSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            writeAllFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPower/dimVolume)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::laserHeatSource::correct()
{
    if (switchedOn())
    {
        calcSource();
    }
    else
    {
        source_.primitiveFieldRef() = 0;
    }
}


// ************************************************************************* //
