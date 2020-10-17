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
    This file is part of stereolithography.

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

#include "laserScanner.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laserScanner::laserScanner(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "laserProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    radius_("radius", dimLength, *this),
    power_("power", dimPower, *this),
    hatchingSpeed_("hatchingSpeed", dimArea/dimTime, *this),
    layerThickness_("layerThickness", dimLength, *this),
    firstLayerThickness_("firstLayerThickness", dimLength, *this)
{
    Info << " -- Laser exposure = " << E().value() << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::laserScanner::height(label i) const
{
    if (i < 1)
    {
        FatalError
            << "Layers are numbered starting from one."
            << exit(FatalError);
    }

    return firstLayerThickness_ + layerThickness_*(i - 1);
}


// ************************************************************************* //
