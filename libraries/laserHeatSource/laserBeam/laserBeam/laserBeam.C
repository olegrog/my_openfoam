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

#include "laserBeam.H"

#include "error.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(laserBeam);
    defineRunTimeSelectionTable(laserBeam, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::laserBeam> Foam::laserBeam::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const laserProperties& laser
)
{
    const word beamType(dict.get<word>("type"));

    Info<< "Selecting laser beam " << beamType << endl;

    const auto cstrIter = dictionaryConstructorTablePtr_->cfind(beamType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "laserBeam",
            beamType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<laserBeam>(cstrIter()(dict, mesh, laser));
}


Foam::laserBeam::laserBeam
(
    const dictionary& dict,
    const fvMesh& mesh,
    const laserProperties& laser
)
:
    mesh_(mesh),
    laser_(laser),
    direction_(dict.get<vector>("direction").normalise())
{
    if (mag(direction_) < SMALL)
    {
        FatalError << "The beam direction is not specified." << exit(FatalError);
    }
}


// ************************************************************************* //
