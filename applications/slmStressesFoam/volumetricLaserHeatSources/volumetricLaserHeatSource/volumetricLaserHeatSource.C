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
    This file is part of slmStressesFoam.

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

#include "volumetricLaserHeatSource.H"

#include "error.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(volumetricLaserHeatSource);
    defineRunTimeSelectionTable(volumetricLaserHeatSource, mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::volumetricLaserHeatSource> Foam::volumetricLaserHeatSource::New
(
    const fvMesh& mesh
)
{
    const IOdictionary dict
    (
        IOobject
        (
            "laserProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false // Do not register
        )
    );

    const word modelType(dict.subDict("source").get<word>("type"));

    Info<< "Selecting volumetric laser heat source model " << modelType << endl;

    const auto cstrIter = meshConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "volumetricLaserHeatSource",
            modelType,
            *meshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<volumetricLaserHeatSource>(cstrIter()(mesh));
}


Foam::volumetricLaserHeatSource::volumetricLaserHeatSource
(
    const fvMesh& mesh
)
:
    laserHeatSource(mesh),
    sourceDict_(subDict("source")),
    A_(sourceDict_.get<scalar>("absorptivity"))
{}


// ************************************************************************* //
