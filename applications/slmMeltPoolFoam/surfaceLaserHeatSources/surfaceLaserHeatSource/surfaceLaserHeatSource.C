/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2021-2023 Oleg Rogozin
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

#include "surfaceLaserHeatSource.H"

#include "error.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(surfaceLaserHeatSource);
    defineRunTimeSelectionTable(surfaceLaserHeatSource, mixtureAdvector);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceLaserHeatSource> Foam::surfaceLaserHeatSource::New
(
    const incompressibleGasMetalMixture& mixture,
    const isoAdvection& advector
)
{
    const IOdictionary dict
    (
        IOobject
        (
            "laserProperties",
            mixture.U().time().constant(),
            mixture.U().db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    const word modelType(dict.subDict("source").get<word>("type"));

    Info<< "Selecting surface laser heat source model " << modelType << endl;

    const auto cstrIter = mixtureAdvectorConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "surfaceLaserHeatSource",
            modelType,
            *mixtureAdvectorConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<surfaceLaserHeatSource>(cstrIter()(mixture, advector));
}


Foam::surfaceLaserHeatSource::surfaceLaserHeatSource
(
    const word& modelType,
    const incompressibleGasMetalMixture& mixture,
    const isoAdvection& advector
)
:
    laserHeatSource(advector.alpha().mesh()),
    elapsedTime_(0),
    startTime_(0),
    modelDict_(subDict("source").subDict(modelType + "Model")),
    mixture_(mixture),
    advector_(advector)
{}


// ************************************************************************* //
