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
    This file is part of solidificationFoam.

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

#include "alloyComponent.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(componentPhase);
    defineRunTimeSelectionTable(componentPhase, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::componentPhase> Foam::componentPhase::New
(
    const dictionaryEntry& entry,
    const word& modelType,
    const fvMesh& mesh,
    const alloyComponent& component
)
{
    Info<< "Selecting componentPhase " << modelType << " for " << entry.keyword()
        << " " << component.keyword() << endl;

    const auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            entry,
            "componentPhase",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<componentPhase>(cstrIter()(entry, mesh, component));
}


Foam::componentPhase::componentPhase
(
    const dictionaryEntry& entry,
    const fvMesh& mesh,
    const alloyComponent& component
)
:
    name_(entry.keyword()),
    diffusion_("diffusion", dimArea/dimTime, entry),
    dict_(entry),
    mesh_(mesh),
    component_(component)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::word Foam::componentPhase::name() const
{
    return keyword() + component_.keyword();
}


// ************************************************************************* //
