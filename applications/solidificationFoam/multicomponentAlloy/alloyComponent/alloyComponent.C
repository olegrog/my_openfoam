/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019-2021 Oleg Rogozin
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::alloyComponent::alloyComponent
(
    const dictionaryEntry& entry,
    const fvMesh& mesh,
    const multicomponentAlloy& alloy
)
:
    volScalarField
    (
        IOobject
        (
            "concentration" + entry.keyword(),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.lookupObject<volScalarField>("phase")
    ),
    fieldByKPtr_(nullptr),
    name_(entry.keyword()),
    molarMass_("molarMass", dimMass/dimMoles, entry),
    alloy_(alloy),
    phases_
    (
        entry.lookup("phases"),
        componentPhase::iNew(entry.dictionary::get<word>("type"), mesh, *this)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField& Foam::alloyComponent::byK()
{
    if (!fieldByKPtr_)
    {
        fieldByKPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    name() + "byK",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                *this
            )
        );
    }

    return *fieldByKPtr_;
}


// ************************************************************************* //
