/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2022 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of movingReferenceFrame.

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

#include "givenFrame.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(givenFrame);
    addToRunTimeSelectionTable(movingReferenceFrame, givenFrame, mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::givenFrame::givenFrame(const fvMesh& mesh)
:
    movingReferenceFrame(mesh),
    dict_(subDict(typeName + "Frame")),
    UrelFunc_(Function1<vector>::New("Urel", dict_, &mesh))
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::givenFrame::evaluate()
{
    const scalar t = mesh_.time().value();
    Urel_.value() = UrelFunc_->value(t);
}


// ************************************************************************* //
