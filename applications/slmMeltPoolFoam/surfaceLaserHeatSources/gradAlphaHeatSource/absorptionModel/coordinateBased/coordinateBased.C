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

#include "coordinateBased.H"

#include "addToRunTimeSelectionTable.H"

#include "laserProperties.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace absorption
    {
        defineTypeName(coordinateBased);
        addToRunTimeSelectionTable(absorptionModel, coordinateBased, dictionary);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::absorption::coordinateBased::A
(
    const volVectorField& gradAlphaM,
    const laserProperties& laser
) const
{
    const fvMesh& mesh = gradAlphaM.mesh();
    const vector& n = laser.beam().direction();

    return max(dimensionedScalar(gradAlphaM.dimensions()), gradAlphaM & n)
        *(1 - (1 - A_)*exp(min(Zero, -(mesh.C() & n)/2/laser.radius())));
}


// ************************************************************************* //
