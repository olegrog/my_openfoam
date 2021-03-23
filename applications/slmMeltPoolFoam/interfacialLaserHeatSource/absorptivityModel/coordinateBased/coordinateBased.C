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

#include "interfacialLaserHeatSource.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace absorptivity
    {
        defineTypeName(coordinateBased);
        addToRunTimeSelectionTable(absorptivityModel, coordinateBased, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::absorptivity::coordinateBased::coordinateBased(const dictionary& dict)
:
    absorptivityModel(dict),
    beamDirection_(dict.get<vector>("beamDirection"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::absorptivity::coordinateBased::value
(
    const volScalarField& alphaM,
    const volVectorField& gradAlphaM,
    const interfacialLaserHeatSource& laser
) const
{
    const fvMesh& mesh = alphaM.mesh();
    return max(dimensionedScalar(gradAlphaM.dimensions()), gradAlphaM & beamDirection_)
        *(1 - (1 - A_)*exp(min(Zero, -(mesh.C() & beamDirection_)/2/laser.radius())));
}


// ************************************************************************* //
