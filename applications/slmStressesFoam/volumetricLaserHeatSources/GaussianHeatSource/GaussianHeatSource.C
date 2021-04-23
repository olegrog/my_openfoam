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

#include "GaussianHeatSource.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(GaussianHeatSource);
    addToRunTimeSelectionTable(volumetricLaserHeatSource, GaussianHeatSource, mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GaussianHeatSource::GaussianHeatSource(const fvMesh& mesh)
:
    volumetricLaserHeatSource(mesh)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GaussianHeatSource::calcSource()
{
    using constant::mathematical::piByTwo;

    const vector& n = beam().direction();
    const dimensionedScalar& R = radius();
    const dimensionedVector& x0 = position();
    const volVectorField& x = mesh_.C();
    const volScalarField z = (x - x0) & n;

    source_ = 2*pos0(z)*A_*beam().I()*exp(-2*magSqr(z/R))/sqrt(piByTwo)/R;
}


// ************************************************************************* //
