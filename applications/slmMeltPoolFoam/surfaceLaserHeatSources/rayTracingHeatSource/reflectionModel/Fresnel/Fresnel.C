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

#include "Fresnel.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace reflection
    {
        defineTypeName(Fresnel);
        addToRunTimeSelectionTable(reflectionModel, Fresnel, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reflection::Fresnel::Fresnel(const dictionary& dict)
:
    reflectionModel(dict),
    n_(dict.get<scalar>("n")),
    k_(dict.get<scalar>("k"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::reflection::Fresnel::rho(const scalar incidentAngle) const
{
    scalar sinTheta = sin(incidentAngle);
    scalar cosTheta = cos(incidentAngle);

    scalar r = sqr(n_) - sqr(k_) - sqr(sinTheta);
    scalar sqrP = (sqrt(sqr(r) + sqr(2*n_*k_)) - r)/2;
    scalar sqrQ = (sqrt(sqr(r) + sqr(2*n_*k_)) + r)/2;
    scalar q = sqrt(sqrQ), p = sqrt(sqrP);

    scalar rhoP =
        (sqr((sqr(n_) - sqr(k_))*cosTheta - q) + sqr(2*n_*k_*cosTheta - p))
       /(sqr((sqr(n_) - sqr(k_))*cosTheta + q) + sqr(2*n_*k_*cosTheta + p));

    scalar rhoN = (sqr(q - cosTheta) + sqrP)/(sqr(q + cosTheta) + sqrP);

    return (rhoP + rhoN)/2;
}


Foam::vector Foam::reflection::Fresnel::R
(
    const vector& i,
    const vector& n
) const
{
    return i - 2*n*(i & n);
}


// ************************************************************************* //
