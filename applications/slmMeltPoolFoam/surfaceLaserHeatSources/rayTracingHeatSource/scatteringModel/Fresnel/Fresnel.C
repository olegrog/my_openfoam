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
#include "constants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace scattering
    {
        defineTypeNameAndDebug(Fresnel, 0);
        addToRunTimeSelectionTable(scatteringModel, Fresnel, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scattering::Fresnel::Fresnel(const dictionary& dict)
:
    scatteringModel(dict),
    n_(dict.get<scalar>("n")),
    k_(dict.get<scalar>("k"))
{
    using constant::mathematical::pi;

    Info<< " -- Normal absorptivity = " << 1 - R(1) << endl;
    DebugInfo
        << " -- Absorptivity (pi/6) = " << 1 - R(cos(pi/6)) << nl
        << " -- Absorptivity (pi/3) = " << 1 - R(cos(pi/3)) << nl
        << " -- Grazing absorptivity = " << 1 - R(0) << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::scattering::Fresnel::reflection
(
    const vector& i,
    const vector& n
) const
{
    return i - 2*n*(i & n);
}


Foam::vector Foam::scattering::Fresnel::refraction
(
    const vector& i,
    const vector& n
) const
{
    scalar sinTheta = sqrt(1 - sqr(i & n));
    scalar r = sqr(n_) - sqr(k_) - sqr(sinTheta);
    scalar p = sqrt((sqrt(sqr(r) + sqr(2*n_*k_)) + r)/2);

    return (i + n*(p - (i & n))).normalise();
}


Foam::scalar Foam::scattering::Fresnel::R(scalar cosTheta) const
{
    scalar sinTheta = sqrt(1 - sqr(cosTheta));
    scalar tanTheta = sinTheta/cosTheta;

    scalar r = sqr(n_) - sqr(k_) - sqr(sinTheta);
    scalar sqrP = (sqrt(sqr(r) + sqr(2*n_*k_)) + r)/2;
    scalar sqrQ = (sqrt(sqr(r) + sqr(2*n_*k_)) - r)/2;
    scalar p = sqrt(sqrP);

    scalar RN = (sqr(p - cosTheta) + sqrQ)/(sqr(p + cosTheta) + sqrQ);
    scalar RP = RN*(sqr(p - sinTheta*tanTheta) + sqrQ)/(sqr(p + sinTheta*tanTheta) + sqrQ);

    return (RP + RN)/2;
}


// ************************************************************************* //
