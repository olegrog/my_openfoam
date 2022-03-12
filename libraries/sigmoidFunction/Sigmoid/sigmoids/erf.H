/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of sigmoidFunction.

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

Class
    Foam::sigmoid::erf

Description
    Error function.

\*---------------------------------------------------------------------------*/

#ifndef sigmoid_erf_H
#define sigmoid_erf_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sigmoid
{

/*---------------------------------------------------------------------------*\
                             Class erf Declaration
\*---------------------------------------------------------------------------*/

struct erf
{
    template<int nPrimes>
    static scalar value(scalar x);
};

template<>
scalar erf::value<0>(scalar x)
{
    using constant::mathematical::pi;
    return Foam::erf(sqrt(pi)/2*x);
}

template<>
scalar erf::value<1>(scalar x)
{
    using constant::mathematical::pi;
    return exp(-sqr(sqrt(pi)/2*x));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sigmoid
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
