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
    This file is part of setZNDsolution.

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

#include "ZNDsystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ZNDsystem::ZNDsystem(const dictionary& dict)
:
    Q_(dict.get<scalar>("Q")),
    E_(dict.get<scalar>("E")),
    k_(dict.get<scalar>("k")),
    g_(dict.get<scalar>("g")),
    D_(sqrt(g_ + (sqr(g_) - 1)/2*Q_) + sqrt((sqr(g_) - 1)/2*Q_))
{
    Info<< "Chapman--Jouguet velocity = " << D_ << endl;
}

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

Foam::scalar _1l(Foam::scalar lambda)
{
    return Foam::max(1 - lambda, Foam::SMALL);
};

} // End namespace

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::ZNDsystem::U(Foam::scalar lambda) const
{
    return (sqr(D_) - g_)/D_/(g_ + 1)*(1 + sqrt(_1l(lambda)));
}


Foam::scalar Foam::ZNDsystem::p(Foam::scalar lambda) const
{
    return (1 + sqr(D_))/(g_ + 1)*(1 + (sqr(D_) - g_)/(1 + sqr(D_))*sqrt(_1l(lambda)));
}


Foam::scalar Foam::ZNDsystem::rho(Foam::scalar lambda) const
{
    return sqr(D_)*(g_ + 1)/(g_*(1 + sqr(D_)) - (sqr(D_) - g_)*sqrt(_1l(lambda)));
}


Foam::label Foam::ZNDsystem::nEqns() const
{
    return 1;
}


void Foam::ZNDsystem::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    dydx[0] = -k_*_1l(y[0])/(U(y[0]) - D_)*exp(-rho(y[0])/p(y[0])*E_);
}


// ************************************************************************* //
