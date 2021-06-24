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

#include "piecewiseLinear.H"

#include "error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseBoundary::piecewiseLinear::piecewiseLinear
(
    const dictionary& dict,
    const alloyComponent& component
)
:
    intervals_(dict.lookup("intervals")),
    liquidus_(component.alloy().liquidus().value())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::phaseBoundary::piecewiseLinear::equilibrium(scalar T) const
{
    for (const auto& interval : intervals_)
    {
        if (interval.Tmin < T && T <= interval.Tmax)
        {
            return interval.equilibrium + (T - liquidus_)/interval.slope;
        }
    }

    FatalErrorInFunction
        << "Temperature " << T << " is out of prescribed intervals." << endl
        << abort(FatalError);

    return 0;
}


Foam::scalar Foam::phaseBoundary::piecewiseLinear::slope(scalar T) const
{
    for (const auto& interval : intervals_)
    {
        if (interval.Tmin < T && T <= interval.Tmax)
        {
            return interval.slope;
        }
    }

    FatalErrorInFunction
        << "Temperature " << T << " is out of prescribed intervals." << endl
        << abort(FatalError);

    return 0;
}


// ************************************************************************* //
