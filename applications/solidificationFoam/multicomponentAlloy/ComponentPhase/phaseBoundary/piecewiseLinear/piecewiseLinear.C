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

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void Foam::phaseBoundary::piecewiseLinear::outOfBounds(scalar T) const
{
    FatalErrorInFunction
        << "Temperature " << T << " is out of prescribed intervals for "
        << phase_.keyword() << " " << component_.keyword() << endl
        << exit(FatalError);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseBoundary::piecewiseLinear::piecewiseLinear
(
    const dictionary& dict,
    const componentPhase& phase
)
:
    phaseBoundaryBase(dict, phase),
    intervals_(dict.lookup("intervals"))
{
    // --- Checks

    const scalar Cthreshold = 1e-3;
    scalar Tprev = 0;
    scalar Cprev = 0;

    for (const auto& interval : intervals_)
    {
        if (Tprev > 0)
        {
            if (mag(Tprev - interval.Tmin) > SMALL)
            {
                FatalErrorInFunction
                    << "Intervals are unsorted: " << Tprev << " != " << interval.Tmin
                    << exit(FatalError);
            }

            if (mag(Cprev - C(interval, interval.Tmin)) > Cthreshold)
            {
                FatalErrorInFunction
                    << "Piecewise approximation has jump "
                    << C(interval, interval.Tmin) - Cprev << " at T = " << Tprev
                    << exit(FatalError);
            }
        }

        Tprev = interval.Tmax;
        Cprev = C(interval, Tprev);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::phaseBoundary::piecewiseLinear::equilibrium(scalar T) const
{
    for (const auto& interval : intervals_)
    {
        if (interval.includes(T))
        {
            return C(interval, T);
        }
    }

    outOfBounds(T);
    return 0;
}


Foam::scalar Foam::phaseBoundary::piecewiseLinear::slope(scalar T) const
{
    for (const auto& interval : intervals_)
    {
        if (interval.includes(T))
        {
            return interval.slope;
        }
    }

    outOfBounds(T);
    return 0;
}


// ************************************************************************* //
