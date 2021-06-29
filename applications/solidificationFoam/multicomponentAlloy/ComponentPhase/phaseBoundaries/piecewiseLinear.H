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

Class
    Foam::phaseBoundary::piecewiseLinear

Description
    Piecewise linear phase boundary.

\*---------------------------------------------------------------------------*/

#ifndef piecewiseLinear_H
#define piecewiseLinear_H

#include "multicomponentAlloy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseBoundary
{

/*---------------------------------------------------------------------------*\
                       Class piecewiseLinear Declaration
\*---------------------------------------------------------------------------*/

class piecewiseLinear
{
    //- Struct containing coeffs for a single temperature interval
    struct phaseBoundaryInterval
    {
        //- Minimum temperature
        const scalar Tmin;

        //- Maximum temperature
        const scalar Tmax;

        //- Equilibrium concentration
        const scalar equilibrium;

        //- Slope along the solidification path
        const scalar slope;

        //- Constructor from dictionary
        phaseBoundaryInterval(const dictionary& dict)
        :
            Tmin(dict.get<scalar>("Tmin")),
            Tmax(dict.get<scalar>("Tmax")),
            equilibrium(dict.get<scalar>("equilibrium")),
            slope(dict.get<scalar>("slope"))
        {}

        //- No empty construct
        phaseBoundaryInterval() = delete;

        //- Return clone (need for PtrList)
        autoPtr<phaseBoundaryInterval> clone() const
        {
            notImplemented("phaseBoundaryInterval::clone() const");
            return autoPtr<phaseBoundaryInterval>();
        }

        //- Selector
        static autoPtr<phaseBoundaryInterval> New(Istream& is)
        {
            dictionaryEntry entry(dictionary::null, is);
            return autoPtr<phaseBoundaryInterval>(new phaseBoundaryInterval(entry));
        }

        // Member Functions

            //- Check if temperature lies within the interval
            bool includes(scalar T) const
            {
                return Tmin < T && T < Tmax*(1 + SMALL);
            }
    };


    // Private Data

        //- List of temperature intervals with coeffs
        PtrList<phaseBoundaryInterval> intervals_;

        //- Reference to the alloy component
        const alloyComponent& component_;

        //- The liquidus temperature
        const scalar liquidus_;


    // Private Member Functions

        //- Return the equilibrium concentration within the interval
        scalar C(const phaseBoundaryInterval& interval, scalar T) const
        {
            return interval.equilibrium + (T - liquidus_)/interval.slope;
        }

        //- Abort if the temperature is out of specified intervals
        void outOfBounds(scalar T) const;

public:

    //- Constructor
    piecewiseLinear
    (
        const dictionary& dict,
        const alloyComponent& component
    );


    // Member Functions

        //- Return the equilibrium concentration
        scalar equilibrium(scalar T) const;

        //- Return the slope along the solidification path
        scalar slope(scalar T) const;
};

} // End namespace phaseBoundary
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
