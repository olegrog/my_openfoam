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
    Foam::phaseBoundary::linearised

Description
    Linearised phase boundary.

\*---------------------------------------------------------------------------*/

#ifndef linearised_H
#define linearised_H

#include "phaseBoundaryBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseBoundary
{

/*---------------------------------------------------------------------------*\
                         Class linearised Declaration
\*---------------------------------------------------------------------------*/

class linearised
:
    public phaseBoundaryBase
{
    // Private Data

        //- Equilibrium concentration
        const scalar equilibrium_;

        //- Slope along the solidification path
        const scalar slope_;

public:

    //- Constructor
    linearised
    (
        const dictionary& dict,
        const componentPhase& phase
    )
    :
        phaseBoundaryBase(dict, phase),
        equilibrium_(dict.get<scalar>("equilibrium")),
        slope_(dict.get<scalar>("slope"))
    {}


    // Member Functions

        //- Return the equilibrium concentration
        scalar equilibrium(scalar T) const
        {
            return equilibrium_ + (T - liquidus_)/slope_;
        }

        //- Return the slope along the solidification path
        scalar slope(scalar) const
        {
            return slope_;
        }
};

} // End namespace phaseBoundary
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
