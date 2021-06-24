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
    Foam::phaseBoundary::tabulated

Description
    Tabulated phase boundary.

\*---------------------------------------------------------------------------*/

#ifndef tabulated_H
#define tabulated_H

#include "TableBase.H"

#include "alloyComponent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseBoundary
{

/*---------------------------------------------------------------------------*\
                         Class tabulated Declaration
\*---------------------------------------------------------------------------*/

class tabulated
:
    public Function1Types::TableBase<scalar>
{

public:

    //- Constructor
    tabulated
    (
        const dictionary& dict,
        const alloyComponent& component
    );


    // Member Functions

        //- Return the equilibrium concentration
        scalar equilibrium(scalar T) const
        {
            return value(T);
        }

        //- Return the slope along the solidification path
        scalar slope(scalar T) const
        {
            NotImplemented;
            return value(T);
        }

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const
        {
            NotImplemented;
            return tmp<Function1<scalar>>(new tabulated(*this));
        }
};

} // End namespace phaseBoundary
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
