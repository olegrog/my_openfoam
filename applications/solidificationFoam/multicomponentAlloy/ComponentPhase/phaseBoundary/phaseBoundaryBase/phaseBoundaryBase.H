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
    Foam::phaseBoundaryBase

Description
    Phase boundary base class.

\*---------------------------------------------------------------------------*/

#ifndef phaseBoundaryBase_H
#define phaseBoundaryBase_H

#include "multicomponentAlloy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class phaseBoundaryBase Declaration
\*---------------------------------------------------------------------------*/

class phaseBoundaryBase
{
protected:

    // Protected Data

        //- Reference to the component phase
        const componentPhase& phase_;

        //- Reference to the alloy component
        const alloyComponent& component_;

        //- Reference to the multicomponent alloy
        const multicomponentAlloy& alloy_;

        //- The liquidus temperature
        const scalar liquidus_;

public:

    //- Constructor
    phaseBoundaryBase
    (
        const dictionary& dict,
        const componentPhase& phase
    )
    :
        phase_(phase),
        component_(phase.component()),
        alloy_(component_.alloy()),
        liquidus_(alloy_.liquidus().value())
    {}
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
