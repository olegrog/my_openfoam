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
    This file is part of movingReferenceFrame.

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
    Foam::givenFrame

Description
    Moving reference frame with a given velocity as a function of time
    using Function1.

SourceFiles
    givenFrame.C

\*---------------------------------------------------------------------------*/

#ifndef givenFrame_H
#define givenFrame_H

#include "Function1.H"

#include "movingReferenceFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class givenFrame Declaration
\*---------------------------------------------------------------------------*/

class givenFrame
:
    public movingReferenceFrame
{
    // Private data

        //- Dictionary
        const dictionary dict_;

        //- Frame velocity as a function of time
        const autoPtr<Function1<vector>> UrelFunc_;


    // Private Member Functions

        //- Evaluate the frame velocity
        void evaluate() final;

public:

    //- Runtime type information
    TypeName("given");

    //- Constructor
    givenFrame(const fvMesh& mesh);

    //- Destructor
    ~givenFrame() final = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
