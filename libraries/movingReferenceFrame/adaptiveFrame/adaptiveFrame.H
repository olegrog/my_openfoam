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
    Foam::adaptiveFrame

Description
    Moving reference frame with an adaptive velocity.

SourceFiles
    adaptiveFrame.C

\*---------------------------------------------------------------------------*/

#ifndef adaptiveFrame_H
#define adaptiveFrame_H

#include "movingReferenceFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class adaptiveFrame Declaration
\*---------------------------------------------------------------------------*/

class adaptiveFrame
:
    public movingReferenceFrame
{
    // Private data

        //- Dictionary
        const dictionary dict_;

        //- Direction of the frame velocity
        const vector direction_;

        //- Adaptation field name
        const word fieldName_;

        //- Reference to the adaptation field
        const volScalarField& field_;

        //- Target mean value of the field
        const scalar meanValue_;

        //- Mean value of the field at the previous time step
        scalar meanValuePrev_;

        //- Factor for deviation from the mean value
        const scalar valueFactor_;

        //- Factor for deviation from its derivative
        const scalar derivativeFactor_;


    // Private Member Functions

        //- Evaluate the frame velocity
        void evaluate() final;

public:

    //- Runtime type information
    TypeName("adaptive");

    //- Constructor
    adaptiveFrame(const fvMesh& mesh);

    //- Destructor
    ~adaptiveFrame() final = default;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
