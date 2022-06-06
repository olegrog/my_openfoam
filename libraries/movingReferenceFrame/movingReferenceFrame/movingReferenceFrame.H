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
    Foam::movingReferenceFrame

Description
    Abstract base class for an inertial moving reference frame.

SourceFiles
    movingReferenceFrame.C

\*---------------------------------------------------------------------------*/

#ifndef movingReferenceFrame_H
#define movingReferenceFrame_H

#include "runTimeSelectionTables.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class movingReferenceFrame Declaration
\*---------------------------------------------------------------------------*/

class movingReferenceFrame
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Dictionary for serialization
        IOdictionary UrelDict_;

        //- Frame velocity
        dimensionedVector Urel_;


private:

    // Private data

        //- Flux due to reference frame
        surfaceScalarField phiRel_;

    // Private Member Functions

        //- Evaluate the frame velocity
        virtual void evaluate() = 0;


public:

    //- Runtime type information
    TypeName("movingReferenceFrame");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        movingReferenceFrame,
        mesh,
        (const fvMesh& mesh),
        (mesh)
    );

    //- Constructor
    movingReferenceFrame(const fvMesh& mesh);

    //- Selector
    static autoPtr<movingReferenceFrame> New(const fvMesh& mesh);

    //- Destructor
    virtual ~movingReferenceFrame() = default;


    // Member functions

        //- Access

            //- Return the relative velocity
            const dimensionedVector& Urel() const
            {
                return Urel_;
            }

            //- Return the relative fluxes
            const surfaceScalarField& phiRel() const
            {
                return phiRel_;
            }

            //- Is the reference frame moving?
            virtual bool moving() const
            {
                return true;
            }

        //- Evaluation

            //- Correct the relative fluxes and others
            void correct(bool oriented = true);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
