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
    This file is part of slmMeltPoolFoam.

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
    Foam::reflection::noReflection

Description
    Dummy reflection model for 'none'

SourceFiles
    noReflection.C

\*---------------------------------------------------------------------------*/

#ifndef reflection_noReflection_H
#define reflection_noReflection_H

#include "reflectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reflection
{

/*---------------------------------------------------------------------------*\
                         Class noReflection Declaration
\*---------------------------------------------------------------------------*/

class noReflection
:
    public reflectionModel
{

public:

    //- Runtime type information
    TypeName("none");

    //- Constructor
    noReflection(const dictionary& dict) : reflectionModel(dict) {}

    //- Destructor
    ~noReflection() final = default;


    // Member Functions

        //- Return reflectivity
        scalar rho(const scalar incidentAngle) const final
        {
            return 0;
        }

        //- Return reflected direction
        vector R(const vector& incident, const vector& n) const final
        {
            return vector::zero;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reflection
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
