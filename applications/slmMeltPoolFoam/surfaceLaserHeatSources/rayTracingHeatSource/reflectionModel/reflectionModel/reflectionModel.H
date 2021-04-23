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
    Foam::reflectionModel

Description
    Abstract base class for radiation reflection and refraction.

\*---------------------------------------------------------------------------*/

#ifndef reflectionModel_H
#define reflectionModel_H

#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class reflectionModel Declaration
\*---------------------------------------------------------------------------*/

class reflectionModel
{
public:

    //- Runtime type information
    TypeName("reflectionModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        reflectionModel,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    //- Constructor
    reflectionModel(const dictionary& dict) {}

    //- Selector
    static autoPtr<reflectionModel> New(const dictionary& dict);

    //- Destructor
    virtual ~reflectionModel() = default;


    // Member Functions

        //- Return reflected direction
        virtual vector R(const vector& incident, const vector& n) const = 0;

        //- Return reflectivity from medium1 to medium2 and a incident angle.
        virtual scalar rho(const scalar incidentAngle) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
