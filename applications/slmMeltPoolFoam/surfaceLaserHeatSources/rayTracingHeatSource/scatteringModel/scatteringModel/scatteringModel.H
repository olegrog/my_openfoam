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
    Foam::scatteringModel

Description
    Abstract base class for radiation scattering.

\*---------------------------------------------------------------------------*/

#ifndef scatteringModel_H
#define scatteringModel_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "scalar.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class scatteringModel Declaration
\*---------------------------------------------------------------------------*/

class scatteringModel
{
    // Private data

        // Scattering threshold
        const scalar threshold_;

public:

    //- Runtime type information
    TypeName("scatteringModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        scatteringModel,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    //- Constructor
    scatteringModel(const dictionary& dict);

    //- Selector
    static autoPtr<scatteringModel> New(const dictionary& dict);

    //- Destructor
    virtual ~scatteringModel() = default;


    // Member Functions

        //- Return scattering threshold
        scalar threshold() const
        {
            return threshold_;
        }

        //- Return direction of reflection
        virtual vector reflection(const vector& i, const vector& n) const = 0;

        //- Return direction of refraction
        virtual vector refraction(const vector& i, const vector& n) const = 0;

        //- Return reflectivity
        virtual scalar R(scalar cosTheta) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
