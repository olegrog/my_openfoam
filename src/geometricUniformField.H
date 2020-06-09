/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is a missed part of OpenFOAM.

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
    Foam::geometricUniformField

Description
    A class representing the concept of a GeometricField of a uniform Type
    used to avoid unnecessary manipulations for objects known at compile-time.

\*---------------------------------------------------------------------------*/

#ifndef geometricUniformField_H
#define geometricUniformField_H

#include "dimensionedTypes.H"
#include "UniformField.H"

#include "uniformFieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class geometricUniformField
:
    public UniformField<Type>
{

    // Private data

        dimensionSet dimensions_;

public:

    // Public typedefs

        typedef UniformField<Type> Internal;
        typedef UniformField<Type> Patch;
        typedef uniformFieldField<Type> Boundary;
        typedef Type cmptType;


    // Constructors

        //- Constructor
        geometricUniformField(const dimensioned<Type>& quantity)
        :
            UniformField<Type>(quantity.value()),
            dimensions_(quantity.dimensions())
        {}


    // Member Operators

        const dimensionSet& dimensions() const
        {
            return dimensions_;
        }

        UniformField<Type> operator()() const
        {
            return *this;
        }

        UniformField<Type> oldTime() const
        {
            return *this;
        }

        uniformFieldField<Type> boundaryField() const
        {
            return uniformFieldField<Type>(this->field());
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
