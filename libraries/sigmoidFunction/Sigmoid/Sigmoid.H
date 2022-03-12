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
    This file is part of sigmoidFunction.

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
    Foam::Sigmoid

Description
    Class template for an arbitrary sigmoid function with internal checks.

SourceFiles
    SigmoidTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Sigmoid_H
#define Sigmoid_H

#include "sigmoidFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Sigmoid Declaration
\*---------------------------------------------------------------------------*/

template<class Function>
class Sigmoid
:
    public sigmoidFunction
{

    // Private Member Functions

        //- Apply the normalised sigmoid function to a field (range between -1 and 1)
        tmp<volScalarField> value0(const tmp<volScalarField>& tfield) const final
        {
            return value<0>(tfield);
        }

        //- Apply the derivative of the normalised sigmoid function to a field
        tmp<volScalarField> value1(const tmp<volScalarField>& tfield) const final
        {
            return value<1>(tfield);
        }

        //- Auxiliary function for an arbitrary derivative
        template<label nPrimes>
        tmp<volScalarField> value(const tmp<volScalarField>& tfield) const;

public:

    //- Runtime type information
    TypeName("Sigmoid");

    //- Constructor
    Sigmoid(scalar a, scalar b);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SigmoidTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
