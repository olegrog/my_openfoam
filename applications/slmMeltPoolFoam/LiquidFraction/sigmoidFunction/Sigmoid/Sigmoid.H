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
    Foam::sigmoid

Description
    Class template for arbitrary sigmoid function with internal checks.

\*---------------------------------------------------------------------------*/

#ifndef Sigmoid_H
#define Sigmoid_H

#include "sigmoidFunction.H"

namespace Foam
{

template<class Function>
class Sigmoid
:
    public sigmoidFunction
{

    // Private Member Functions

        //- Apply the normalized sigmoid function to a field (range between -1 and 1)
        tmp<volScalarField> value0(const tmp<volScalarField>& tfield) const override;

        //- Apply the derivative of the normalized sigmoid function to a field
        tmp<volScalarField> value1(const tmp<volScalarField>& tfield) const override;

        //- Auxiliary function for an arbitrary derivative
        template<label nPrimes>
        tmp<volScalarField> value(const tmp<volScalarField>& tfield) const;

public:

    //- Runtime type information
    TypeName("Sigmoid");

    //- Constructor
    Sigmoid(scalar a, scalar b);
};


} // End namespace Foam

#ifdef NoRepository
#   include "Sigmoid.C"
#endif

#endif
