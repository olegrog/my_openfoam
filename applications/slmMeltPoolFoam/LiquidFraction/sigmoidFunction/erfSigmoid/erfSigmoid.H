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
    Foam::erfSigmoid

Description
    Error function.

\*---------------------------------------------------------------------------*/

#ifndef erfSigmoid_H
#define erfSigmoid_H

#include "sigmoidFunction.H"

namespace Foam
{

class erfSigmoid
:
    public sigmoidFunction
{

    // Private Member Functions

        //- Apply the normalized sigmoid function to a field (range between -1 and 1)
        tmp<volScalarField> value(const tmp<volScalarField>& tfield) const;

        //- Apply the first derivative of the normalized sigmoid function to a field
        tmp<volScalarField> derivative1(const tmp<volScalarField>& tfield) const;

public:

    //- Runtime type information
    TypeName("erf");

    //- Constructor
    erfSigmoid(scalar a, scalar b) : sigmoidFunction(a, b) {}
};

} // End namespace Foam

#endif
