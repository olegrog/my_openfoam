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

public:

    //- Runtime type information
    TypeName("erf");

    // Member Functions

        //- Apply the sigmoid function to a field
        virtual tmp<volScalarField> operator()(const volScalarField& field) const;

        //- Apply the first derivative of sigmoid function to a field
        virtual tmp<volScalarField> der(const volScalarField& field) const;
};

} // End namespace Foam

#endif
