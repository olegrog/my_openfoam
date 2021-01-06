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

\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

#include "Sigmoid.H"

#include "tanh.H"
#include "erf.H"
#include "cut.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSigmoidFunction(sigmoidType, functionType)                         \
                                                                               \
namespace Foam                                                                 \
{                                                                              \
    typedef sigmoidType<sigmoid::functionType> sigmoidType##functionType;      \
                                                                               \
    addNamedToRunTimeSelectionTable                                            \
    (                                                                          \
        sigmoidFunction,                                                       \
        sigmoidType##functionType,                                             \
        interval,                                                              \
        functionType                                                           \
    );                                                                         \
                                                                               \
    defineTemplateTypeNameWithName                                             \
    (                                                                          \
        sigmoidType##functionType,                                             \
        #sigmoidType"<"#functionType">"                                        \
    );                                                                         \
}

makeSigmoidFunction(Sigmoid, tanh)
makeSigmoidFunction(Sigmoid, erf)
makeSigmoidFunction(Sigmoid, cut)

// ************************************************************************* //
