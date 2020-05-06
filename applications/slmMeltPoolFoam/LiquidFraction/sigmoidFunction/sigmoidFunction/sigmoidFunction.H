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
    Foam::sigmoidFunction

Description
    Abstract base class for sigmoid function.

\*---------------------------------------------------------------------------*/

#ifndef sigmoidFunction_H
#define sigmoidFunction_H

#include "volFields.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class sigmoidFunction
{

public:

    //- Runtime type information
    TypeName("sigmoidFunction");

    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        sigmoidFunction,
        empty,
        (),
        ()
    );

    //- Constructor
    sigmoidFunction() = default;

    //- Selector
    static autoPtr<sigmoidFunction> New(const dictionary& dict);

    //- Destructor
    virtual ~sigmoidFunction() = default;

    // Member Functions

        //- Apply the sigmoid function to a field
        virtual tmp<volScalarField> operator()(const volScalarField& field) const = 0;

        //- Apply the first derivative of sigmoid function to a field
        virtual tmp<volScalarField> der(const volScalarField& field) const = 0;
};

} // End namespace Foam

#endif
