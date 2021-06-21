/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020-2021 Oleg Rogozin
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
    Foam::sigmoidFunction

Description
    Abstract base class for sigmoid function.

\*---------------------------------------------------------------------------*/

#ifndef sigmoidFunction_H
#define sigmoidFunction_H

#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class sigmoidFunction Declaration
\*---------------------------------------------------------------------------*/

class sigmoidFunction
{
    // Private data

        //- Minimum value of sigmoid
        const scalar a_;

        //- Maximum value of sigmoid
        const scalar b_;

    // Private Member Functions

        //- Apply the normalised sigmoid function to a field (range between -1 and 1)
        virtual tmp<volScalarField> value0(const tmp<volScalarField>& tfield) const = 0;

        //- Apply the derivative of the normalised sigmoid function to a field
        virtual tmp<volScalarField> value1(const tmp<volScalarField>& tfield) const = 0;

public:

    //- Runtime type information
    TypeName("sigmoidFunction");

    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        sigmoidFunction,
        interval,
        (scalar a, scalar b),
        (a, b)
    );

    //- Constructor
    sigmoidFunction(scalar a, scalar b) : a_(a), b_(b) {}

    //- Selector
    static autoPtr<sigmoidFunction> New(const dictionary& dict, scalar a, scalar b);

    //- Destructor
    virtual ~sigmoidFunction() = default;

    // Member Functions

        //- Apply the sigmoid function to a field
        tmp<volScalarField> value(const volScalarField& field) const;

        //- Apply the derivative of the sigmoid function to a field
        tmp<volScalarField> derivative(const volScalarField& field) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
