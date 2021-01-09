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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "sigmoidFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(sigmoidFunction);
    defineRunTimeSelectionTable(sigmoidFunction, interval);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::sigmoidFunction> Foam::sigmoidFunction::New
(
    const dictionary& dict,
    scalar a,
    scalar b
)
{
    const word sigmoidType(dict.get<word>("sigmoid"));

    Info<< "Selecting sigmoid function " << sigmoidType << endl;

    const auto cstrIter = intervalConstructorTablePtr_->cfind(sigmoidType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "sigmoid",
            sigmoidType,
            *intervalConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<sigmoidFunction>(cstrIter()(a, b));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::sigmoidFunction::value
(
    const volScalarField& field
) const
{
    return (b_ - a_)/2*value0(2*field/(b_ - a_)) + (a_ + b_)/2;
}


Foam::tmp<Foam::volScalarField> Foam::sigmoidFunction::derivative
(
    const volScalarField& field
) const
{
    return value1(2*field/(b_ - a_));
}


// ************************************************************************* //
