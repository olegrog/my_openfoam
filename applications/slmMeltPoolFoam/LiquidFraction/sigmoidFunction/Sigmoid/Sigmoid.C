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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "sigFpe.H"

#include "Sigmoid.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Function>
Foam::Sigmoid<Function>::Sigmoid(scalar a, scalar b)
:
    sigmoidFunction(a, b)
{
    const bool sigActive = sigFpe::active();
    // SIGFPE can be trapped, but here we have to disable this feature
    if (sigActive) sigFpe::unset();
    const List<Pair<scalar>> actualVsExpected
    {
        { Function::template value<0>(-VGREAT), -1},
        { Function::template value<0>(0), 0},
        { Function::template value<0>(VGREAT), 1},
        { Function::template value<1>(-VGREAT), 0},
        { Function::template value<1>(0), 1},
        { Function::template value<1>(VGREAT), 0},
    };
    if (sigActive) sigFpe::set();

    for (const auto& pair : actualVsExpected)
    {
        if (notEqual(pair.first(), pair.second()))
        {
            FatalErrorInFunction
                << "wrong implemented sigmoid function since "
                << pair.first() << " != " << pair.second()
                << exit(FatalError);
        }
    }
}

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class Function>
template<int nPrimes>
Foam::tmp<Foam::volScalarField> Foam::Sigmoid<Function>::value
(
    const tmp<volScalarField>& tfield
) const
{
    const volScalarField& field = tfield();

    // This constructor allocates storage for the field but does not set values.
    auto tres =
        tmp<volScalarField>::New
        (
            IOobject
            (
                "sigmoid" + string(nPrimes, '\'') + '(' + field.name() + ')',
                field.instance(),
                field.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            field.mesh(),
            field.dimensions()
        );

    volScalarField& res = tres.ref();

    forAll(field, cellI)
    {
        res[cellI] = Function::template value<nPrimes>(field[cellI]);
    }
    forAll(field.boundaryField(), patchi)
    {
        forAll(field.boundaryField()[patchi], faceI)
        {
            res.boundaryFieldRef(false)[patchi][faceI] = Function::template value<nPrimes>
                (field.boundaryField()[patchi][faceI]);
        }
    }

    return tres;
}


template<class Function>
Foam::tmp<Foam::volScalarField> Foam::Sigmoid<Function>::value0
(
    const tmp<volScalarField>& tfield
) const
{
    return value<0>(tfield);
}

template<class Function>
Foam::tmp<Foam::volScalarField> Foam::Sigmoid<Function>::value1
(
    const tmp<volScalarField>& tfield
) const
{
    return value<1>(tfield);
}


// ************************************************************************* //
