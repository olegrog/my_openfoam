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

#include "error.H"
#include "sigFpe.H"

#include "generateGeometricField.H"

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
                << "Wrong implemented sigmoid function since "
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

    return generateGeometricField<volScalarField>
    (
        "sigmoid" + string(nPrimes, '\'') + '(' + field.name() + ')',
        field.mesh(),
        dimless,
        [](scalar field)
        {
            return Function::template value<nPrimes>(field);
        },
        field
    );
}


// ************************************************************************* //
