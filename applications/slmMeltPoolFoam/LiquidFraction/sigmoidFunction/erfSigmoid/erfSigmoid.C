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

#include "erfSigmoid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeName(erfSigmoid);

    addToRunTimeSelectionTable
    (
        sigmoidFunction,
        erfSigmoid,
        interval
    );
}

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::erfSigmoid::value
(
    const tmp<volScalarField>& tfield
) const
{
    using constant::mathematical::pi;
    return erf(sqrt(pi)/2*tfield);
}

Foam::tmp<Foam::volScalarField> Foam::erfSigmoid::derivative1
(
    const tmp<volScalarField>& tfield
) const
{
    using constant::mathematical::pi;
    return exp(-sqr(sqrt(pi)/2*tfield));
};

// ************************************************************************* //
