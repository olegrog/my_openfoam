/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020-2020 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of thermalDebindingFoam.

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

#include "polymerComponent.H"

#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polymerComponent::polymerComponent
(
    const word& name,
    const dictionary& polymerComponentDict,
    const fvMesh& mesh
)
:
    volScalarField
    (
        IOobject
        (
            "massFraction" + name,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("initialValue", dimless, polymerComponentDict),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    name_(name),
    polymerComponentDict_(polymerComponentDict),
    rateConstant_("rateConstant", inv(dimTime), polymerComponentDict_),
    activationEnergy_("activationEnergy", dimEnergy/dimMoles, polymerComponentDict_)
{}


// ************************************************************************* //
