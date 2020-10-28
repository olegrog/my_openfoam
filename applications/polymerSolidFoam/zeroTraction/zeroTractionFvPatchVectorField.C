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
    This file is part of polymerSolidFoam.

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

#include "zeroTractionFvPatchVectorField.H"

#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zeroTractionFvPatchVectorField::zeroTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = Zero;
}


Foam::zeroTractionFvPatchVectorField::zeroTractionFvPatchVectorField
(
    const zeroTractionFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper)
{}


Foam::zeroTractionFvPatchVectorField::zeroTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vectorField("gradient", dict, p.size());
}


Foam::zeroTractionFvPatchVectorField::zeroTractionFvPatchVectorField
(
    const zeroTractionFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf)
{}


Foam::zeroTractionFvPatchVectorField::zeroTractionFvPatchVectorField
(
    const zeroTractionFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zeroTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& E = patch().lookupPatchField<volScalarField, scalar>("E");
    const fvPatchField<scalar>& nu = patch().lookupPatchField<volScalarField, scalar>("nu");

    scalarField mu(E/(2.0*(1.0 + nu)));
    scalarField lambda(nu*E/((1.0 + nu)*(1.0 - 2.0*nu)));

    scalarField twoMuLambda(2*mu + lambda);

    vectorField n(patch().nf());

    const fvPatchField<symmTensor>& sigma =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigma");

    gradient() = snGrad() - (n & sigma)/twoMuLambda;

    fixedGradientFvPatchVectorField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

makePatchTypeField
(
    fvPatchVectorField,
    zeroTractionFvPatchVectorField
);

} // End namespace Foam


// ************************************************************************* //
