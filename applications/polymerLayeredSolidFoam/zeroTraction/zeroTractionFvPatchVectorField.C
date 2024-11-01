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

    const vectorField n(patch().nf());
    const fvPatchField<scalar>& mu = patch().lookupPatchField<volScalarField, scalar>("mu");
    const fvPatchField<scalar>& lambda = patch().lookupPatchField<volScalarField, scalar>("lambda");
    const fvPatchField<symmTensor>& sigma =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigma");
    const fvPatchField<vector>& boundaryForce =
        patch().lookupPatchField<volVectorField, vector>("boundaryForce");

    gradient() = snGrad() - ((n & sigma) - boundaryForce)/(2*mu + lambda);

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void Foam::zeroTractionFvPatchVectorField::write(Ostream& os) const
{
    // write gradient for restarting
    fixedGradientFvPatchVectorField::write(os);
    // write value for more accurate vizualization
    writeEntry("value", os);
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
