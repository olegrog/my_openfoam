/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

#include "snitPressureFvPatchField.H"
#include "volFields.H"                      // for volScalarField
#include "addToRunTimeSelectionTable.H"     // for makePatchTypeField

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

snitPressureFvPatchField::snitPressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


snitPressureFvPatchField::snitPressureFvPatchField
(
    const snitPressureFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


snitPressureFvPatchField::snitPressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF)
{
    // initialize as zero gradient, because transportProperties are inaccessible
    operator==(patchInternalField());
}


snitPressureFvPatchField::snitPressureFvPatchField
(
    const snitPressureFvPatchField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf)
{}


snitPressureFvPatchField::snitPressureFvPatchField
(
    const snitPressureFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void snitPressureFvPatchField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const dictionary& transportProperties =
        db().lookupObject<IOdictionary>("transportProperties");
    const scalar gamma7 =
        dimensionedScalar("gamma7", transportProperties).value();
    const fvPatchField<scalar>& T =
        patch().lookupPatchField<volScalarField, scalar>("T");

    gradient() = -0.25 * gamma7 / T * pow3(T.snGrad());
    fixedGradientFvPatchScalarField::updateCoeffs();
}

void snitPressureFvPatchField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
}

makePatchTypeField
(
    fvPatchScalarField,
    snitPressureFvPatchField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
