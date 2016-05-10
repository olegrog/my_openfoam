/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------

\*---------------------------------------------------------------------------*/

#include "extrapolatedGradientFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

extrapolatedGradientFvPatchScalarField::
extrapolatedGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


extrapolatedGradientFvPatchScalarField::
extrapolatedGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF)
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}


extrapolatedGradientFvPatchScalarField::
extrapolatedGradientFvPatchScalarField
(
    const extrapolatedGradientFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


extrapolatedGradientFvPatchScalarField::
extrapolatedGradientFvPatchScalarField
(
    const extrapolatedGradientFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf)
{}


extrapolatedGradientFvPatchScalarField::
extrapolatedGradientFvPatchScalarField
(
    const extrapolatedGradientFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void extrapolatedGradientFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    word field = this->dimensionedInternalField().name();
    
    vectorField nf = patch().nf();
    
    volVectorField gradField 
        = fvc::grad(db().lookupObject<volScalarField>(field));

    gradient() 
        = (nf & gradField.boundaryField()[patch().index()]
        .patchInternalField()());
    
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void extrapolatedGradientFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    extrapolatedGradientFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
