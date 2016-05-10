/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------

\*---------------------------------------------------------------------------*/

#include "extrapolatedSlipGradientFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

extrapolatedSlipGradientFvPatchVectorField::
extrapolatedSlipGradientFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF)
{}


extrapolatedSlipGradientFvPatchVectorField::
extrapolatedSlipGradientFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


extrapolatedSlipGradientFvPatchVectorField::
extrapolatedSlipGradientFvPatchVectorField
(
    const extrapolatedSlipGradientFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(ptf, p, iF, mapper)
{}


extrapolatedSlipGradientFvPatchVectorField::
extrapolatedSlipGradientFvPatchVectorField
(
    const extrapolatedSlipGradientFvPatchVectorField& ptf
)
:
    fixedGradientFvPatchVectorField(ptf)
{}


extrapolatedSlipGradientFvPatchVectorField::
extrapolatedSlipGradientFvPatchVectorField
(
    const extrapolatedSlipGradientFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void extrapolatedSlipGradientFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    word field = this->dimensionedInternalField().name();
    
    vectorField nf = patch().nf();
    
    volTensorField gradField 
        = fvc::grad(db().lookupObject<volVectorField>(field));

    gradient() 
        = (nf & gradField.boundaryField()[patch().index()]
        .patchInternalField()());
    
    fixedGradientFvPatchVectorField::updateCoeffs();
}

void extrapolatedSlipGradientFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const symmTensorField tang = I - sqr(patch().nf());
    vectorField::operator=
    (
        tang & (this->patchInternalField() + gradient()/this->patch().deltaCoeffs())
    );

    fvPatchVectorField::evaluate();
}

void extrapolatedSlipGradientFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    extrapolatedSlipGradientFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
