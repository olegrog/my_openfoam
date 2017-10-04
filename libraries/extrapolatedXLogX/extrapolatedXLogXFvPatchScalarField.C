/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------

\*---------------------------------------------------------------------------*/

#define NoConstructFromTmp
#include "extrapolatedXLogXFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"
#undef NoConstructFromTmp

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

extrapolatedXLogXFvPatchScalarField::
extrapolatedXLogXFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


extrapolatedXLogXFvPatchScalarField::
extrapolatedXLogXFvPatchScalarField
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


extrapolatedXLogXFvPatchScalarField::
extrapolatedXLogXFvPatchScalarField
(
    const extrapolatedXLogXFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


extrapolatedXLogXFvPatchScalarField::
extrapolatedXLogXFvPatchScalarField
(
    const extrapolatedXLogXFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf)
{}


extrapolatedXLogXFvPatchScalarField::
extrapolatedXLogXFvPatchScalarField
(
    const extrapolatedXLogXFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void extrapolatedXLogXFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    word field = this->internalField().name();
    
    vectorField nf = patch().nf();
    
    volVectorField gradField 
        = fvc::grad(db().lookupObject<volScalarField>(field));

    extrapolatedXLogXFvPatchScalarField::operator==(
        this->internalField() - 2*(nf & gradField.boundaryField()[patch().index()]
                .patchInternalField()())
    );
    gradient() 
        = 2*(nf & gradField.boundaryField()[patch().index()]
        .patchInternalField()());
    
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void extrapolatedXLogXFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    extrapolatedXLogXFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
