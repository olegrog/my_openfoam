/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::extrapolatedSlipGradientFvPatchVectorField

Description

SourceFiles
    extrapolatedSlipGradientFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef extrapolatedSlipGradientFvPatchVectorFields_H
#define extrapolatedSlipGradientFvPatchVectorFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class extrapolatedGradientFvPatch Declaration
\*---------------------------------------------------------------------------*/

class extrapolatedSlipGradientFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{
    // Private data

public:

    //- Runtime type information
    TypeName("extrapolatedSlipGradient");


    // Constructors

        //- Construct from patch and internal field
        extrapolatedSlipGradientFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        extrapolatedSlipGradientFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  extrapolatedSlipGradientFvPatchVectorField onto a new patch
        extrapolatedSlipGradientFvPatchVectorField
        (
            const extrapolatedSlipGradientFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        extrapolatedSlipGradientFvPatchVectorField
        (
            const extrapolatedSlipGradientFvPatchVectorField&
        );

        //- Construct and return a clone
        tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new extrapolatedSlipGradientFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        extrapolatedSlipGradientFvPatchVectorField
        (
            const extrapolatedSlipGradientFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new extrapolatedSlipGradientFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        void updateCoeffs();

        //- Evaluate the patch field
        void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::blocking
        );

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
