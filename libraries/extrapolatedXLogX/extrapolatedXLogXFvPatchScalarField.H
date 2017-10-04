/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::extrapolatedXLogXFvPatchScalarField

Description

SourceFiles
    extrapolatedXLogXFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef extrapolatedXLogXFvPatchScalarFields_H
#define extrapolatedXLogXFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class extrapolatedXLogXFvPatch Declaration
\*---------------------------------------------------------------------------*/

class extrapolatedXLogXFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

public:

    //- Runtime type information
    TypeName("extrapolatedXLogX");


    // Constructors

        //- Construct from patch and internal field
        extrapolatedXLogXFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        extrapolatedXLogXFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  extrapolatedXLogXFvPatchScalarField onto a new patch
        extrapolatedXLogXFvPatchScalarField
        (
            const extrapolatedXLogXFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        extrapolatedXLogXFvPatchScalarField
        (
            const extrapolatedXLogXFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new extrapolatedXLogXFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        extrapolatedXLogXFvPatchScalarField
        (
            const extrapolatedXLogXFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new extrapolatedXLogXFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
