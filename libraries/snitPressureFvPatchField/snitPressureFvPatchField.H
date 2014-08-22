/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
-------------------------------------------------------------------------------
Class
    Foam::snitPressureFvPatchField

Group
    grpWallBoundaryConditions

Description
    This boundary condition applies a calculating-gradient condition on the
    pressure for using in conjunction with the snitFoam solver.

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            snitPressureGradient;
    }
    \endverbatim

SourceFiles
    snitPressureFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef snitPressureFvPatchField_H
#define snitPressureFvPatchField_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class snitPressureFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class snitPressureFvPatchField
:
    public fixedGradientFvPatchScalarField
{
public:

    //- Runtime type information
    TypeName("snitPressureGradient");


// Constructors

    //- Construct from patch and internal field
    snitPressureFvPatchField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    snitPressureFvPatchField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping the given snitPressureFvPatchField
    //  onto a new patch
    snitPressureFvPatchField
    (
        const snitPressureFvPatchField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    snitPressureFvPatchField
    (
        const snitPressureFvPatchField&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField > clone() const
    {
        return tmp<fvPatchScalarField >
        (
            new snitPressureFvPatchField(*this)
        );
    }

    //- Construct as copy setting internal field reference
    snitPressureFvPatchField
    (
        const snitPressureFvPatchField&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField > clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField >
        (
            new snitPressureFvPatchField(*this, iF)
        );
    }


// Member functions

    // Evaluation functions

    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
