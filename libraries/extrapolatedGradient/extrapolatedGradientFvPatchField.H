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

Class
    Foam::extrapolatedGradientFvPatchField

Description
    This boundary condition means that the second normal derivative is equal to
    zero.

    Credit to Eugene de Villiers for extrapolatedGradientFvPatchScalarField,
    which can be downloaded from
    https://www.cfd-online.com/Forums/openfoam/89731-second-derivative-zero-boundary-condition.html

SourceFiles
    extrapolatedGradientFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef extrapolatedGradientFvPatchField_H
#define extrapolatedGradientFvPatchField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class extrapolatedGradientFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class extrapolatedGradientFvPatchField
:
    public fixedGradientFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("extrapolatedGradient");


    // Constructors

        //- Construct from patch and internal field
        extrapolatedGradientFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        extrapolatedGradientFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  extrapolatedGradientFvPatchField<Type> onto a new patch
        extrapolatedGradientFvPatchField
        (
            const extrapolatedGradientFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        extrapolatedGradientFvPatchField
        (
            const extrapolatedGradientFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new extrapolatedGradientFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        extrapolatedGradientFvPatchField
        (
            const extrapolatedGradientFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new extrapolatedGradientFvPatchField<Type>(*this, iF)
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

#ifdef NoRepository
    #include "extrapolatedGradientFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
