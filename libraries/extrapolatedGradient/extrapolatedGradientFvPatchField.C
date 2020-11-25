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

\*---------------------------------------------------------------------------*/

#include "extrapolatedGradientFvPatchField.H"
//#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
//#include "volFields.H"
//#include "uniformDimensionedFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::extrapolatedGradientFvPatchField<Type>::
extrapolatedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::extrapolatedGradientFvPatchField<Type>::
extrapolatedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<Type>(p, iF)
{
    fvPatchField<Type>::operator=(this->patchInternalField());
    this->gradient() = Zero;
}


template<class Type>
Foam::extrapolatedGradientFvPatchField<Type>::
extrapolatedGradientFvPatchField
(
    const extrapolatedGradientFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::extrapolatedGradientFvPatchField<Type>::
extrapolatedGradientFvPatchField
(
    const extrapolatedGradientFvPatchField<Type>& ptf
)
:
    fixedGradientFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::extrapolatedGradientFvPatchField<Type>::
extrapolatedGradientFvPatchField
(
    const extrapolatedGradientFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::extrapolatedGradientFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    word fieldName = this->internalField().name();

    vectorField nf = this->patch().nf();

    typedef GeometricField<Type, fvPatchField, volMesh> volField;
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> volGradField;

    volGradField gradField = fvc::grad(this->db().template lookupObject<volField>(fieldName));

    this->gradient() =
    (
        nf & gradField.boundaryField()[this->patch().index()].patchInternalField()()
    );

    fixedGradientFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::extrapolatedGradientFvPatchField<Type>::write(Ostream& os) const
{
    fixedGradientFvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
