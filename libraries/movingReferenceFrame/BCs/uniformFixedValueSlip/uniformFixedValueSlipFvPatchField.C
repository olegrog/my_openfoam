/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2022 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of movingReferenceFrame.

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

#include "uniformFixedValueSlipFvPatchField.H"

#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformFixedValueSlipFvPatchField<Type>::uniformFixedValueSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(p, iF),
    uniformValue_(nullptr)
{}


template<class Type>
Foam::uniformFixedValueSlipFvPatchField<Type>::uniformFixedValueSlipFvPatchField
(
    const uniformFixedValueSlipFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchField<Type>(p, iF),   // Don't map
    uniformValue_(ptf.uniformValue_.clone(p.patch()))
{
    if (mapper.direct() && !mapper.hasUnmapped())
    {
        // Use mapping instead of re-evaluation
        this->map(ptf, mapper);
    }
    else
    {
        // Evaluate since value not mapped
        this->evaluate();
    }
}


template<class Type>
Foam::uniformFixedValueSlipFvPatchField<Type>::uniformFixedValueSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    transformFvPatchField<Type>(p, iF, dict),
    uniformValue_(PatchFunction1<Type>::New(p.patch(), "uniformValue", dict))
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        this->evaluate();
    }
}


template<class Type>
Foam::uniformFixedValueSlipFvPatchField<Type>::uniformFixedValueSlipFvPatchField
(
    const uniformFixedValueSlipFvPatchField<Type>& ptf
)
:
    transformFvPatchField<Type>(ptf),
    uniformValue_(ptf.uniformValue_.clone(this->patch().patch()))
{}


template<class Type>
Foam::uniformFixedValueSlipFvPatchField<Type>::uniformFixedValueSlipFvPatchField
(
    const uniformFixedValueSlipFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(ptf, iF),
    uniformValue_(ptf.uniformValue_.clone(this->patch().patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformFixedValueSlipFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    transformFvPatchField<Type>::autoMap(mapper);
    uniformValue_().autoMap(mapper);

    if (uniformValue_().constant())
    {
        // If mapper is not dependent on time we're ok to evaluate
        this->evaluate();
    }
}


template<class Type>
void Foam::uniformFixedValueSlipFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    transformFvPatchField<Type>::rmap(ptf, addr);

    const uniformFixedValueSlipFvPatchField& tiptf =
        refCast<const uniformFixedValueSlipFvPatchField>(ptf);

    uniformValue_().rmap(tiptf.uniformValue_(), addr);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::uniformFixedValueSlipFvPatchField<Type>::snGrad() const
{
    const vectorField nHat(this->patch().nf());
    const Field<Type> pif(this->patchInternalField());
    const scalar t = this->db().time().timeOutputValue();

    return
    (
        nHat*(nHat & uniformValue_->value(t))
      + transform(I - sqr(nHat), pif)
      - pif
    )*this->patch().deltaCoeffs();
}


template<class Type>
void Foam::uniformFixedValueSlipFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const vectorField nHat(this->patch().nf());
    const scalar t = this->db().time().timeOutputValue();

    Field<Type>::operator=
    (
        nHat*(nHat & uniformValue_->value(t))
      + transform(I - sqr(nHat), this->patchInternalField())
    );

    transformFvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::uniformFixedValueSlipFvPatchField<Type>::snGradTransformDiag() const
{
    const vectorField nHat(this->patch().nf());
    vectorField diag(nHat.size());

    diag.replace(vector::X, mag(nHat.component(vector::X)));
    diag.replace(vector::Y, mag(nHat.component(vector::Y)));
    diag.replace(vector::Z, mag(nHat.component(vector::Z)));

    return transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
}


template<class Type>
void Foam::uniformFixedValueSlipFvPatchField<Type>::write(Ostream& os) const
{
    transformFvPatchField<Type>::write(os);
    uniformValue_->writeData(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
