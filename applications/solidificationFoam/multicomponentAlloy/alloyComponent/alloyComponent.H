/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::alloyComponent

Description
    Concentration of the single alloy component.
    Used as part of the multicomponentAlloy for phase-field simulations.

SourceFiles
    alloyComponent.C

\*---------------------------------------------------------------------------*/

#ifndef alloyComponent_H
#define alloyComponent_H

#include "volFields.H"
#include "dictionaryEntry.H"

namespace Foam {

class alloyComponent : public volScalarField {
public:
    // Constructors

    //- Construct from components
    alloyComponent(
        const word& name,
        const dictionary& alloyComponentDict,
        const fvMesh& mesh,
        const dimensionedScalar& meltingTemp
    );

    //- Return a pointer to a new alloy component created on freestore from Istream
    class iNew {
    public:
        iNew(
            const fvMesh& mesh,
            const dimensionedScalar& meltingTemp
        ) :
            mesh_(mesh),
            meltingTemp_(meltingTemp)
        {}
        autoPtr<alloyComponent> operator()(Istream& is) const {
            dictionaryEntry ent(dictionary::null, is);
            return autoPtr<alloyComponent>(new alloyComponent(
                ent.keyword(), ent, mesh_, meltingTemp_));
        }
    private:
        const fvMesh& mesh_;
        const dimensionedScalar& meltingTemp_;
    };

    //- Return clone
    autoPtr<alloyComponent> clone() const;

    // Member Functions

    const word& name() const { return name_; }

    const word& keyword() const { return name_; }

    //- Return limited concentration
    tmp<volScalarField> limitedConcentration() const {
        const volScalarField& concentration = *this;
        return Foam::min(Foam::max(concentration, scalar(0)), scalar(1));
    }

    //- Return the equilibrium concentration
    tmp<volScalarField> equilibrium(
        const volScalarField& phase,
        const volScalarField& T
    ) const {
        return (equilibriumS_ + (T - meltingTemp_) / slopeS_) * phase
             + (equilibriumL_ + (T - meltingTemp_) / slopeL_) * (1 - phase);
    }

    //- Return the diffusion coefficient
    tmp<volScalarField> diffusion(const volScalarField& phase) const {
        return diffusionS_ * phase + diffusionL_ * (1 - phase);
    }

    //- Return the slopes along the solidification path (solid and liquid)
    template<label Phase> dimensionedScalar slope() const;
    template<> dimensionedScalar slope<0>() const { return slopeS_; }
    template<> dimensionedScalar slope<1>() const { return slopeL_; }

    //- Return the diffusion coefficients (solid and liquid)
    template<label Phase> dimensionedScalar diffusion() const;
    template<> dimensionedScalar diffusion<0>() const { return diffusionS_; }
    template<> dimensionedScalar diffusion<1>() const { return diffusionL_; }

    //- Return the difference of the equilibrium concentrations
    dimensionedScalar deltaA() const {
        return equilibriumS_ - equilibriumL_;
    }

    dimensionedScalar molarVolume() const {
        return atomicMass_ / densityMelting_;
    }

private:
    word name_;
    dictionary alloyComponentDict_;
    dimensionedScalar densityMelting_;
    dimensionedScalar atomicMass_;
    dimensionedScalar equilibriumS_;
    dimensionedScalar equilibriumL_;
    dimensionedScalar slopeS_;
    dimensionedScalar slopeL_;
    dimensionedScalar diffusionS_;
    dimensionedScalar diffusionL_;
    dimensionedScalar meltingTemp_;
};

} // End namespace Foam

#endif
