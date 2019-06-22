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
    Foam::multicomponentAlloy

Description
    Multicomponent alloy.

SourceFiles
    multicomponentAlloy.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentAlloy_H
#define multicomponentAlloy_H

#include "alloyComponent/alloyComponent.H"
#include "PtrDictionary.H"

namespace Foam {

class multicomponentAlloy : public IOdictionary {
public:

    // Constructors
    multicomponentAlloy(const fvMesh& mesh);

    // Member Functions

    //- Return the alloy components
    PtrDictionary<alloyComponent>& components() {
        return components_;
    }

    //- Return the solidus temperature
    dimensionedScalar solidus() const {
        return solidus_;
    }

    //- Return the liquidus temperature
    dimensionedScalar liquidus() const {
        return liquidus_;
    }

    //- Return the interface energy
    dimensionedScalar interfaceEnergy() const {
        return interfaceEnergy_;
    }

    //- Return the dimensionless undercooling
    dimensionedScalar undercooling(const dimensionedScalar& T) const {
        return (liquidus_ - T) / (liquidus_ - solidus_);
    }

    //- Return the dimensionless undercooling
    tmp<volScalarField> undercooling(const volScalarField& T) const {
        return (liquidus_ - T) / (liquidus_ - solidus_);
    }

    //- Return the relaxation time
    dimensionedScalar relaxationTime() const;

    //- Return the mean capillary length
    dimensionedScalar capillaryLength() const;

    //- Return the mean diffusion coefficient
    dimensionedScalar diffusionL() const;

    //- Return the chemical driving force of the solidification
    tmp<volScalarField> chemicalDrivingForce(
        const volScalarField& phase,
        const volScalarField& T
    ) const;

    //- Return the partition function
    tmp<volScalarField> partition(const volScalarField& phi) const {
        return phi + (1 - phi) * factorS_ / factorL_;
    }

    dimensionedScalar partitionPrime() const {
        return 1 - factorS_ / factorL_;
    }

private:
    template<label Phase>
    dimensionedScalar calcFactor() const {
        PtrDictionary<alloyComponent>::const_iterator iter = components_.begin();
        dimensionedScalar result = iter().deltaA() / iter().slope<Phase>();
        for (++iter; iter != components_.end(); ++iter) {
            result += iter().deltaA() / iter().slope<Phase>();
        }
        return dimensionedScalar(entropyChange_ / result);
    }

private:
    const dimensionedScalar entropyChange_, interfaceEnergy_;
    const dimensionedScalar solidus_, liquidus_;
    PtrDictionary<alloyComponent> components_;
    const dimensionedScalar factorS_;
    const dimensionedScalar factorL_;
};

} // End namespace Foam

#endif
