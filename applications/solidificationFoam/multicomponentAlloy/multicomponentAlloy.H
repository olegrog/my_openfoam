/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019-2020 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of solidificationFoam.

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
    Foam::multicomponentAlloy

Description
    Multicomponent alloy.

SourceFiles
    multicomponentAlloy.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentAlloy_H
#define multicomponentAlloy_H

#include "volFields.H"
#include "PtrDictionary.H"

#include "alloyComponent.H"

namespace Foam
{

class multicomponentAlloy
:
    public IOdictionary
{
    // Private data

        const dimensionedScalar entropyChange_, interfaceEnergy_;
        const dimensionedScalar solidus_, liquidus_;
        PtrDictionary<alloyComponent> components_;
        const dimensionedScalar factorS_;
        const dimensionedScalar factorL_;

    // Private Member functions

        template<label Phase>
        dimensionedScalar calcFactor() const
        {
            auto iter = components_.begin();

            dimensionedScalar result = iter().deltaA()/iter().slope<Phase>()*iter().molarVolume();

            for (++iter; iter != components_.end(); ++iter)
            {
                result += iter().deltaA()/iter().slope<Phase>()*iter().molarVolume();
            }

            return dimensionedScalar(entropyChange_/result);
        }

public:

    //- Constructor
    multicomponentAlloy(const fvMesh& mesh);

    //- Destructor
    ~multicomponentAlloy() = default;

    // Member Functions

        //- Return the alloy components
        PtrDictionary<alloyComponent>& components()
        {
            return components_;
        }

        // Access

            //- Return the solidus temperature
            dimensionedScalar solidus() const
            {
                return solidus_;
            }

            //- Return the liquidus temperature
            dimensionedScalar liquidus() const
            {
                return liquidus_;
            }

            //- Return the interface energy
            dimensionedScalar interfaceEnergy() const
            {
                return interfaceEnergy_;
            }

            //- Return the dimensionless undercooling
            template<class T1>
            auto undercooling(const T1& T) const -> decltype(T+T)
            {
                return (liquidus_ - T)/(liquidus_ - solidus_);
            }

            //- Return the relaxation time
            dimensionedScalar relaxationTime() const;

            //- Return the mean capillary length
            dimensionedScalar capillaryLength() const;

            //- Return the mean diffusion coefficient
            dimensionedScalar diffusionL() const;

            //- Return the chemical driving force of the solidification
            tmp<volScalarField> chemicalDrivingForce
            (
                const volScalarField& phase,
                const volScalarField& T
            ) const;

            //- Return the partition function
            tmp<volScalarField> partition(const volScalarField& phi) const
            {
                return phi + (1 - phi)*factorS_/factorL_;
            }

            dimensionedScalar partitionPrime() const
            {
                return 1 - factorS_/factorL_;
            }

};

} // End namespace Foam

#endif
