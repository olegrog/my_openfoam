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

namespace Foam
{

class alloyComponent
:
    public volScalarField
{
    // Private data

        word name_;
        dictionary alloyComponentDict_;
        dimensionedScalar rhoMelting_;
        dimensionedScalar molarMass_;
        dimensionedScalar equilibriumS_;
        dimensionedScalar equilibriumL_;
        dimensionedScalar slopeS_;
        dimensionedScalar slopeL_;
        dimensionedScalar diffusionS_;
        dimensionedScalar diffusionL_;
        dimensionedScalar Tmelting_;

public:
    // Constructors

        //- Construct from components
        alloyComponent
        (
            const word& name,
            const dictionary& alloyComponentDict,
            const fvMesh& mesh,
            const dimensionedScalar& Tmelting
        );

        //- Return a pointer to a new alloy component created on freestore from Istream
        class iNew
        {
            const fvMesh& mesh_;
            const dimensionedScalar& Tmelting_;

        public:

            iNew
            (
                const fvMesh& mesh,
                const dimensionedScalar& Tmelting
            )
            :
                mesh_(mesh),
                Tmelting_(Tmelting)
            {}

            autoPtr<alloyComponent> operator()(Istream& is) const
            {
                dictionaryEntry ent(dictionary::null, is);
                return autoPtr<alloyComponent>
                (
                    new alloyComponent(ent.keyword(), ent, mesh_, Tmelting_)
                );
            }
        };

        //- Return clone (need for PtrList)
        autoPtr<alloyComponent> clone() const
        {
            notImplemented("alloyComponent::clone() const");
            return autoPtr<alloyComponent>();
        }

    // Member Functions

        // Access

            const word& name() const
            {
                return name_;
            }

            const word& keyword() const
            {
                return name_;
            }

                //- Return limited concentration
            tmp<volScalarField> limitedConcentration() const
            {
                const volScalarField& concentration = *this;
                return Foam::min(Foam::max(concentration, scalar(0)), scalar(1));
            }

            //- Return the equilibrium concentration
            template<class T1, class T2>
            auto equilibrium(const T1& phase, const T2& T) const -> decltype(T+T)
            {
                return (equilibriumS_ + (T - Tmelting_)/slopeS_)*phase
                     + (equilibriumL_ + (T - Tmelting_)/slopeL_)*(1 - phase);
            }

            //- Return the diffusion coefficient
            tmp<volScalarField> diffusion(const volScalarField& phase) const
            {
                return diffusionS_ * phase + diffusionL_ * (1 - phase);
            }

            //- Return the slopes along the solidification path (solid and liquid)
            template<label Phase>
            dimensionedScalar slope() const
            {
                return Phase ? slopeL_ : slopeS_;
            }

            //- Return the diffusion coefficients (solid and liquid)
            template<label Phase>
            dimensionedScalar diffusion() const
            {
                return Phase ? diffusionL_ : diffusionS_;
            }

            //- Return the difference of the equilibrium concentrations
            dimensionedScalar deltaA() const
            {
                return equilibriumS_ - equilibriumL_;
            }

            //- Return the molar volume
            dimensionedScalar molarVolume() const
            {
                return molarMass_ / rhoMelting_;
            }

};

} // End namespace Foam

#endif
