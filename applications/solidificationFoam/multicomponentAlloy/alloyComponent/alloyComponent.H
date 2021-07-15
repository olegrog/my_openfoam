/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019-2021 Oleg Rogozin
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
    Abstract base class for concentration of a single alloy component.
    Used as part of the multicomponentAlloy for phase-field simulations.

SourceFiles
    alloyComponent.C

\*---------------------------------------------------------------------------*/

#ifndef alloyComponent_H
#define alloyComponent_H

#include "volFields.H"
#include "dictionaryEntry.H"
#include "PtrDictionary.H"

#include "componentPhase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class multicomponentAlloy;

/*---------------------------------------------------------------------------*\
                        Class alloyComponent Declaration
\*---------------------------------------------------------------------------*/

class alloyComponent
:
    public volScalarField
{
    // Private Data

        //- Field divided by partition coefficient
        autoPtr<volScalarField> fieldByKPtr_;

        //- Name of the component
        const word name_;

        //- Molar mass [kg/mol]
        const dimensionedScalar molarMass_;

        //- Reference to the alloy
        const multicomponentAlloy& alloy_;

        //- Pointer to the dictionary of alloy components
        PtrDictionary<componentPhase> phases_;

public:

    // Constructors

        //- Construct from components
        alloyComponent
        (
            const dictionaryEntry& entry,
            const fvMesh& mesh,
            const multicomponentAlloy& alloy
        );

        //- Helper for construction of alloyComponent PtrList/PtrDictionary
        class iNew
        {
            const fvMesh& mesh_;
            const multicomponentAlloy& alloy_;

        public:

            iNew
            (
                const fvMesh& mesh,
                const multicomponentAlloy& alloy
            )
            :
                mesh_(mesh),
                alloy_(alloy)
            {}

            autoPtr<alloyComponent> operator()(Istream& is) const
            {
                dictionaryEntry entry(dictionary::null, is);
                return autoPtr<alloyComponent>
                (
                    new alloyComponent(entry, mesh_, alloy_)
                );
            }
        };

        //- Return clone (need for PtrList)
        autoPtr<alloyComponent> clone() const
        {
            notImplemented("alloyComponent::clone() const");
            return autoPtr<alloyComponent>();
        }


    //- Selector
    static autoPtr<alloyComponent> New
    (
        const dictionaryEntry& entry,
        const fvMesh& mesh,
        const multicomponentAlloy& alloy
    );


    // Member Functions

        // Access

            //- Return name of the alloy component (need for PtrDictionary)
            const word& keyword() const
            {
                return name_;
            }

            //- Return the molar mass [kg/mol]
            const dimensionedScalar& molarMass() const
            {
                return molarMass_;
            }

            //- Return the phase diagram
            const componentPhase& phase(const word& phaseName) const
            {
                return phases_[phaseName];
            }

            //- Return reference to the alloy
            const multicomponentAlloy& alloy() const
            {
                return alloy_;
            }

            //- Return the field divided by the partition coefficient
            volScalarField& byK();

        // Evaluation

            //- Return the equilibrium concentration []
            template<class T1, class T2>
            auto equilibrium(const T1& solidFraction, const T2& T) const -> decltype(T+T)
            {
                return solidFraction*phases_["solid"].equilibrium(T)
                    + (1 - solidFraction)*phases_["liquid"].equilibrium(T);
            }

            //- Return the difference of the equilibrium concentrations []
            template<class T1>
            auto delta(const T1& T) const -> decltype(T+T)
            {
                return phases_["solid"].equilibrium(T) - phases_["liquid"].equilibrium(T);
            }

            //- Return the derivative of difference of the equilibrium concentrations [1/K]
            template<class T1>
            auto deltaPrime(const T1& T) const -> decltype(T+T)
            {
                return phases_["solid"].equilibriumPrime(T)
                    - phases_["liquid"].equilibriumPrime(T);
            }

            //- Return the derivative of the equilibrium concentrations w.r.t temperature [1/K]
            template<class T1, class T2>
            auto equilibriumPrimeT(const T1& solidFraction, const T2& T) const -> decltype(T+T)
            {
                return solidFraction/phases_["solid"].slope(T)
                    + (1 - solidFraction)/phases_["liquid"].slope(T);
            }

            //- Return the diffusion coefficient [m2/s]
            tmp<volScalarField> diffusion
            (
                const volScalarField& solidFraction,
                const volScalarField& partition
            ) const
            {
                return phases_["solid"].diffusion()*solidFraction
                    + phases_["liquid"].diffusion()*(1 - solidFraction)*partition;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
