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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class multicomponentAlloy Declaration
\*---------------------------------------------------------------------------*/

class multicomponentAlloy
:
    public IOdictionary
{
    // Private data

        //- Molar mass of the base component [kg/mol]
        const dimensionedScalar baseMolarMass_;

        //- Latent heat of fusion [J/kg]
        const dimensionedScalar latentHeat_;

        //- Interface energy [J/m^2]
        const dimensionedScalar interfaceEnergy_;

        //- Liquidus temperature [K]
        const dimensionedScalar liquidus_;

        //- Density of the solid [kg/m^3]
        const dimensionedScalar rhoSolid_;

        //- Entropy change [J/m^3/K]
        const dimensionedScalar entropyChange_;

        //- List of phases
        const wordList phaseNames_;

        //- Pointer to the dictionary of alloy components
        PtrDictionary<alloyComponent> components_;

        //- Solidus temperature [K]
        const dimensionedScalar solidus_;

    // Private Member functions

        //- Calculate a thermodynamic factor [J/m^3]
        //
        //  \tparam Phase  Index of phase: 0 - liquid, 1 - solid
        template<class T1>
        auto factor(const word& phaseName, const T1& T) const -> decltype(T+T);

        //- Calculate \sum A^2
        template<class T1>
        auto sumSqrA(const T1& T) const -> decltype(T+T);

        //- Calculate \sum D_L
        dimensionedScalar sumDL() const;

        //- Calculate \sum A^2/D_L
        template<class T1>
        auto sumSqrAperDL(const T1& T) const -> decltype(T+T);

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

            //- Return the solidus temperature [K]
            const dimensionedScalar& solidus() const
            {
                return solidus_;
            }

            //- Return the liquidus temperature [K]
            const dimensionedScalar& liquidus() const
            {
                return liquidus_;
            }

            //- Return the interface energy [J/m^2]
            const dimensionedScalar& interfaceEnergy() const
            {
                return interfaceEnergy_;
            }

            //- Return the Gibbs--Thomson coefficient [Km]
            dimensionedScalar GibbsThomson() const
            {
                return interfaceEnergy_/entropyChange_;
            }

        // Evaluation

            //- Return the molar mass of alloy [kg/mol]
            dimensionedScalar molarMass() const;

            //- Return the solidification interval [K]
            dimensionedScalar deltaTemp() const
            {
                return sumSqrA(liquidus_)*factor("solid", liquidus_)/entropyChange_;
            }

            //- Return the dimensionless undercooling []
            template<class T1>
            auto undercooling(const T1& T) const -> decltype(T+T)
            {
                return (liquidus_ - T)/(liquidus_ - solidus_);
            }

            //- Return the mean capillary length [m]
            template<class T1>
            auto capillaryLength(const T1& T) const -> decltype(T+T)
            {
                return interfaceEnergy_/factor("liquid", T)/sumDL()/sumSqrAperDL(T);
            }

            //- Return the mean diffusion coefficient in liquid [m^2/s]
            template<class T1>
            auto diffusionL(const T1& T) const -> decltype(T+T)
            {
                return sumSqrA(T)/sumSqrAperDL(T);
            }

            //- Return the total sum of the growth restriction factors [Js/m^5]
            template<class T1>
            auto sumRestrictionFactors(const T1& T) const -> decltype(T+T)
            {
                return factor("liquid", T)*sumSqrAperDL(T);
            }

            //- Return the chemical driving force [J/m^3]
            tmp<volScalarField> chemicalDrivingForce
            (
                const volScalarField& phase,
                const volScalarField& T
            ) const;

            //- Return the partition function []
            tmp<volScalarField> partition
            (
                const volScalarField& phase,
                const volScalarField& T
            ) const
            {
                return phase + (1 - phase)*factor("solid", T)/factor("liquid", T);
            }

            //- Return the first derivative of the partition function []
            tmp<volScalarField> partitionPrime(const volScalarField& T) const
            {
                return 1 - factor("solid", T)/factor("liquid", T);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multicomponentAlloyTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
