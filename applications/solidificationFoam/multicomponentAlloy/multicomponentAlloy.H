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

        //- Pointer to the dictionary of alloy components
        PtrDictionary<alloyComponent> components_;

        //- Solidus temperature [K]
        const dimensionedScalar solidus_;

    // Private Member functions

        //- Calculate a thermodynamic factor [J/m^3]
        //
        //  \tparam Phase  Index of phase: 0 - liquid, 1 - solid
        template<label Phase>
        dimensionedScalar factor() const;

        //- Calculate \sum A^2
        dimensionedScalar sumSqrA() const;

        //- Calculate \sum D_L
        dimensionedScalar sumD() const;

        //- Calculate \sum A^2/D_L
        dimensionedScalar sumSqrAperD() const;

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

            //- Return the molar mass of alloy [kg/mol]
            dimensionedScalar molarMass() const;

            //- Return the solidus temperature [K]
            dimensionedScalar solidus() const
            {
                return solidus_;
            }

            //- Return the liquidus temperature [K]
            dimensionedScalar liquidus() const
            {
                return liquidus_;
            }

            //- Return the solidification interval [K]
            dimensionedScalar deltaTemp() const
            {
                return sumSqrA()*factor<1>()/entropyChange_;
            }

            //- Return the interface energy [J/m^2]
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

            //- Return the mean capillary length [m]
            dimensionedScalar capillaryLength() const
            {
                return interfaceEnergy_/factor<0>()/sumD()/sumSqrAperD();
            }

            //- Return the Gibbs--Thomson coefficient [Km]
            dimensionedScalar GibbsThomson() const
            {
                return interfaceEnergy_/entropyChange_;
            }

            //- Return the mean diffusion coefficient [m^2/s]
            dimensionedScalar diffusionL() const
            {
                return sumSqrA()/sumSqrAperD();
            }

            //- Return the total sum of the growth restriction factors [Js/m^5]
            dimensionedScalar sumRestrictionFactors() const
            {
                return factor<0>()*sumSqrAperD();
            }

            //- Return the chemical driving force [J/m^3]
            tmp<volScalarField> chemicalDrivingForce
            (
                const volScalarField& phase,
                const volScalarField& T
            ) const;

            //- Return the partition function
            tmp<volScalarField> partition(const volScalarField& phase) const
            {
                return phase + (1 - phase)*factor<1>()/factor<0>();
            }

            //- Return the first derivative of the partition function
            dimensionedScalar partitionPrime() const
            {
                return 1 - factor<1>()/factor<0>();
            }

};

} // End namespace Foam

#ifdef NoRepository
#   include "multicomponentAlloyTemplates.C"
#endif

#endif
