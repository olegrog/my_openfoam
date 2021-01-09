/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020-2021 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of gasMetalThermalProperties.

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
    gasMetalThermo

Description
    Thermophysical model for gas-metal mixture.
    This class is used only for calculating fields, but not for storing them.

SourceFiles
    gasMetalThermo.C
    gasMetalThermoTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef gasMetalThermo_H
#define gasMetalThermo_H

#include "IOdictionary.H"
#include "Polynomial.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class gasMetalThermo Declaration
\*---------------------------------------------------------------------------*/

class gasMetalThermo
:
    public IOdictionary
{

    struct kCp
    {
        //- Thermal conductivity polynomial coeffs
        const Polynomial<2> k;

        //- Specific heat at constant pressure polynomial coeffs
        const Polynomial<2> Cp;

        //- Constructor
        kCp(const dictionary& dict)
        :
            k(dict.lookup(coeffsName("k"))),
            Cp(dict.lookup(coeffsName("Cp")))
        {}
    };

    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Dictionaries
        const dictionary metalDict_, solidDict_, liquidDict_, gasDict_;

        //- Properties of all phases
        const kCp solid_, liquid_, gas_;

        //- Temperature of melting [K]
        const scalar Tmelting_;

        //- Temperature of boiling [K]
        const scalar Tboiling_;

        //- Latent heat of fusion [J/kg]
        const scalar Hfusion_;

        //- Latent heat of vapourisation [J/kg]
        const scalar Hvapour_;

        //- Molecular weight [g/mol]
        const scalar metalW_, gasW_;

    // Private Member Functions

        //- Coeffs name
        static word coeffsName(const char* name)
        {
            return word(name) + "Coeffs";
        }

        //- Return metal phase fraction (liquid/vapour) [-]
        template<class T1, class T2>
        tmp<volScalarField> metalPhaseFraction
        (
            const T1& h,
            const T1& hAtPhaseTransition,
            const T2& metalFraction,
            const scalar latentHeat,
            const word& phaseName
        ) const;

        //- Return derivative of a metal phase fraction w.r.t. enthalpy [kg/J]
        template<class T1, class T2>
        tmp<volScalarField> metalPhaseFractionPrimeEnthalpy
        (
            const T1& h,
            const T1& hAtPhaseTransition,
            const T2& metalFraction,
            const scalar latentHeat,
            const word& phaseName
        ) const;

        //- Return derivative of a metal phase fraction w.r.t. gas fraction [-]
        template<class T1, class T2>
        tmp<volScalarField> metalPhaseFractionPrimeGasFraction
        (
            const T1& h,
            const T1& hAtPhaseTransition,
            const T2& metalFraction,
            const scalar latentHeat,
            const scalar hPrePhaseTransitionPrime,
            const word& phaseName
        ) const;

        //- No copy construct
        gasMetalThermo(const gasMetalThermo&) = delete;

        //- No copy assignment
        void operator=(const gasMetalThermo&) = delete;


public:

    //- Constructor
    gasMetalThermo(const fvMesh& mesh);

    //- Destructor
    ~gasMetalThermo() = default;

    // Member Functions

        // Access

            //- Return latent heat of metal fusion [J/kg]
            dimensionedScalar Hfusion() const
            {
                return dimensionedScalar("Hfusion", dimEnergy/dimMass, Hfusion_);
            }

            //- Return latent heat of metal vapourisation [J/kg]
            dimensionedScalar Hvapour() const
            {
                return dimensionedScalar("Hvapour", dimEnergy/dimMass, Hvapour_);
            }

            //- Return temperature of melting [K]
            dimensionedScalar Tmelting() const
            {
                return dimensionedScalar("Tmelting", dimTemperature, Tmelting_);
            }

            //- Return temperature of boiling [K]
            dimensionedScalar Tboiling() const
            {
                return dimensionedScalar("Tboiling", dimTemperature, Tboiling_);
            }

            //- Return metal molar mass [kg/mol]
            dimensionedScalar metalM() const
            {
                return dimensionedScalar("metalM", dimMass/dimMoles, 1e-3*metalW_);
            }

            //- Return gas molar mass [kg/mol]
            dimensionedScalar gasM() const
            {
                return dimensionedScalar("gasM", dimMass/dimMoles, 1e-3*gasW_);
            }


        // Evaluation

            //- Heat capacity at constant pressure [J/kg/K]
            template<class T1, class T2, class T3>
            tmp<volScalarField> Cp
            (
                const T1& T,
                const T2& liquidFraction,
                const T3& gasFraction
            ) const;

            //- Thermal conductivity [W/m/K]
            template<class T1, class T2, class T3>
            tmp<volScalarField> k
            (
                const T1& T,
                const T2& liquidFraction,
                const T3& gasFraction
            ) const;

            //- Specific enthalpy [J/kg]
            template<class T1, class T2, class T3>
            tmp<volScalarField> h
            (
                const T1& T,
                const T2& liquidFraction,
                const T2& vapourFraction,
                const T3& gasFraction,
                const word& name = "h"
            ) const;

            //- Specific enthalpy at the melting temperature [J/kg]
            template<class T1>
            tmp<volScalarField> hAtMelting(const T1& gasFraction) const;

            //- Specific enthalpy at the boiling temperature [J/kg]
            template<class T1>
            tmp<volScalarField> hAtBoiling(const T1& gasFraction) const;

            //- Derivative of sensible enthalpy w.r.t. gas fraction [J/kg]
            template<class T1>
            tmp<volScalarField> sensibleEnthalpyPrimeGasFraction(const T1& T) const;

            //- Temperature [K]
            template<class T1, class T2, class T3>
            tmp<volScalarField> T
            (
                const T1& h,
                const T1& hAtMelting,
                const T2& liquidFraction,
                const T2& vapourFraction,
                const T3& gasFraction
            ) const;

            //- Liquid fraction [-]
            template<class T1, class T2>
            tmp<volScalarField> liquidFraction
            (
                const T1& h,
                const T1& hAtMelting,
                const T2& metalFraction
            ) const
            {
                return metalPhaseFraction(h, hAtMelting, metalFraction, Hfusion_, "liquid");
            }

            //- Vapour fraction [-]
            template<class T1, class T2>
            tmp<volScalarField> vapourFraction
            (
                const T1& h,
                const T1& hAtBoiling,
                const T2& metalFraction
            ) const
            {
                return metalPhaseFraction(h, hAtBoiling, metalFraction, Hvapour_, "vapour");
            }

            //- Derivative of liquid fraction w.r.t. enthalpy [kg/J]
            template<class T1, class T2>
            tmp<volScalarField> liquidFractionPrimeEnthalpy
            (
                const T1& h,
                const T1& hAtMelting,
                const T2& metalFraction
            ) const
            {
                return metalPhaseFractionPrimeEnthalpy
                (
                    h, hAtMelting, metalFraction, Hfusion_, "liquid"
                );
            }

            //- Derivative of vapour fraction w.r.t. enthalpy [kg/J]
            template<class T1, class T2>
            tmp<volScalarField> vapourFractionPrimeEnthalpy
            (
                const T1& h,
                const T1& hAtBoiling,
                const T2& metalFraction
            ) const
            {
                return metalPhaseFractionPrimeEnthalpy
                (
                    h, hAtBoiling, metalFraction, Hvapour_, "vapour"
                );
            }

            //- Derivative of liquid fraction w.r.t. gas fraction [-]
            template<class T1, class T2>
            tmp<volScalarField> liquidFractionPrimeGasFraction
            (
                const T1& h,
                const T1& hAtMelting,
                const T2& metalFraction
            ) const
            {
                const scalar hPreMeltingPrime =
                    gas_.Cp.integral(0, Tmelting_) - solid_.Cp.integral(0, Tmelting_);

                return metalPhaseFractionPrimeGasFraction
                (
                    h, hAtMelting, metalFraction, Hfusion_, hPreMeltingPrime, "liquid"
                );
            }

            //- Derivative of vapour fraction w.r.t. gas fraction [-]
            template<class T1, class T2>
            tmp<volScalarField> vapourFractionPrimeGasFraction
            (
                const T1& h,
                const T1& hAtBoiling,
                const T2& metalFraction
            ) const
            {
                const scalar hPreBoilingPrime =
                     gas_.Cp.integral(0, Tboiling_) - solid_.Cp.integral(0, Tmelting_)
                   - Hfusion_ - liquid_.Cp.integral(Tmelting_, Tboiling_);

                return metalPhaseFractionPrimeGasFraction
                (
                    h, hAtBoiling, metalFraction, Hvapour_, hPreBoilingPrime, "vapour"
                );
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "gasMetalThermoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
