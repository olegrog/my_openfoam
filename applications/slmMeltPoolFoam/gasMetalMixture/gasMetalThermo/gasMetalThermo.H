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
    This file is part of slmMeltPoolFoam.

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

#include "LiquidFraction.H"

namespace Foam
{

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

        //- No copy construct
        gasMetalThermo(const gasMetalThermo&) = delete;

        //- No copy assignment
        void operator=(const gasMetalThermo&) = delete;


public:

    //- Constructor
    gasMetalThermo
    (
        const word& metalPhase,
        const word& gasPhase,
        const fvMesh& mesh
    );

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
                const T3& gasFraction
            ) const;

            //- Derivative of enthalpy w.r.t. gas fraction [J/kg]
            template<class T1, class T2>
            tmp<volScalarField> hPrimeGasFraction
            (
                const T1& T,
                const T2& liquidFractionPrimeGasFraction
            ) const;

            //- Specific enthalpy at the melting temperature [J/kg]
            template<class T1>
            tmp<volScalarField> hAtMelting(const T1& gasFraction) const;

            //- Temperature [K]
            template<class T1, class T2, class T3>
            tmp<volScalarField> T
            (
                const T1& h,
                const T1& hAtMelting,
                const T2& liquidFraction,
                const T3& gasFraction
            ) const;

        //- Construct a LiquidFraction object
        LiquidFraction liquidFraction
        (
            const volScalarField& h,
            const volScalarField& hAtMelting,
            const volScalarField& metalFraction
        ) const;
};


} // End namespace Foam

#ifdef NoRepository
#   include "gasMetalThermoTemplates.C"
#endif

#endif
