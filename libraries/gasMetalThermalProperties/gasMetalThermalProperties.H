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
    Foam::gasMetalThermalProperties

Description
    Base class for gas-metal mixture with thermal properties.

SourceFiles
    gasMetalThermalProperties.C
    gasMetalThermalPropertiesTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef gasMetalThermalProperties_H
#define gasMetalThermalProperties_H

#include "surfaceInterpolate.H"

#include "gasMetalThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class gasMetalThermalProperties Declaration
\*---------------------------------------------------------------------------*/

template<class Mixture>
class gasMetalThermalProperties
{
protected:

    // Protected data

        //- Thermophysical model
        const gasMetalThermo thermo_;

        //- Flag to write auxiliary fields
        const bool writeAllFields_;

        // References

            // Reference to the mixture
            const Mixture& mixture_;

            // Reference to the metal fraction
            const volScalarField& alphaM_;

            // Reference to the gas fraction
            const volScalarField& alphaG_;

        // Fields

            //- Gradient of the metal fraction [1/m]
            volVectorField gradAlphaM_;

            //- Temperature [K]
            volScalarField T_;

            //- Specific enthalpy [J/kg]
            volScalarField h_;

            //- Specific enthalpy at the melting temperature [J/kg]
            volScalarField hAtMelting_;

            //- Liquid fraction in the gas-metal mixture [-]
            volScalarField liquidFraction_;

            //- Derivative of the liquid fraction w.r.t. enthalpy [kg/J]
            volScalarField liquidFractionPrimeEnthalpy_;

            //- Heat capacity at constant pressure [J/kg/K]
            volScalarField Cp_;

            //- Thermal conductivity [W/m/K]
            volScalarField k_;

            //- Derivative of the sensible enthalpy w.r.t. gas fraction [W/m/K]
            volScalarField HsPrimeAlphaG_;

private:

    // Private data

        //- Flag to indicate whether to update redistribution_
        mutable bool updatedRedistr_;

        //- Heat source redistribution function
        mutable tmp<volScalarField> tredistribution_;


    // Private Member Functions

        //- Calculate metal fractions
        void calcMetalFractions();

        //- Calculate redistribution function for surface heat sources
        void calcRedistribution() const;

        //- No copy construct
        gasMetalThermalProperties(const gasMetalThermalProperties&) = delete;

        //- No copy assignment
        void operator=(const gasMetalThermalProperties&) = delete;

public:

    //- Runtime type information
    TypeName("gasMetalThermalProperties");

    //- Constructor
    gasMetalThermalProperties
    (
        const fvMesh& mesh,
        const Mixture& mixture
    );

    //- Destructor
    virtual ~gasMetalThermalProperties() = default;


    // Member Functions

        // Access

            //- Return the thermophysical model
            const gasMetalThermo& thermo() const
            {
                return thermo_;
            }

            //- Return gradient of the metal fraction [1/m]
            const volVectorField& gradAlphaM() const
            {
                return gradAlphaM_;
            }

            //- Return temperature [K]
            const volScalarField& T() const
            {
                return T_;
            }

            //- Return enthalpy [J/kg]
            volScalarField& h()
            {
                return h_;
            }

            //- Return liquid fraction [-]
            const volScalarField& liquidFraction() const
            {
                return liquidFraction_;
            }

            //- Return derivative of the liquid fraction w.r.t. enthalpy [kg/J]
            const volScalarField& liquidFractionPrimeEnthalpy() const
            {
                return liquidFractionPrimeEnthalpy_;
            }

            //- Return heat capacity at constant pressure [J/kg/K]
            const volScalarField& Cp() const
            {
                return Cp_;
            }

            //- Return face-interpolated heat capacity at constant pressure [J/kg/K]
            tmp<surfaceScalarField> Cpf() const
            {
                return fvc::interpolate(Cp_, "interpolate(Cp)");
            }

            //- Return thermal conductivity [W/m/K]
            const volScalarField& k() const
            {
                return k_;
            }

            //- Return face-interpolated thermal conductivity [W/m/K]
            tmp<surfaceScalarField> kf() const
            {
                return fvc::interpolate(k_, "interpolate(k)");
            }

            //- Return face-interpolated derivative of sensible enthalpy w.r.t. gas fraction [W/m/K]
            tmp<surfaceScalarField> HsPrimeAlphaGf() const
            {
                return fvc::interpolate(HsPrimeAlphaG_, "interpolate(HsPrimeAlphaG)");
            }

            //- Return redistribution function for surface heat sources
            const volScalarField& surfaceHeatSourceRedistribution() const;

        // Evaluation

            //- Correct alpha-dependent fields
            virtual void correct();

            //- Correct enthalpy-dependent fields
            virtual void correctThermo();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "gasMetalThermalPropertiesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
