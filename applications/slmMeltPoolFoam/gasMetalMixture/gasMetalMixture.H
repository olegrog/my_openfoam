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
    Foam::gasMetalMixture

Description
    Transport model for gas-metal mixture.

SourceFiles
    gasMetalMixture.C

\*---------------------------------------------------------------------------*/

#ifndef gasMetalMixture_H
#define gasMetalMixture_H

#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "Function1.H"

#include "gasMetalThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class gasMetalMixture Declaration
\*---------------------------------------------------------------------------*/

class gasMetalMixture
:
    public immiscibleIncompressibleTwoPhaseMixture
{
    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Thermophysical model
        const gasMetalThermo thermo_;

        //- Flag to write auxiliary fields
        const bool writeProperties_;

        // Fields

            //- Temperature [K]
            volScalarField T_;

            //- Specific enthalpy [J/kg]
            volScalarField h_;

            //- Specific enthalpy at the melting temperature [J/kg]
            volScalarField hAtMelting_;

            //- Liquid fraction with its derivatives
            LiquidFraction liquidFraction_;

            //- Derivative of specific enthalpy w.r.t. gas fraction [J/kg]
            volScalarField hPrimeGasFraction_;

            //- Heat capacity at constant pressure [J/kg/K]
            volScalarField Cp_;

            //- Thermal conductivity [W/m/K]
            volScalarField k_;

            //- Density [kg/m^3]
            volScalarField rho_;

        // Constants

            //- Surface tension as a function of temperature
            autoPtr<Function1<scalar>> pSigma_;

            //- Derivative of the surface tension w.r.t temperature [kg/s^2/K]
            const dimensionedScalar dSigmaDT_;

            //- Coefficient related to the permeability of mushy layer [kg/m^3/s]
            const dimensionedScalar mushyCoeff_;

    // Private Member Functions

        //- No copy construct
        gasMetalMixture(const gasMetalMixture&) = delete;

        //- No copy assignment
        void operator=(const gasMetalMixture&) = delete;

public:

    // Static data

        // Declare name of the class and its debug switch
        ClassName("gasMetalMixture");

    //- Constructor
    gasMetalMixture
    (
        const volVectorField& U,
        const surfaceScalarField& phi
    );

    //- Destructor
    virtual ~gasMetalMixture() = default;

    // Member Functions

        // Access

            //- Return the thermophysical model
            const gasMetalThermo& thermo() const
            {
                return thermo_;
            }

            //- Return the LiquidFraction object
            const LiquidFraction& liquidFraction() const
            {
                return liquidFraction_;
            }

            //- Return temperature [K]
            const volScalarField& T() const
            {
                return T_;
            }

            //- Return derivative of enthalpy w.r.t. gas fraction [J/kg]
            const volScalarField& hPrimeGasFraction() const
            {
                return hPrimeGasFraction_;
            }

            //- Return heat capacity at constant pressure [J/kg/K]
            const volScalarField& Cp() const
            {
                return Cp_;
            }

            //- Return thermal conductivity [W/m/K]
            const volScalarField& k() const
            {
                return k_;
            }

            //- Return enthalpy [J/kg]
            volScalarField& h()
            {
                return h_;
            }

            //- Return density [kg/m^3]
            volScalarField& rho()
            {
                return rho_;
            }

        // Evaluation

            //- Return the thermal Marangoni force [N/m^3]
            tmp<volVectorField> marangoniForce(const volVectorField& gradT) const;

            //- Return the solid-phase damping coefficient [N/m^3]
            tmp<volScalarField> solidPhaseDamping() const;

            //- Correct active internal fields
            void correct() override;

            //- Correct passive internal fields
            void finalCorrect();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
