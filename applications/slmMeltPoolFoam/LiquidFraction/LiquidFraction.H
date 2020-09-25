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
    Foam::LiquidFraction

Description
    Class for calculation of liquid fraction in a gas-metal mixture.

SourceFiles
    LiquidFraction.C

\*---------------------------------------------------------------------------*/

#ifndef LiquidFraction_H
#define LiquidFraction_H

#include "volFields.H"

#include "sigmoidFunction.H"

namespace Foam
{

class LiquidFraction
{
    // Private data

        //- The sigmoid function
        autoPtr<sigmoidFunction> pSigmoid_;

        // References to fields

            //- Metal fraction
            const volScalarField& metalFraction_;

            //- Enthalpy
            const volScalarField& h_;

            //- Enthalpy at the melting temperature (depends on gas fraction)
            const volScalarField& hAtMelting_;

        // Dimensioned constants

            //- Latent heat of fusion
            const dimensionedScalar Hfusion_;

            //- Derivative of enthalpy w.r.t. gas fraction at the melting point
            const dimensionedScalar hAtMeltingPrime_;

        // Internal fields

            //- Liquid fraction in gas-metal mixture
            volScalarField liquidFraction_;

            //- Liquid fraction in metal
            volScalarField inMetal_;

            //- Derivative w.r.t. gas fraction
            volScalarField dGasFraction_;

            //- Derivative w.r.t. enthalpy
            volScalarField dEnthalpy_;

            //- Melt indicator (the maximum liquid fraction in time)
            volScalarField wasMelted_;

public:

    //- Constructor
    LiquidFraction
    (
        const volScalarField& metalFraction,
        const volScalarField& h,
        const volScalarField& hAtFusion,
        const dimensionedScalar& hFusion,
        const dimensionedScalar& hAtFusionPrime,
        const dictionary& dict
    );

    //- Destructor
    ~LiquidFraction() = default;

    // Member Functions

        //- Correct active internal fields
        void correct();

        //- Correct passive internal fields
        void finalCorrect();

        //- Return the liquid fraction in gas-metal mixture
        const volScalarField& operator()() const
        {
            return liquidFraction_;
        }

        //- Return the derivative w.r.t. gas fraction
        const volScalarField& dGasFraction() const
        {
            return dGasFraction_;
        }

        //- Return the derivative w.r.t. enthalpy
        const volScalarField& dEnthalpy() const
        {
            return dEnthalpy_;
        }

        //- Return the liquid fraction in the metal
        const volScalarField& inMetal() const
        {
            return inMetal_;
        }
};


} // End namespace Foam

#endif
