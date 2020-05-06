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
    Class for calculation of liquid fraction in a gas-metal medium.

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
:
    public volScalarField
{
    // Private data

        //- Reference to the gas fraction field
        const volScalarField& metalFraction_;

        //- Reference to the enthalpy field
        const volScalarField& enthalpy_;

        //- Reference to the enthalpy of fusion
        const dimensionedScalar& enthalpyFusion_;

        //- Reference to the enthalpy at fusion point
        const volScalarField& enthalpyAtFusion_;

        //- Reference to the first derivative of enthalpy at fusion point
        const dimensionedScalar& enthalpyAtFusionPrime_;

        //- The sigmoid function
        autoPtr<sigmoidFunction> sigmoid_;

        //- The first derivative with respect to gas fraction
        volScalarField dAlpha_;

        //- The first derivative with respect to enthalpy
        volScalarField dEnthalpy_;

        //- The melt indicator (the maximum liquid fraction in time)
        volScalarField wasMelted_;

public:

    //- Constructor
    LiquidFraction
    (
        const volScalarField& metalFraction,
        const volScalarField& enthalpy,
        const volScalarField& enthalpyAtFusion,
        const dimensionedScalar& enthalpyFusion,
        const dimensionedScalar& enthalpyAtFusionPrime,
        const dictionary& thermalProperties
    );

    //- Destructor
    ~LiquidFraction() = default;

    // Member Functions

        //- Update the active internal fields
        void update();

        //- Update the passive internal fields
        void finalUpdate();

        //- Return the first derivative with respect to gas fraction
        tmp<volScalarField> dAlpha() const
        {
            return dAlpha_;
        }

        //- Return the first derivative with respect to enthalpy
        tmp<volScalarField> dEnthalpy() const
        {
            return dEnthalpy_;
        }
};

} // End namespace Foam

#endif
