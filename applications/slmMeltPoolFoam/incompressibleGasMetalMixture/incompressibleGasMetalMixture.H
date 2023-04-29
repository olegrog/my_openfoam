/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020-2023 Oleg Rogozin
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
    Foam::incompressibleGasMetalMixture

Description
    Transport model for incompressible immiscible gas-metal mixture.

SourceFiles
    incompressibleGasMetalMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleGasMetalMixture_H
#define incompressibleGasMetalMixture_H

#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "Function1.H"

#include "gasMetalThermalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class incompressibleGasMetalMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleGasMetalMixture
:
    public immiscibleIncompressibleTwoPhaseMixture,
    public gasMetalThermalProperties<incompressibleGasMetalMixture>
{
    // Private data

        // Constants

            //- Surface tension as a function of temperature [J/m^2]
            autoPtr<Function1<scalar>> sigmaPtr_;

            //- Coefficient related to the permeability of mushy layer [kg/m^3/s]
            const dimensionedScalar mushyCoeff_;

            //- Critical temperature for metal vapour [K]
            const dimensionedScalar Tcritical_;


    // Private Member Functions

        //- Calculates the Marangoni coefficient
        scalar dSigmaDT(scalar T) const;

        //- Returns the Marangoni coefficient [J/m^2/K]
        tmp<volScalarField> dSigmaDT() const;

        //- No copy construct
        incompressibleGasMetalMixture(const incompressibleGasMetalMixture&) = delete;

        //- No copy assignment
        void operator=(const incompressibleGasMetalMixture&) = delete;

public:

    //- Runtime type information
    TypeName("incompressibleGasMetalMixture");

    //- Constructor
    incompressibleGasMetalMixture
    (
        const volVectorField& U,
        const surfaceScalarField& phi
    );

    //- Destructor
    ~incompressibleGasMetalMixture() override = default;


    // Member Functions

        // Evaluation

            //- Return the thermal Marangoni force [N/m^3]
            tmp<volVectorField> marangoniForce() const;

            //- Return the solid-phase damping coefficient [N/m^3]
            tmp<volScalarField> solidPhaseDamping() const;

            //- Return the vapour pressure [Pa]
            tmp<volScalarField> vapourPressure(const dimensionedScalar& p0) const;

            //- Correct alpha-dependent fields
            void correct() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
