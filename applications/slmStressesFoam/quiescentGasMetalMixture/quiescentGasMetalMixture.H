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
    This file is part of slmStressesFoam.

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
    Foam::quiescentGasMetalMixture

Description
    Transport model for quiescent gas-metal mixture.

SourceFiles
    quiescentGasMetalMixture.C

\*---------------------------------------------------------------------------*/

#ifndef quiescentGasMetalMixture_H
#define quiescentGasMetalMixture_H

#include "IOdictionary.H"
#include "twoPhaseMixture.H"

#include "gasMetalThermalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class quiescentGasMetalMixture Declaration
\*---------------------------------------------------------------------------*/

class quiescentGasMetalMixture
:
    public IOdictionary,
    public twoPhaseMixture,
    public gasMetalThermalProperties
{
    // Private data

        // Phase1 density [kg/m^3]
        dimensionedScalar rho1_;

        // Phase2 density [kg/m^3]
        dimensionedScalar rho2_;

    // Private Member Functions

        //- No copy construct
        quiescentGasMetalMixture(const quiescentGasMetalMixture&) = delete;

        //- No copy assignment
        void operator=(const quiescentGasMetalMixture&) = delete;

public:

    // Static data

        // Declare name of the class and its debug switch
        ClassName("quiescentGasMetalMixture");

    //- Constructor
    quiescentGasMetalMixture(const fvMesh& mesh);

    //- Destructor
    virtual ~quiescentGasMetalMixture() = default;

    // Member functions

        //- Return phase1 density [kg/m^3]
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return phase2 density [kg/m^3]
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
