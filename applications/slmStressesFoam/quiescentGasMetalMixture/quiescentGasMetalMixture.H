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
    Transport model for a quiescent gas-metal mixture.

SourceFiles
    quiescentGasMetalMixture.C

\*---------------------------------------------------------------------------*/

#ifndef quiescentGasMetalMixture_H
#define quiescentGasMetalMixture_H

#include "quiescentTwoPhaseMixture.H"
#include "gasMetalThermalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class quiescentGasMetalMixture Declaration
\*---------------------------------------------------------------------------*/

class quiescentGasMetalMixture
:
    public quiescentTwoPhaseMixture,
    public gasMetalThermalProperties<quiescentGasMetalMixture>
{

    // Private Member Functions

        //- No copy construct
        quiescentGasMetalMixture(const quiescentGasMetalMixture&) = delete;

        //- No copy assignment
        void operator=(const quiescentGasMetalMixture&) = delete;

public:

    //- Runtime type information
    TypeName("quiescentGasMetalMixture");

    //- Constructor
    quiescentGasMetalMixture(const fvMesh& mesh);

    //- Destructor
    ~quiescentGasMetalMixture() override = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
