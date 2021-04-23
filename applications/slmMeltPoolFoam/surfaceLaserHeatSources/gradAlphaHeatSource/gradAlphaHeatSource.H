/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2021 Oleg Rogozin
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
    Foam::gradAlphaHeatSource

Description
    The surface laser heat source calculated from the metal fraction gradient.

SourceFiles
    gradAlphaHeatSource.C

\*---------------------------------------------------------------------------*/

#ifndef gradAlphaHeatSource_H
#define gradAlphaHeatSource_H

#include "surfaceLaserHeatSource.H"
#include "absorptionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class gradAlphaHeatSource Declaration
\*---------------------------------------------------------------------------*/

class gradAlphaHeatSource
:
    public surfaceLaserHeatSource
{
    // Private data

        //- Dictionary
        const dictionary dict_;

        //- The absorption model
        autoPtr<absorptionModel> absorptionModelPtr_;


    // Private Member Functions

        //- Calculate the heat source
        void calcSource() final;

public:

    //- Runtime type information
    TypeName("gradAlpha");

    //- Constructor
    gradAlphaHeatSource
    (
        const incompressibleGasMetalMixture& mixture,
        const isoAdvection& advector
    );

    //- Destructor
    ~gradAlphaHeatSource() final = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
