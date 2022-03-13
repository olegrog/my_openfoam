/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                  laserDTRM | Copyright (C) 2017-2019 OpenCFD Ltd.
       rayTracingHeatSource | Copyright (C) 2021-2022 Oleg Rogozin
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
    Foam::rayTracingHeatSource

Description
    The surface laser heat source calculated from ray tracing.

SourceFiles
    rayTracingHeatSource.C

\*---------------------------------------------------------------------------*/

#ifndef rayTracingHeatSource_H
#define rayTracingHeatSource_H

#include "Cloud.H"
#include "Random.H"

#include "surfaceLaserHeatSource.H"
#include "scatteringModel.H"
#include "rayTracingParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class rayTracingHeatSource Declaration
\*---------------------------------------------------------------------------*/

class rayTracingHeatSource
:
    public surfaceLaserHeatSource
{
    // Private data

        //- Dictionary
        const dictionary dict_;

        //- Total number of rays
        const label nRays_;

        //- Ratio of the cut-off radius to the beam radius
        const scalar RcutByR_;

        //- Flag to write OBJ file for inspection in ParaView
        const bool writeOBJ_;

        //- Generator of random numbers
        Random random_;

        //- Flag to use the subcell data obtained from isoAdvector
        const bool useSubCellData_;

        //- The scattering model
        autoPtr<scatteringModel> scatteringModelPtr_;

    // Private Member Functions

        //- Write particle paths in the Wavefront OBJ file format
        void writeOBJ(const Cloud<rayTracingParticle>& cloud);

        //- Calculate the heat source
        void calcSource() final;

public:

    //- Runtime type information
    TypeName("rayTracing");

    //- Constructor
    rayTracingHeatSource
    (
        const incompressibleGasMetalMixture& mixture,
        const isoAdvection& advector
    );

    //- Destructor
    ~rayTracingHeatSource() final = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
