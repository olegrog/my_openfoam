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
    Foam::surfaceLaserHeatSource

Description
    Abstract base class for a surface heat source based on isoAdvection.

SourceFiles
    surfaceLaserHeatSource.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceLaserHeatSource_H
#define surfaceLaserHeatSource_H

#include "runTimeSelectionTables.H"
#include "isoAdvection.H"

#include "laserHeatSource.H"
#include "incompressibleGasMetalMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class surfaceLaserHeatSource Declaration
\*---------------------------------------------------------------------------*/

class surfaceLaserHeatSource
:
    public laserHeatSource
{

    // Private data

        //- Total time spent calculating the heat source
        scalar elapsedTime_;

        //- Time stamp for starting the heat source calculation
        scalar startTime_;


protected:

    // Protected data

        //- Heat source model dictionary
        const dictionary modelDict_;

        //- Reference to the mixture
        const incompressibleGasMetalMixture& mixture_;

        //- Reference to isoAdvector
        const isoAdvection& advector_;


    // Protected Member Functions

        //- Start timer (should be called in calcSource)
        void startTimer()
        {
            startTime_ = mesh_.time().elapsedCpuTime();
        }

        //- Stop timer (should be called in calcSource)
        void stopTimer()
        {
            elapsedTime_ += mesh_.time().elapsedCpuTime() - startTime_;
        }


public:

    //- Runtime type information
    TypeName("surfaceLaserHeatSource");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        surfaceLaserHeatSource,
        mixtureAdvector,
        (
            const incompressibleGasMetalMixture& mixture,
            const isoAdvection& advector
        ),
        (mixture, advector)
    );

    //- Constructor
    surfaceLaserHeatSource
    (
        const word& modelType,
        const incompressibleGasMetalMixture& mixture,
        const isoAdvection& advector
    );

    //- Selector
    static autoPtr<surfaceLaserHeatSource> New
    (
        const incompressibleGasMetalMixture& mixture,
        const isoAdvection& advector
    );

    //- Destructor
    ~surfaceLaserHeatSource() override = default;


    // Member functions

        //- Return total time spent calculating the heat source
        scalar elapsedTime() const
        {
            return elapsedTime_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
