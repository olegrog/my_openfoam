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
    This file is part of laserHeatSource.

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
    Foam::laserHeatSource

Description
    Abstract base class for a laser heat source.

SourceFiles
    laserHeatSource.C

\*---------------------------------------------------------------------------*/

#ifndef laserHeatSource_H
#define laserHeatSource_H

#include "volFields.H"

#include "laserProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class laserHeatSource Declaration
\*---------------------------------------------------------------------------*/

class laserHeatSource
:
    public laserProperties
{
    // Private Member Functions

        //- No copy construct
        laserHeatSource(const laserHeatSource&) = delete;

        //- No copy assignment
        void operator=(const laserHeatSource&) = delete;

        //- Calculate the heat source
        virtual void calcSource() = 0;

protected:

    // Protected data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Flag to write auxiliary fields
        const bool writeAllFields_;

        //- Heat source [W/m3]
        volScalarField source_;

public:

    //- Runtime type information
    TypeName("laserHeatSource");

    //- Constructor
    laserHeatSource(const fvMesh& mesh);

    //- Destructor
    ~laserHeatSource() = default;


    // Member Functions

        //- Return the laser heat source
        const volScalarField& operator()() const
        {
            return source_;
        }

        //- Correct the laser heat source
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
