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
    Foam::volumetricLaserHeatSource

Description
    Abstract base class for a volumetric heat source.

SourceFiles
    volumetricLaserHeatSource.C

\*---------------------------------------------------------------------------*/

#ifndef volumetricLaserHeatSource_H
#define volumetricLaserHeatSource_H

#include "runTimeSelectionTables.H"

#include "laserHeatSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class volumetricLaserHeatSource Declaration
\*---------------------------------------------------------------------------*/

class volumetricLaserHeatSource
:
    public laserHeatSource
{
protected:

    // Protected data

        //- Source dictionary
        const dictionary sourceDict_;

        //- Absorptivity coefficient
        const scalar A_;

public:

    //- Runtime type information
    TypeName("volumetricLaserHeatSource");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        volumetricLaserHeatSource,
        mesh,
        (const fvMesh& mesh),
        (mesh)
    );

    //- Constructor
    volumetricLaserHeatSource(const fvMesh& mesh);

    //- Selector
    static autoPtr<volumetricLaserHeatSource> New(const fvMesh& mesh);

    //- Destructor
    ~volumetricLaserHeatSource() override = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
