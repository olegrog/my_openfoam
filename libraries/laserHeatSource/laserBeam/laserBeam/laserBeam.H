/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2021-2022 Oleg Rogozin
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
    Foam::laserBeam

Description
    Abstract base class for perfectly collimated laser beam properties.

SourceFiles
    laserBeam.C

\*---------------------------------------------------------------------------*/

#ifndef laserBeam_H
#define laserBeam_H

#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class laserProperties;

/*---------------------------------------------------------------------------*\
                       Class laserBeam Declaration
\*---------------------------------------------------------------------------*/

class laserBeam
{

protected:

    // Protected data

        //- Reference to laser properties
        const fvMesh& mesh_;

        //- Reference to laser properties
        const laserProperties& laser_;

        //- Beam direction
        const vector direction_;

public:

    //- Runtime type information
    TypeName("laserBeam");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        laserBeam,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const laserProperties& laser
        ),
        (dict, mesh, laser)
    );

    //- Constructor
    laserBeam
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const laserProperties& laser
    );

    //- Selector
    static autoPtr<laserBeam> New
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const laserProperties& laser
    );

    //- Destructor
    virtual ~laserBeam() = default;


    // Member Functions

        //- Return the direction of the laser beam
        const vector& direction() const
        {
            return direction_;
        }

        //- Intensity as a function of radial coordinate [W/m2]
        virtual scalar I(scalar r) const = 0;

        //- Integral of the intensity over the disk of radius r [W]
        virtual scalar integralI(scalar r) const = 0;

        //- Intensity field [W/m2]
        virtual tmp<volScalarField> I() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
