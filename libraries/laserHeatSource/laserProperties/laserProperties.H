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
    Foam::laserProperties

Description
    Properties of the laser source and scanning system.

SourceFiles
    laserProperties.C

\*---------------------------------------------------------------------------*/

#ifndef laserProperties_H
#define laserProperties_H

#include "IOdictionary.H"
#include "dimensionedTypes.H"
#include "Time.H"
#include "Function1.H"

#include "laserBeam.H"
#include "movingReferenceFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class laserProperties Declaration
\*---------------------------------------------------------------------------*/

class laserProperties
:
    public IOdictionary
{

private:

    // Private data

        //- Reference to the time
        const Time& time_;

        //- Laser power as a function of time
        autoPtr<Function1<scalar>> powerPtr_;

        //- Laser radius as a function of time
        autoPtr<Function1<scalar>> radiusPtr_;

        //- Laser scanning velocity as a function of time
        autoPtr<Function1<vector>> velocityPtr_;

        //- Laser position at zero time
        const dimensionedVector coordStart_;

        //- Time when the laser is turned off
        const dimensionedScalar timeStop_;

        //- Direction of the laser beam
        autoPtr<laserBeam> beamPtr_;

        //- Dictionary for serialization
        IOdictionary laserDict_;

        //- Current laser position
        dimensionedVector position_;

        //- Pointer to moving reference frame
        const movingReferenceFrame* MRF_;


    // Private Member Functions

        //- No copy construct
        laserProperties(const laserProperties&) = delete;

        //- No copy assignment
        void operator=(const laserProperties&) = delete;

public:

    //- Runtime type information
    TypeName("laserProperties");

    //- Constructor
    laserProperties(const fvMesh& mesh);

    //- Destructor
    ~laserProperties() = default;


    // Member Functions

        //- Access

            //- Return the laser power
            dimensionedScalar power() const
            {
                return dimensionedScalar("power", dimPower, powerPtr_->value(time_.value()));
            }

            //- Return the laser radius
            dimensionedScalar radius() const
            {
                return dimensionedScalar("radius", dimLength, radiusPtr_->value(time_.value()));
            }

            //- Return the laser velocity
            dimensionedVector velocity() const;

            //- Return the laser coordinates
            dimensionedVector position() const
            {
                return position_;
            }

            //- Return true when laser source is swithed on
            bool switchedOn() const
            {
                return time_ < timeStop_;
            }

            //- Return reference to laser beam
            const laserBeam& beam() const
            {
                return *beamPtr_;
            }


        //- Evaluation

            //- Update the laser coordinates
            void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
