/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of stereolithography.

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
    Foam::laserScanner

Description
    Properties of laser 3D scanning system.

SourceFiles
    laserScanner.C

\*---------------------------------------------------------------------------*/

#ifndef laserScanner_H
#define laserScanner_H

#include "IOdictionary.H"
#include "dimensionedScalar.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class laserScanner Declaration
\*---------------------------------------------------------------------------*/

class laserScanner
:
    public IOdictionary
{
    // Private data

        //- Radius of the beam [m]
        const dimensionedScalar radius_;

        //- Power [W]
        const dimensionedScalar power_;

        //- Speed of surface hatching [m^2/s]
        const dimensionedScalar hatchingSpeed_;

        //- Thickness of layers [m]
        const dimensionedScalar layerThickness_;

        //- Thickness of the first layer [m]
        const dimensionedScalar firstLayerThickness_;

public:

    //- Constructor
    laserScanner(const fvMesh& mesh);

    //- Destructor
    ~laserScanner() = default;

    // Member Functions

        // Access

            //- Return the radius of the beam [m]
            dimensionedScalar R() const
            {
                return radius_;
            }

            //- Return the exposure on the surface of each layer [J/m^2]
            dimensionedScalar E() const
            {
                return power_/hatchingSpeed_;
            }

        // Evaluation

            //- Return the hatching speed
            dimensionedScalar V(const dimensionedScalar& L) const
            {
                return hatchingSpeed_/L;
            }

            //- Return the coordinate of i-th layer
            dimensionedScalar height(label i) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
