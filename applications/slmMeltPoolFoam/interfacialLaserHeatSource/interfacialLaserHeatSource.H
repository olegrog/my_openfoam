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
    Foam::interfacialLaserHeatSource

Description
    Laser heat source with selectable absorptivity model.

SourceFiles
    interfacialLaserHeatSource.C

\*---------------------------------------------------------------------------*/

#ifndef interfacialLaserHeatSource_H
#define interfacialLaserHeatSource_H

#include "IOdictionary.H"
#include "volFields.H"

#include "absorptivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class interfacialLaserHeatSource Declaration
\*---------------------------------------------------------------------------*/

class interfacialLaserHeatSource
:
    public IOdictionary
{
    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the time
        const Time& time_;

        //- Power distribution
        volScalarField powerDistribution_;

        //- Local absorptivity
        volScalarField absorptivity_;

        //- Reference to the metal fraction
        const volScalarField& alphaM_;

        //- Reference to the metal fraction gradient
        const volVectorField& gradAlphaM_;

        //- Laser radius
        const dimensionedScalar radius_;

        //- Laser power
        const dimensionedScalar power_;

        //- Laser scanning velocity
        const dimensionedVector velocity_;

        //- Laser position at zero time
        const dimensionedVector coordStart_;

        //- Time when the laser is turned off
        const dimensionedScalar timeStop_;

        //- Pointer to the absorptivity model
        autoPtr<absorptivityModel> absorptivityModelPtr_;

    // Private Member Functions

        //- Return a Gaussian profile
        tmp<volScalarField> surfaceGaussian
        (
            const volVectorField& x,
            dimensionedVector x0
        ) const;

        //- No copy construct
        interfacialLaserHeatSource(const interfacialLaserHeatSource&) = delete;

        //- No copy assignment
        void operator=(const interfacialLaserHeatSource&) = delete;

public:

    // Static data

        // Declare name of the class and its debug switch
        ClassName("interfacialLaserHeatSource");

    //- Constructor
    interfacialLaserHeatSource
    (
        const volScalarField& alphaM,
        const volVectorField& gradAlphaM
    );

    //- Destructor
    ~interfacialLaserHeatSource() = default;

    // Member Functions

        // Access

            //- Return the laser power
            const dimensionedScalar& power() const
            {
                return power_;
            }

            //- Return the laser radius
            const dimensionedScalar& radius() const
            {
                return radius_;
            }

            //- Return the laser coordinates
            dimensionedVector position() const
            {
                return coordStart_ + velocity_*time_;
            }

            //- Return the laser heat source
            tmp<volScalarField> value() const
            {
                return absorptivity_*powerDistribution_;
            }

        // Evaluation

            //- Update the time-dependent laser properties
            void update();

            //- Correct the alpha-dependent absorptivity
            void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
