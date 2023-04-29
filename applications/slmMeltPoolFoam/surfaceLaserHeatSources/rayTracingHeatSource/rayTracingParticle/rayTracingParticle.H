/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
               DTRMParticle | Copyright (C) 2017-2019 OpenCFD Ltd
         rayTracingParticle | Copyright (C) 2021-2023 Oleg Rogozin
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
    Foam::rayTracingParticle

Description
    Particle for ray tracing

SourceFiles
    rayTracingParticleI.H
    rayTracingParticle.C
    rayTracingParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef rayTracingParticle_H
#define rayTracingParticle_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "volFieldsFwd.H"
#include "scatteringModel.H"
#include "reconstructionSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class rayTracingParticle;

// Forward declaration of friend functions
Ostream& operator<<(Ostream&, const rayTracingParticle&);

/*---------------------------------------------------------------------------*\
                        Class rayTracingParticle Declaration
\*---------------------------------------------------------------------------*/

class rayTracingParticle
:
    public particle
{
    // Private data

        //- Initial position [m]
        point p0_;

        //- Target position [m]
        point p1_;

        //- Initial fraction of energy source [-]
        scalar dQ0_;

        //- Remaining fraction of energy source [-]
        scalar dQ_;

        //- Is the particle being absorbed?
        bool isBeingAbsorbed_;

        //- Interaction type
        unsigned type_;


public:

    friend class Cloud<rayTracingParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Private data

            //- Flag to use the reconstructed interface
            const bool useSubCellData_;

            //- Reconstruction scheme of a geometric VoF
            const reconstructionSchemes& surf_;

            //- Scattering model
            const scatteringModel& scatteringModel_;

            //- Alpha metal field [-]
            const volScalarField& alphaM_;

            //- Gradient of the alpha metal field [1/m]
            const volVectorField& gradAlphaM_;

            //- Heat source term [W/m3]
            volScalarField& Q_;

            //- List of cells marked with some bits
            Map<bitSet> markedCells_;


    // Private Member Functions

        //- No copy construct
        trackingData(const trackingData&) = delete;

        //- No copy assignment
        void operator=(const trackingData&) = delete;

    public:

        enum cellProperty
        {
            ABSORBING, SCATTERING, NO_INTERFACE
        };

        //- Constructor
        inline trackingData
        (
            Cloud<rayTracingParticle>& cloud,
            bool useSubCellData,
            const reconstructionSchemes& surf,
            const scatteringModel& scattering,
            const volScalarField& alphaM,
            const volVectorField& gradAlphaM,
            volScalarField& Q
        );

        //- Destructor
        ~trackingData();


        // Member functions

            inline bool useSubCellData() const;
            inline const reconstructionSchemes& surf() const;
            inline const scatteringModel& scattering() const;
            inline scalar alphaM(label cellI) const;
            inline vector gradAlphaM(label cellI) const;
            inline scalar& Q(label cellI);
            inline void markCell(label cellI, cellProperty property);
    };


    // Static Data Members

        //- Runtime type information
        TypeName("rayTracingParticle");

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;

        //- String representation of properties
        AddToPropertyList
        (
            particle,
            " p0"
          + " p1"
          + " dQ0"
          + " dQ"
          + " isBeingAbsorbed";
        );


    // Constructors

        //- Construct from components, with searching for tetFace and
        //  tetPt unless disabled by doCellFacePt = false.
        rayTracingParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const vector& targetPosition,
            scalar dQ,
            label cellI,
            bool isBeingAbsorbed
        );

        //- Construct from Istream (for parallel transfer)
        rayTracingParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as copy
        rayTracingParticle(const rayTracingParticle& p);


        //- Factory class to read-construct particles (for parallel transfer)
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh) : mesh_(mesh) {}

            autoPtr<rayTracingParticle> operator()(Istream& is) const
            {
                return autoPtr<rayTracingParticle>
                (
                    new rayTracingParticle(mesh_, is, true)
                );
            }
        };


    // Access

        //- Return the initial position
        const point& p0() const
        {
            return p0_;
        }


    // Tracking

        //- Move
        bool move(Cloud<rayTracingParticle>&, trackingData&, scalar);


    // Member Functions

        //- Overridable function to handle the particle hitting a processorPatch
        void hitProcessorPatch(Cloud<rayTracingParticle>&, trackingData&);

        //- Overridable function to handle the particle hitting a wallPatch
        void hitWallPatch(Cloud<rayTracingParticle>&, trackingData&);

        //- Overridable function to handle the particle hitting a patch.
        //  Executed before other patch-hitting functions.
        bool hitPatch(Cloud<rayTracingParticle>&, trackingData&)
        {
            return false;
        }

        //- Write individual parcel properties to stream (for parallel transfer)
        void writeProperties
        (
            Ostream& os,
            const wordRes& filters,
            const word& delim,
            bool namesOnly = false
        ) const;


    // Ostream operator

        friend Ostream& operator<<(Ostream& os, const rayTracingParticle& p);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rayTracingParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
