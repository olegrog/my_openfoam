/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
               DTRMParticle | Copyright (C) 2017-2019 OpenCFD Ltd
         rayTracingParticle | Copyright (C) 2021 Oleg Rogozin
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
#include "interpolationCell.H"
#include "volFieldsFwd.H"
#include "reflectionModel.H"
#include "interpolationCellPoint.H"

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

        //- Initial position
        point p0_;

        //- Target position
        point p1_;

        //- Initial radiation intensity [1/m2]
        scalar I0_;

        //- Radiation intensity [1/m2]
        scalar I_;

        //- Area of radiation
        scalar dA_;

        //- Is particle transmissive?
        bool isTransmissive_;


public:

    friend class Cloud<rayTracingParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Interpolators for continuous phase fields

            const interpolationCell<scalar>& AInterp_;
            const interpolationCellPoint<vector>& gradAlphaMInterp_;

            //- Reflected cells
            const boolField& reflectingCells_;

            //- Reflection model
            const reflectionModel& reflection_;

            //- Heat source term
            volScalarField& Q_;


    public:

        // Constructors

            inline trackingData
            (
                Cloud<rayTracingParticle>& spc,
                const interpolationCell<scalar>& AInterp,
                const interpolationCellPoint<vector>& gradAlphaMInterp,
                const boolField&,
                const reflectionModel&,
                volScalarField& Q
            );


        // Member functions

            inline const interpolationCell<scalar>& AInterp() const;
            inline const interpolationCellPoint<vector>& gradAlphaMInterp() const;
            inline const boolField& reflectingCells() const;
            inline const reflectionModel& reflection() const;

            inline scalar& Q(label celli);
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
          + " I0"
          + " I"
          + " dA"
          + " isTransmissive";
        );


    // Constructors

        //- Construct from components, with searching for tetFace and
        //  tetPt unless disabled by doCellFacePt = false.
        rayTracingParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const vector& targetPosition,
            const scalar I,
            const label cellI,
            const scalar dA,
            const bool isTransmissive
        );

        //- Construct from components
        rayTracingParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const vector& position,
            const vector& targetPosition,
            const scalar I,
            const scalar dA,
            const bool isTransmissive
        );

        //- Construct from Istream
        rayTracingParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as copy
        rayTracingParticle(const rayTracingParticle& p);


        //- Factory class to read-construct particles used for parallel transfer
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

        //- Return const access to the initial position
        inline const point& p0() const;

        //- Return const access to the target position
        inline const point& p1() const;

        //- Return const access to the initial intensity
        inline scalar I0() const;

        //- Return const access to the current intensity
        inline scalar I() const;

        //- Return const access dA
        inline scalar dA() const;


    // Edit

        //- Return access to the target position
        inline point& p1();

        //- Return access to the initial intensity
        inline scalar& I0();

        //- Return access to the current intensity
        inline scalar& I();

        //- Return access to dA
        inline scalar& dA();

        //- Return access to reflectedId
        inline label& reflectedId();


    // Tracking

        //- Move
        bool move(Cloud<rayTracingParticle>& , trackingData&, const scalar);


    // Member Functions

        //- Overridable function to handle the particle hitting a processorPatch
        void hitProcessorPatch
        (
            Cloud<rayTracingParticle>&,
            trackingData& td
        );

        //- Overridable function to handle the particle hitting a wallPatch
        void hitWallPatch
        (
            Cloud<rayTracingParticle>&,
            trackingData& td
        );

        bool hitPatch
        (
            Cloud<rayTracingParticle>&,
            trackingData& td
        );


        // I-O

            //- Write individual parcel properties to stream
            void writeProperties
            (
                Ostream& os,
                const wordRes& filters,
                const word& delim,
                const bool namesOnly = false
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
