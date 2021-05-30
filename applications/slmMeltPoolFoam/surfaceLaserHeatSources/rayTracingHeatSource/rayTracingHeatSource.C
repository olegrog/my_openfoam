/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                  laserDTRM | Copyright (C) 2017-2019 OpenCFD Ltd.
       rayTracingHeatSource | Copyright (C) 2021 Oleg Rogozin
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

\*---------------------------------------------------------------------------*/

#include "rayTracingHeatSource.H"

#include "addToRunTimeSelectionTable.H"
#include "constants.H"
#include "Random.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rayTracingHeatSource, 0);
    addToRunTimeSelectionTable(surfaceLaserHeatSource, rayTracingHeatSource, mixtureAdvector);

    defineTemplateTypeNameAndDebugWithName
    (
        Cloud<rayTracingParticle>,
        "rayTracingCloud",
        0
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rayTracingHeatSource::rayTracingHeatSource
(
    const incompressibleGasMetalMixture& mixture,
    const isoAdvection& advector
)
:
    surfaceLaserHeatSource(typeName, mixture, advector),
    nTheta_(modelDict_.get<label>("nTheta")),
    nr_(modelDict_.get<label>("nr")),
    maxr_(modelDict_.get<label>("maxr")),
    writeOBJ_(modelDict_.getOrDefault("writeOBJ", false)),
    scatteringModelPtr_(scatteringModel::New(modelDict_.subDict("scattering")))
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::rayTracingHeatSource::writeOBJ(const Cloud<rayTracingParticle>& cloud)
{
    const fileName outputFile
    (
        type() / word::printf("ParticlePath_%06d.obj", mesh_.time().timeIndex())
    );

    // Collect rays (straight lines) from all the processors
    List<DynamicList<Pair<point>>> allRays(Pstream::nProcs());
    for (const rayTracingParticle& p : cloud)
    {
        allRays[Pstream::myProcNo()].append(Pair<point>(p.p0(), p.position()));
    }
    Pstream::gatherList(allRays);

    if (Pstream::master())
    {
        mkDir(outputFile.path());
        OBJstream os(outputFile);

        forAll(allRays, proci)
        {
            for (const Pair<point>& p : allRays[proci])
            {
                os.write(linePointRef(p.first(), p.second()));
            }
        }

        os.flush();
    }
}


void Foam::rayTracingHeatSource::calcSource()
{
    using constant::mathematical::pi;

    startTimer();

    // 1. Prepare calculations

    const scalar dr = maxr_*radius().value()/nr_;
    const scalar dTheta = 2*pi/nTheta_;
    const scalar maxTrackLength = mesh_.bounds().mag();

    source_.primitiveFieldRef() = 0;
    label nMissed = 0;
    scalar totalPower = 0;
    scalar totalArea = 0;

    const scalar laserPower = power().value();
    const point laserPos = position().value();
    const vector laserDir = beam().direction();

    DebugInfo
        << "Laser position = " << laserPos << nl
        << "Laser direction = " << laserDir << endl;

    // A unit vector normal to the laser direction
    vector vHatInPlane = Zero;
    {
        Random rnd(1234);
        scalar magr = 0;

        while (magr < VSMALL)
        {
            const vector v = rnd.sample01<vector>();
            vHatInPlane = v - (v & laserDir)*laserDir;
            magr = mag(vHatInPlane);
        }
    }
    vHatInPlane.normalise();

    // 2. Generate a cloud of rayTracing particles

    Cloud<rayTracingParticle> cloud(mesh_, "cloud", IDLList<rayTracingParticle>());
    particle::particleCount_ = 0;

    for (label ri = 0; ri < nr_; ri++)
    {
        // TODO: nonuniform step along radial axis (should be implemented in laserBeam)
        const scalar r1 = SMALL + dr*ri;
        const scalar r2 = r1 + dr;
        const scalar rP = (r1 + r2)/2;
        const point localR = rP*vHatInPlane;

        for (label thetai = 0; thetai < nTheta_; thetai++)
        {
            const scalar theta1 = SMALL + dTheta*thetai;
            const scalar theta2 = theta1 + dTheta;
            const scalar thetaP = (theta1 + theta2)/2;

            const quaternion Q(laserDir, thetaP);

            const point localPos = Q.R() & localR;
            const point p0 = laserPos + localPos;
            const point p1 = p0 + maxTrackLength*laserDir;

            const scalar Ip = beam().I(p0)/laserPower; // [1/m2]
            const scalar dA = (sqr(r2) - sqr(r1))*dTheta/2; // [m2]
            const scalar dQ = Ip*dA; // fraction of the energy source transmitted by the ray

            totalPower += dQ*laserPower;
            totalArea += dA;

            const label cellI = mesh_.findCell(p0);

            if (cellI != -1)
            {
                // Create a new particle
                auto pPtr = autoPtr<rayTracingParticle>::New(mesh_, p0, p1, dQ, cellI, false);

                // Add to cloud
                cloud.addParticle(pPtr.release());
            }

            if (returnReduce(cellI, maxOp<label>()) == -1)
            {
                if (++nMissed <= 10)
                {
                    WarningInFunction
                        << "Cannot find owner cell for point at " << p0 << endl;
                }
            }
        }
    }

    if (nMissed)
    {
        Info<< "Seeding missed " << nMissed << " locations" << endl;
    }

    DebugInfo
        << "Total power of the laser = " << totalPower << nl
        << "Total area of the laser = " << totalArea << endl;

    // 3. Prepare the tracking data

    rayTracingParticle::trackingData td
    (
        cloud,
        const_cast<isoAdvection&>(advector_).surf(),
        scatteringModelPtr_,
        mixture_.alpha1(),
        mixture_.gradAlphaM(),
        source_
    );

    // 4. Evolve the cloud

    DebugInfo
        << "Initial number of particles = " << returnReduce(cloud.size(), sumOp<label>()) << endl;

    cloud.move(cloud, td, 0);

    DebugInfo
        << "Final number of particles = " << returnReduce(cloud.size(), sumOp<label>()) << endl;

    // 5. Dimensionalise the laser heat source

    source_.primitiveFieldRef() *= laserPower/mesh_.V();

    // 6. Dump particles trajectories

    if (writeOBJ_ && mesh_.time().writeTime())
    {
        writeOBJ(cloud);
    }

    DebugInfo << "Total power absorbed = " << fvc::domainIntegrate(source_).value() << endl;

    stopTimer();
}


// ************************************************************************* //
