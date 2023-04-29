/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                  laserDTRM | Copyright (C) 2017-2019 OpenCFD Ltd.
       rayTracingHeatSource | Copyright (C) 2021-2023 Oleg Rogozin
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
    nRays_(modelDict_.get<label>("nRays")),
    RcutByR_(modelDict_.get<scalar>("RcutByR")),
    writeOBJ_(modelDict_.getOrDefault("writeOBJ", false)),
    random_(modelDict_.getOrDefault("seed", 1234)),
    useSubCellData_(modelDict_.getOrDefault("useSubCellData", true)),
    scatteringModelPtr_(scatteringModel::New(modelDict_.subDict("scattering")))
{
    // 1. Check that threshold is low enough

    const scalar threshold = scatteringModelPtr_->threshold();
    const scalar densityRatio = (mixture.rho2()/mixture.rho1()).value();

    if (threshold > densityRatio)
    {
        FatalError
            << "Scattering threshold = " << threshold << " is not low enough." << nl
            << "It should be less than the density ratio = " << densityRatio
            << exit(FatalError);
    }

    // 2. Check that RcutByR is appropriate

    const label nRings = floor(sqrt(0. + nRays_));
    const scalar Rcut = RcutByR_*radius().value();
    const scalar P = power().value();
    const scalar dr = Rcut/nRings;
    const scalar lastIdr = (beam().integralI(Rcut) - beam().integralI(Rcut - dr))/P;
    const scalar cutFraction = 1 - beam().integralI(Rcut)/P;

    if (lastIdr*nRays_ < 3)
    {
        FatalError
            << "`RcutByR` is too large for the specified number of rays." << nl
            << "Fraction of the energy source at " << RcutByR_ << "*R = " << lastIdr
            << exit(FatalError);
    }

    if (cutFraction > 0.05)
    {
        Warning
            << "Fraction of the truncated energy source = " << cutFraction << endl;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::rayTracingHeatSource::writeOBJ(const Cloud<rayTracingParticle>& cloud)
{
    const fileName outputFile
    (
        type() / word::printf("ParticlePath_%08d.obj", mesh_.time().timeIndex())
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
            DebugInfo
                << " -- processor " << proci << " contains "
                << allRays[proci].size() << " particles" << endl;

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

    const label nRings = floor(sqrt(0. + nRays_));
    const scalar Rcut = RcutByR_*radius().value();
    const scalar dr = Rcut/nRings;
    const scalar maxTrackLength = mesh_.bounds().mag();

    source_.primitiveFieldRef() = 0;
    label nMissed = 0;
    scalar totalPower = 0;
    scalar totalArea = 0;

    const scalar laserPower = power().value();
    const point laserPos = position().value();
    const vector laserDir = beam().direction();
    const scalar cutFraction = 1 - beam().integralI(Rcut)/laserPower;

    DebugInfo
        << "Laser position = " << laserPos << nl
        << "Laser direction = " << laserDir << endl;

    // 2. Generate a cloud of rayTracing particles

    Cloud<rayTracingParticle> cloud(mesh_, "cloud", IDLList<rayTracingParticle>());
    particle::particleCount_ = 0;

    for (label ri = 0; ri < nRings; ri++)
    {
        // A unit vector normal to the laser direction
        vector vHatInPlane = Zero;
        {
            scalar magr = 0;

            while (magr < VSMALL)
            {
                const vector v = random_.sample01<vector>();
                vHatInPlane = v - (v & laserDir)*laserDir;
                magr = mag(vHatInPlane);
            }
        }
        vHatInPlane.normalise();

        const scalar r1 = dr*ri;
        const scalar r2 = r1 + dr;
        const scalar rP = (r1 + r2)/2;
        const point localR = rP*vHatInPlane;
        const scalar Idr = (beam().integralI(r2) - beam().integralI(r1))/laserPower; // []
        const label nTheta = ceil(nRays_*Idr);
        const scalar dTheta = 2*pi/nTheta;
        const scalar dA = (sqr(r2) - sqr(r1))*dTheta/2; // [m2]

        // Fraction of the energy source transmitted by the ray
        const scalar dQ = Idr/nTheta/(1 - cutFraction);

        for (label thetai = 0; thetai < nTheta; thetai++)
        {
            const scalar theta1 = dTheta*thetai;
            const scalar theta2 = theta1 + dTheta;
            const scalar thetaP = (theta1 + theta2)/2;

            const quaternion Q(laserDir, thetaP);

            const point localPos = Q.R() & localR;
            const point p0 = laserPos + localPos;
            const point p1 = p0 + maxTrackLength*laserDir;

            // const scalar Ip = beam().I(p0)/laserPower; // [1/m2]
            // dQ is approximately equal to Ip*dA

            const label cellI = mesh_.findCell(p0);

            if (cellI != -1)
            {
                if (mixture_.alpha1()[cellI] > SMALL)
                {
                    Warning
                        << "Particle is created in cell #" << cellI << " with alphaM = "
                        << mixture_.alpha1()[cellI] << endl;
                }

                totalPower += dQ*laserPower;
                totalArea += dA;

                // Create a new particle
                auto pPtr = autoPtr<rayTracingParticle>::New(mesh_, p0, p1, dQ, cellI, false);

                // Add to cloud
                cloud.addParticle(pPtr.release());
            }

            if (returnReduce(cellI, maxOp<label>()) == -1)
            {
                if (++nMissed <= 10)
                {
                    Warning
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
        << "Total emitted laser power = " << totalPower << " of " << laserPower << nl
        << "Total emission area = " << totalArea << " of " << pi*sqr(Rcut) << endl;

    // 3. Prepare the tracking data

    rayTracingParticle::trackingData td
    (
        cloud,
        useSubCellData_,
        const_cast<isoAdvection&>(advector_).surf(),
        scatteringModelPtr_,
        mixture_.alpha1(),
        mixture_.gradAlphaM(),
        source_
    );

    // 4. Evolve the cloud

    DebugInfo
        << "Initial number of particles = " << returnReduce(cloud.size(), sumOp<label>()) << endl;

    const bool oldThrowingErr = FatalError.throwing(true);
    try
    {
        cloud.move(cloud, td, 0);
    }
    catch (const error& err)
    {
        // Dump fields and particles in case of crash
        const_cast<Time&>(cloud.time()).writeNow();
        writeOBJ(cloud);
        returnReduceOr(true);   // MPI barrier
        const_cast<error&>(err).abort();
    }
    FatalError.throwing(oldThrowingErr);

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
