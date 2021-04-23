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
#include "Cloud.H"
#include "constants.H"
#include "Random.H"

#include "rayTracingParticle.H"

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
    reflectionModelPtr_(reflectionModel::New(modelDict_.subDict("reflection")))
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::rayTracingHeatSource::calcSource()
{
    using constant::mathematical::pi;

    source_.primitiveFieldRef() = 0;
    Cloud<rayTracingParticle> cloud(mesh_, "cloud", IDLList<rayTracingParticle>());

    // Initialise the rayTracing particles
    const vector lPosition = position().value();
    const vector lDir = beam().direction();

    DebugInfo
        << "Laser position: " << lPosition << nl
        << "Laser direction: " << lDir << endl;

    vector rArea = Zero;
    {
        Random rnd(1234);
        scalar magr = 0.0;

        while (magr < VSMALL)
        {
            vector v = rnd.sample01<vector>();
            rArea = v - (v & lDir)*lDir;
            magr = mag(rArea);
        }
    }
    rArea.normalise();

    scalar dr = maxr_*radius().value()/nr_;
    scalar dTheta = 2*pi/nTheta_;
    scalar maxTrackLength = mesh_.bounds().mag();
    label nMissed = 0;
    point p0 = lPosition;
    scalar power = 0;
    scalar area = 0;
    point p1(p0);

    for (label ri = 0; ri < nr_; ri++)
    {
        // TODO: nonuniform step along radial axis (should be implemented in laserBeam)
        scalar r1 = SMALL + dr*ri;
        scalar r2 = r1 + dr;
        scalar rP = (r1 + r2)/2;
        vector localR = rP*rArea;

        scalar theta0 = 0;
        for (label thetai = 0; thetai < nTheta_; thetai++)
        {
            scalar theta1 = theta0 + SMALL + dTheta*thetai;
            scalar theta2 = theta1 + dTheta;
            scalar thetaP = (theta1 + theta2)/2;

            quaternion Q(lDir, thetaP);

            vector localPos = Q.R() & localR;
            p0 = lPosition + localPos;
            p1 = p0 + maxTrackLength*lDir/2;

            scalar Ip = beam().I(p0);
            scalar dAi = (sqr(r2) - sqr(r1))*dTheta/2;

            power += Ip*dAi;
            area += dAi;

            label cellI = mesh_.findCell(p0);

            if (cellI != -1)
            {
                // Create a new particle
                auto pPtr = autoPtr<rayTracingParticle>::New(mesh_, p0, p1, Ip, cellI, dAi, false);

                // Add to cloud
                cloud.addParticle(pPtr.release());
            }

            if (returnReduce(cellI, maxOp<label>()) == -1)
            {
                if (++nMissed <= 10)
                {
                    WarningInFunction
                        << "Cannot find owner cell for focalPoint at "
                        << p0 << endl;
                }
            }
        }
    }

    if (nMissed)
    {
        Info<< "Seeding missed " << nMissed << " locations" << endl;
    }

    DebugInfo
        << "Total power of the laser: " << power << nl
        << "Total area of the laser: " << area << nl
        << endl;


    tmp<volScalarField> treflectingCells = volScalarField::New
    (
        "reflectingCellsVol",
        mesh_,
        dimensionedScalar()
    );
    volScalarField& reflectingCellsVol = treflectingCells.ref();

    boolField reflectingCells(mesh_.nCells(), false);
    volScalarField A = mag(mixture_.gradAlphaM());
    const interpolationCell<scalar> AInterp(A);

    forAll(A, cellI)
    {
        if (A[cellI] > ROOTSMALL)
        {
            reflectingCells[cellI] = true;
            reflectingCellsVol[cellI] = 1;
        }
    }

    auto gradAlphaMInterpPtr =
        autoPtr<interpolationCellPoint<vector>>::New(mixture_.gradAlphaM());

    rayTracingParticle::trackingData td
    (
        cloud,
        AInterp,
        gradAlphaMInterpPtr,
        reflectingCells,
        *reflectionModelPtr_,
        source_
    );

    Info<< "Move particles..."
        << returnReduce(cloud.size(), sumOp<label>()) << endl;

    cloud.move(cloud, td, mesh_.time().deltaTValue());

    // Normalize by cell volume
    source_.primitiveFieldRef() /= mesh_.V();

    // TODO: remove this
    const volScalarField& redistribution = mixture_.surfaceHeatSourceRedistribution();
    source_ *= redistribution;

    if (debug)
    {
        Info<< "Final number of particles..."
            << returnReduce(cloud.size(), sumOp<label>()) << endl;

        OFstream osRef(type() + ":particlePath.obj");
        label vertI = 0;

        List<pointField> positions(Pstream::nProcs());
        List<pointField> p0(Pstream::nProcs());

        DynamicList<point>  positionsMyProc;
        DynamicList<point>  p0MyProc;

        for (const rayTracingParticle& p : cloud)
        {
            positionsMyProc.append(p.position());
            p0MyProc.append(p.p0());
        }

        positions[Pstream::myProcNo()].transfer(positionsMyProc);
        p0[Pstream::myProcNo()].transfer(p0MyProc);

        Pstream::gatherList(positions);
        Pstream::scatterList(positions);
        Pstream::gatherList(p0);
        Pstream::scatterList(p0);

        for (label proci = 0; proci < Pstream::nProcs(); ++proci)
        {
            const pointField& pos = positions[proci];
            const pointField& pfinal = p0[proci];
            forAll(pos, i)
            {
                meshTools::writeOBJ(osRef, pos[i]);
                vertI++;
                meshTools::writeOBJ(osRef, pfinal[i]);
                vertI++;
                osRef << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }

        osRef.flush();

        Info << "Total energy absorbed: " << fvc::domainIntegrate(source_).value() << endl;

        if (mesh_.time().outputTime())
        {
            reflectingCellsVol.write();
        }
    }
}


// ************************************************************************* //
