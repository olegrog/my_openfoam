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

\*---------------------------------------------------------------------------*/

#include "rayTracingParticle.H"
#include "constants.H"
#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rayTracingParticle, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rayTracingParticle::rayTracingParticle
(
    const polyMesh& mesh,
    const vector& position,
    const vector& targetPosition,
    const scalar dQ,
    const label cellI,
    const bool isBeingAbsorbed
)
:
    particle(mesh, position, cellI),
    p0_(position),
    p1_(targetPosition),
    dQ0_(dQ),
    dQ_(dQ),
    isBeingAbsorbed_(isBeingAbsorbed),
    type_(0)
{}


Foam::rayTracingParticle::rayTracingParticle(const rayTracingParticle& p)
:
    particle(p),
    p0_(p.p0_),
    p1_(p.p1_),
    dQ0_(p.dQ0_),
    dQ_(p.dQ_),
    isBeingAbsorbed_(p.isBeingAbsorbed_),
    type_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::rayTracingParticle::move
(
    Cloud<rayTracingParticle>& cloud,
    trackingData& td,
    scalar
)
{
    // Reset particle properties
    td.switchProcessor = false;
    td.keepParticle = true;

    // Do nothing if particle is terminated (parallel transfer case)
    if (1 - stepFraction() < SMALL)
    {
        return td.keepParticle;
    }

    // Parameter surfCellTol is used in plicRDF.C and gradAlpha.C
    const scalar alphaTol = td.surf().modelDict().getOrDefault<scalar>("surfCellTol", 1e-8);
    const scalar maxTrackLength = mesh().bounds().mag();
    const scalar smallLength = SMALL*maxTrackLength;
    const vector s = p1_ - p0_;
    const label maxCloudSize = 1e5;
    const label maxIter = 1e5;
    label iter = 0;

    const word procId = Pstream::parRun() ? word::printf("<%d>", origProc()) : "";
    DebugPout
        << "Particle #" << origId() << procId << " is moving along " << s/mag(s)
        << " and transmitting " << dQ_ << " of energy" << endl;

    if (stepFraction() > SMALL)
    {
        DebugPout
            << " --- just comes from another processor with stepFraction = "
            << stepFraction() << endl;
    }

    if (dQ0_ - dQ_ > SMALL)
    {
        DebugPout << " --- dQ0_ = " << dQ0_ << endl;
    }

    // --- Helper functions

    // Function for translating energy from the particle to the heat source in the cell
    auto absorb = [&td, this](scalar cellI, scalar ds, scalar absorptionPath)
    {
        const scalar cellSize = cbrt(mesh().cellVolumes()[cellI]);
        scalar absorptionFraction = absorptionPath*ds/cellSize;

        // Do not absorb more than alpha to prevent local overheating
        absorptionFraction = min(absorptionFraction, td.alphaM(cellI));

        // Do not absorb more energy than the particle has
        scalar deltaQ = min(dQ0_*absorptionFraction, dQ_);

        if (deltaQ < SMALL) return;

        td.Q(cellI) += deltaQ;
        dQ_ -= deltaQ;

        td.markCell(cellI, trackingData::ABSORBING);

        DebugPout
            << " --- absorbed energy = " << deltaQ
            << ", absorptionFraction = " << absorptionFraction
            << ", alphaM = " << td.alphaM(cellI) << endl;
    };

    // Function for stopping particle advancement
    auto terminate = [&td, this]()
    {
        stepFraction() = 1;
        DebugPout << " +++ particle is terminated" << endl;
    };

    // Function for adding a new particle to the cloud
    auto addParticle = [&td, &cloud, this, maxTrackLength]
    (
        scalar cellI,
        scalar energyFraction,
        const point& startingPoint,
        const vector& incidentDir,
        const vector& nHat,
        bool isRefracted
    )
    {
        const word type = isRefracted ? "refract" : "reflect";
        const vector direction =
            isRefracted
          ? td.scattering().refraction(incidentDir, nHat)
          : td.scattering().reflection(incidentDir, nHat);
        const point finalPoint = startingPoint + direction*maxTrackLength;

        auto pPtr = autoPtr<rayTracingParticle>::New
        (
            mesh(), startingPoint, finalPoint, energyFraction*dQ_, cellI, isRefracted
        );
        pPtr->reset(); // to set stepFraction = 0

        DebugPout
            << " --- particle #" << pPtr->origId() << " is " << type
            << "ed into direction " << direction << endl;

        cloud.addParticle(pPtr.release());

        // In some rare cases (e.g., grazing incidence), the reflection of the rays is looped
        if (cloud.size() == maxCloudSize - 100)
        {
            SeriousError << "Too many particles have been spawn!" << endl;
            debug = 1;
        }
        if (cloud.size() > maxCloudSize)
        {
            FatalError << "Too many particles have been spawn!" << exit(FatalError);
        }
    };

    // Function returning nHat directed toward the metal using VoF subcell information
    auto getVoFNormal = [&td](scalar cellI)
    {
        vector nHat = td.surf().normal()[cellI];
        nHat.normalise();
        DebugPout << " --- nHat (VoF) = " << nHat << endl;
        return nHat;
    };

    // Function returning nHat directed toward the metal using gradAlphaM
    auto getGradAlphaNormal = [this, &td](scalar cellI)
    {
        vector nHat = td.gradAlphaM(cellI);
        nHat.normalise();
        DebugPout << " --- nHat (gradAlphaM) = " << nHat << endl;

        if (mag(nHat) < SMALL)
        {
            Warning
                << "Particle #" << origId() << " in cell #" << cellI
                << ": mag(nHat) = " << mag(nHat) << ", alphaM = " << td.alphaM(cellI)
                << ", dQ = " << dQ_ << ", previous interaction type = " << type_ << endl;
        }
        return nHat;
    };

    // Function returning the path to the ray-interface intersection point
    // relative to the total path traversed by the particle inside the cell
    auto getPathToInterface = [&td](scalar cellI, const point& p0, const vector& dsv)
    {
        const vector C = td.surf().centre()[cellI];
        const vector normal = td.surf().normal()[cellI];
        const plane interfacePlane(C, normal, true); // true means to normalise the normal
        const scalar ds = mag(dsv);
        const scalar pathToInterface = interfacePlane.normalIntersect(p0, dsv/ds)/ds;
        DebugPout << " --- pathToInterface = " << pathToInterface << endl;
        return pathToInterface;
    };

    // Function returning the distance from the point to the interface
    auto getDistanceToInterface = [&td](scalar cellI, const point& p0)
    {
        const vector C = td.surf().centre()[cellI];
        vector nHat = td.surf().normal()[cellI];
        nHat.normalise();
        return (C - p0) & nHat;
    };

    // Function returning true if the entire traversed path goes inside the transparent medium
    auto pathGoesInsideGas = []
    (
        scalar cosTheta,
        scalar distanceToInterface,
        scalar pathToInterface
    )
    {
        // NB: SMALL is too small to capture grazing incidence
        return
            (cosTheta < -SMALL && pathToInterface <= 0) // outcoming && in the previous cell
         || (cosTheta > SMALL && pathToInterface >= 1) // incoming && in the next cell
         || (mag(cosTheta) < ROOTSMALL && distanceToInterface > 0); // parallel && out of metal
    };

    // Function returning true if the entire traversed path goes inside the absorbing medium
    auto pathGoesInsideMetal = []
    (
        scalar cosTheta,
        scalar distanceToInterface,
        scalar pathToInterface
    )
    {
        return
            (cosTheta < -SMALL && pathToInterface >= 1) // outcoming && in the next cell
         || (cosTheta > SMALL && pathToInterface <= 0) // incoming && in the previous cell
         || (mag(cosTheta) < ROOTSMALL && distanceToInterface <= 0); // parallel && out of gas
    };

    // --- Particle tracking loop
    do
    {
        const point p0 = position();
        const label cellI = cell();

        // Move the particle through the cell
        trackToFace(s, 1 - stepFraction());

        const vector dsv = position() - p0;
        const scalar ds = mag(dsv);

        if (iter == maxIter - 100)
        {
            SeriousError << "Too many iterations for particle #" << origId() << endl;
            debug = 1;
        }

        if (ds < smallLength)
        {
            DebugPout << " --- ds = 0, stepFraction = " << stepFraction() << endl;
            hitFace(s, cloud, td);
            continue;
        }

        const vector incidentDir = dsv/ds;

        // A. Do nothing if the cell consists of transparent medium only
        if (td.alphaM(cellI) < alphaTol)
        {
            if (isBeingAbsorbed_)
            {
                type_ = 1;
                DebugPout
                    << " +++ (1) particle to be absorbed enters the gas cell #" << cellI << endl;
                isBeingAbsorbed_ = false;
            }

            hitFace(s, cloud, td);
        }
        // B. Replace scattering by full absorption for low-energy particles
        else if (dQ_ < td.scattering().threshold())
        {
            type_ = 2;
            DebugPout << " +++ (2) low-energy particle" << endl;
            absorb(cellI, ds, 1);

            if (dQ_ < SMALL)
            {
                terminate();
                break;
            }

            hitFace(s, cloud, td);
        }
        // C. Absorb if the particle is marked as being absorbed and inside the absorbing medium
        else if (isBeingAbsorbed_)
        {
            // A particle to be absorbed always starts from the ray-interface intersection point
            // and penetrates the absorbing medium; therefore, absorb energy as much as possible.
            if (stepFraction()*mag(s) - ds < smallLength)
            {
                type_ = 3;
                DebugPout << " +++ (3) the first absorbing cell" << endl;
                absorb(cellI, ds, 1);
            }
            else
            {
                const vector nHat = getVoFNormal(cellI);

                // Absorb proportionally to the traversed path in the absorbing medium
                if (mag(nHat) > SMALL && td.useSubCellData())
                {
                    const scalar cosTheta = incidentDir & nHat;
                    const scalar distanceToInterface = getDistanceToInterface(cellI, p0);
                    const scalar pathToInterface = getPathToInterface(cellI, p0, dsv);

                    // Unmark the particle as being absorbed if it moves through the gas only
                    if (pathGoesInsideGas(cosTheta, distanceToInterface, pathToInterface))
                    {
                        type_ = 4;
                        DebugPout << " +++ (4) particle is no longer absorbed" << endl;
                        isBeingAbsorbed_ = false;
                    }
                    else
                    {
                        scalar pathInMetal = cosTheta < 0 ? pathToInterface : 1 - pathToInterface;

                        if (pathGoesInsideMetal(cosTheta, distanceToInterface, pathToInterface))
                        {
                            pathInMetal = 1;
                        }

                        type_ = 5;
                        DebugPout
                            << " +++ (5) absorbing cell with interface" << nl
                            << " --- pathInMetal = " << pathInMetal << endl;
                        absorb(cellI, ds, pathInMetal);
                    }
                }
                // Absorb proportionally to the metal fraction in the absence of subcell data
                else
                {
                    type_ = 6;
                    DebugPout << " +++ (6) absorbing cell without interface" << endl;
                    absorb(cellI, ds, 1);
                }
            }

            // Terminate the particle when all the energy is exhausted
            if (dQ_ < SMALL)
            {
                terminate();
                break;
            }

            hitFace(s, cloud, td);
        }
        // D. Scatter if the particle is not being absorbed and hit the gas-metal interface
        else
        {
            scalar transmissivity = 0;
            scalar pathToInterface = 0;
            vector nHat = getVoFNormal(cellI);
            scalar cosTheta = incidentDir & nHat; // positive for incoming rays

            // 1. Determine the ray-interface intersection point and the related quantities
            if (mag(nHat) > SMALL && td.useSubCellData())
            {
                const scalar distanceToInterface = getDistanceToInterface(cellI, p0);
                pathToInterface = getPathToInterface(cellI, p0, dsv);

                // Do not scatter rays that has not reached the absorbing medium
                if (pathGoesInsideGas(cosTheta, distanceToInterface, pathToInterface))
                {
                    type_ = 7;
                    DebugPout << " +++ (7) particle does not hit the interface" << endl;
                    hitFace(s, cloud, td);
                    continue;
                }
                // Scatter rays that enter the absorbing medium inside the cell
                else if (cosTheta > SMALL && pathToInterface > SMALL)
                {
                    type_ = 8;
                    DebugPout << " +++ (8) particle hit the interface inside the cell" << endl;
                }
                // Scatter rays that enter the absorbing medium at the face
                else
                {
                    type_ = 9;
                    DebugPout
                        << " +++ (9) particle penetrates into the metal at the face" << nl
                        << " --- alphaM = " << td.alphaM(cellI) << endl;
                    pathToInterface = 0;
                    // Recalculate the normal since the VoF normal cannot be used directly
                    nHat = getGradAlphaNormal(cellI);
                }
            }
            else
            {
                // Recalculate the normal since the VoF normal is absent
                nHat = getGradAlphaNormal(cellI);

                // Partially scatter rays that enter a cell that contains metal and no interface
                type_ = 10;
                DebugPout
                    << " +++ (10) particle scatters inside the cell without interface" << nl
                    << " --- alphaM = " << td.alphaM(cellI) << endl;
                transmissivity = 1 - td.alphaM(cellI);
                pathToInterface = 1 - td.alphaM(cellI);

                // Warn if the required subcell information is absent
                const scalar alphaTol2 = min(1e-2, 1e3*alphaTol);
                if (alphaTol2 < td.alphaM(cellI) && td.alphaM(cellI) < 1 - alphaTol2)
                {
                    td.markCell(cellI, trackingData::NO_INTERFACE);
                }
            }

            if (mag(nHat) < SMALL)
            {
                // Use normal incidence in case of extraordinary situations
                cosTheta = 1;
                FatalError << "Extraordinary event! Please report about it!" << exit(FatalError);
            }
            else
            {
                // Update the nHat-dependent quantities
                cosTheta = incidentDir & nHat;
            }

            // Do not scatter outcoming rays
            if (cosTheta < SMALL)
            {
                DebugPout
                    << " --- cosTheta = " << cosTheta << nl
                    << " +++ particle is outcoming" << endl;
                hitFace(s, cloud, td);
                continue;
            }

            // Do not scatter grazing incidence
            if (mag(cosTheta) < ROOTSMALL)
            {
                DebugPout
                    << " --- cosTheta = " << cosTheta << nl
                    << " +++ particle is grazing along the interface" << endl;
                hitFace(s, cloud, td);
                continue;
            }

            const point intersectionP = p0 + pathToInterface*dsv;

            DebugPout
                << "Particle #" << origId() << " is ready for scattering in cell #" << cellI << nl
                << " --- stepFraction = " << stepFraction() << nl
                << " --- alphaM = " << td.alphaM(cellI) << nl
                << " --- cosTheta = " << cosTheta << nl
                << " --- incidentDir = " << incidentDir << nl
                << " --- p0 = " << p0 << nl
                << " --- intersectionP = " << intersectionP << nl
                << " --- p1 = " << position() << endl;

            // 2. Compute the reflectivity and absorptivity coefficients
            const scalar R = (1 - transmissivity)*td.scattering().R(cosTheta);
            const scalar A = (1 - transmissivity)*(1 - R);
            DebugPout << " --- reflectivity = " << R << ", absorptivity = " << A << endl;

            // 3. Add new particles to the cloud and redistribute the energy to them
            addParticle(cellI, R, intersectionP, incidentDir, nHat, false);
            addParticle(cellI, A, intersectionP, incidentDir, nHat, true);
            dQ_ *= transmissivity;
            td.markCell(cellI, trackingData::SCATTERING);

            // 4. Transmit if the energy is not exhausted otherwise terminate
            if (transmissivity > SMALL)
            {
                DebugPout
                    << " --- transmissivity = " << transmissivity << nl
                    << " +++ particle transmits " << dQ_ << " of energy further" << endl;
                hitFace(s, cloud, td);
            }
            else
            {
                // Return particle to the intersection point before
                trackToFace(intersectionP - position(), 1);
                terminate();
                break;
            }
        }

    } while (!td.switchProcessor && stepFraction() < 1 && ++iter < maxIter);

    if (iter == maxIter)
    {
        FatalError << "Too many iterations for particle #" << origId() << exit(FatalError);
    }

    return td.keepParticle;
}


void Foam::rayTracingParticle::hitProcessorPatch(Cloud<rayTracingParticle>&, trackingData& td)
{
    td.switchProcessor = true;

    DebugPout << " +++ particle hit a processor patch #" << patch() << endl;
}


void Foam::rayTracingParticle::hitWallPatch(Cloud<rayTracingParticle>&, trackingData& td)
{
    stepFraction() = 1;

    DebugPout
        << " +++ particle hit a wall patch with " << dQ_ << " of energy" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rayTracingParticle::trackingData::~trackingData()
{
    label nAbsorbing = 0;
    label nScattering = 0;

    forAllConstIters(markedCells_, iter)
    {
        const label cellI = iter.key();
        const bitSet& bits = iter.val();

        bits.test(ABSORBING) && nAbsorbing++;
        bits.test(SCATTERING) && nScattering++;

        if (bits.test(NO_INTERFACE))
        {
            Warning
                << "There is no reconstructed interface in cell #" << cellI
                << ", where alphaM = " << alphaM(cellI) << endl;
        }
    }

    DebugPout
        << "Cells: " << nAbsorbing << " is absorbing, " << nScattering << " is scattering" << endl;
}


// ************************************************************************* //
