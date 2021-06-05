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
    isBeingAbsorbed_(isBeingAbsorbed)
{}


Foam::rayTracingParticle::rayTracingParticle(const rayTracingParticle& p)
:
    particle(p),
    p0_(p.p0_),
    p1_(p.p1_),
    dQ0_(p.dQ0_),
    dQ_(p.dQ_),
    isBeingAbsorbed_(p.isBeingAbsorbed_)
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
    auto absorb = [&td, this](scalar cellI, scalar absorptionFraction)
    {
        scalar deltaQ = min(dQ0_*absorptionFraction, dQ_);
        if (deltaQ < SMALL) return;

        td.Q(cellI) += deltaQ;
        dQ_ -= deltaQ;

        DebugPout
            << " --- absorbed energy = " << deltaQ
            << ", absorptionFraction = " << absorptionFraction
            << ", alphaM = " << td.alphaM(cellI) << endl;
    };

    // Function for stopping propagation of the particle
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
        if (energyFraction > td.scattering().threshold())
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
        }

        return energyFraction*dQ_;
    };

    // Function returning the unit normal vector directed to the gas using VoF subcell information
    auto getVoFNormal = [&td](scalar cellI)
    {
        vector nHat = -td.surf().normal()[cellI];
        nHat.normalise();
        DebugPout << " --- nHat (VoF) = " << nHat << endl;
        return nHat;
    };

    // Function returning the unit normal vector directed to the gas using gradAlphaM
    auto getGradAlphaNormal = [&td](scalar cellI)
    {
        vector nHat = -td.gradAlphaM(cellI);
        nHat.normalise();
        DebugPout << " --- nHat (gradAlphaM) = " << nHat << endl;
        return nHat;
    };

    // Function returning the distance to the ray-interface intersection point in terms of
    // the path traversed by the particle inside the cell
    auto getPathToInterface = [&td](scalar cellI, const point& p0, const vector& dsv)
    {
        vector C = td.surf().centre()[cellI];
        vector normal = -td.surf().normal()[cellI];
        plane interfacePlane(C, normal, true); // true means to normalise the normal
        scalar ds = mag(dsv);
        scalar pathToInterface = interfacePlane.normalIntersect(p0, dsv/ds)/ds;
        DebugPout << " --- pathToInterface = " << pathToInterface << endl;
        return pathToInterface;
    };

    // Function returning true if the entire traversed path goes inside the transparent medium
    auto pathGoesInsideGas = [](scalar cosTheta, scalar pathToInterface)
    {
        return
            (cosTheta < -SMALL && pathToInterface < 0) // outcoming && in the previous cell
         || (cosTheta > SMALL && pathToInterface > 1); // incoming && in the next cell
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
                DebugPout
                    << " +++ (1) particle to be absorbed enters the gas cell #" << cellI << endl;
                isBeingAbsorbed_ = false;
            }

            hitFace(s, cloud, td);
            continue;
        }

        // B. Replace scattering by full absorption for low-energy particles
        if (dQ_ < td.scattering().threshold() && dQ_ < td.alphaM(cellI))
        {
            DebugPout << " +++ (2) low-energy particle" << endl;
            absorb(cellI, 1);

            if (dQ_ < SMALL)
            {
                terminate();
                break;
            }
            else
            {
                hitFace(s, cloud, td);
                continue;
            }
        }

        // C. Absorb if the particle is marked as being absorbed and inside an absorbing medium
        if (isBeingAbsorbed_)
        {
            // A particle to be absorbed always starts from the ray-interface intersection point
            // and penetrates the absorbing medium; therefore, absorb energy as much as possible.
            if (stepFraction()*mag(s) - ds < smallLength)
            {
                DebugPout << " +++ (3) the first absorbing cell" << endl;
                const scalar cellSize = cbrt(mesh().cellVolumes()[cellI]);
                absorb(cellI, min(td.alphaM(cellI), ds/cellSize));
            }
            else
            {
                const vector nHat = getVoFNormal(cellI);

                // Absorb proportionally to the traversed path in the absorbing medium
                if (mag(nHat) > SMALL)
                {
                    const scalar cosTheta = -incidentDir & nHat;
                    const scalar pathToInterface = getPathToInterface(cellI, p0, dsv);

                    // Unmark the particle as being absorbed if it propagates through the gas only
                    if (pathGoesInsideGas(cosTheta, pathToInterface))
                    {
                        DebugPout << " +++ (4) particle is no longer absorbed" << endl;
                        isBeingAbsorbed_ = false;
                    }
                    else
                    {
                        const scalar pathToInterfaceInMetal =
                            cosTheta < 0 ? pathToInterface : 1 - pathToInterface;
                        DebugPout
                            << " +++ (5) absorbing cell with interface" << nl
                            << " --- pathToInterfaceInMetal = " << pathToInterfaceInMetal << endl;
                        absorb(cellI, min(td.alphaM(cellI), pathToInterfaceInMetal));
                    }
                }
                // Absorb proportionally to the metal fraction in the absence of subcell data
                else
                {
                    DebugPout << " +++ (6) absorbing cell without interface" << endl;
                    absorb(cellI, td.alphaM(cellI));
                }
            }

            // Terminate propagation when all energy is exhausted
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
            scalar deltaQ = 0;
            scalar transmissivity = 0;
            scalar pathToInterface = 0;
            vector nHat = getVoFNormal(cellI);
            scalar cosTheta = -incidentDir & nHat; // positive for incoming rays

            // 1. Determine the ray-interface intersection point and the related quantities
            if (mag(nHat) > SMALL)
            {
                pathToInterface = getPathToInterface(cellI, p0, dsv);

                // Do not scatter rays that has not reached the absorbing medium
                if (pathGoesInsideGas(cosTheta, pathToInterface))
                {
                    DebugPout << " +++ (7) particle does not hit the interface" << endl;
                    hitFace(s, cloud, td);
                    continue;
                }
                // Scatter rays that enter the absorbing medium inside the cell
                else if (cosTheta > SMALL && pathToInterface > SMALL)
                {
                    DebugPout << " +++ (8) particle hit the interface inside the cell" << endl;
                }
                // Scatter rays that enter the absorbing medium at the face
                else
                {
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

                // Scatter rays that enter the absorbing medium at the face
                if (td.alphaM(cellI) > 1 - alphaTol)
                {
                    DebugPout << " +++ (10) particle is inside the fully metal cell" << endl;
                }
                // Partially scatter rays that enter a cell that contains metal and no interface
                else
                {
                    DebugPout
                        << " +++ (11) particle is inside the cell that contains metal" << nl
                        << " --- alphaM = " << td.alphaM(cellI) << endl;
                    transmissivity = 1 - td.alphaM(cellI);
                    // Assume that the scattering occurs at the middle of the cell path
                    pathToInterface = 0.5;
                }

                // Warn if the required subcell information is absent
                const scalar alphaTol2 = min(1e-2, 1e3*alphaTol);
                if (alphaTol2 < td.alphaM(cellI) && td.alphaM(cellI) < 1 - alphaTol2)
                {
                    Warning
                        << "There is no reconstructed interface in cell #" << cellI
                        << ", where alphaM = " << td.alphaM(cellI) << endl;
                }
            }

            // Crash if the interface normal failed to be determined
            if (mag(nHat) < SMALL)
            {
                FatalError
                    << "Particle #" << origId() << " in cell #" << cellI << nl
                    << "mag(nHat) = " << mag(nHat) << ", alphaM = " << td.alphaM(cellI)
                    << exit(FatalError);
            }

            // Update the nHat-dependent quantities
            cosTheta = -incidentDir & nHat;

            // Do not scatter outcoming rays
            if (cosTheta < SMALL)
            {
                DebugPout
                    << " --- cosTheta = " << cosTheta << nl
                    << " +++ particle is outcoming" << endl;
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

            if (transmissivity > SMALL)
            {
                DebugPout
                    << " --- transmissivity = " << transmissivity << endl;
            }

            // 2. Compute the reflectivity and absorptivity coefficients
            const scalar R = (1 - transmissivity)*td.scattering().R(cosTheta);
            const scalar A = (1 - transmissivity)*(1 - R);
            DebugPout << " --- reflectivity = " << R << ", absorptivity = " << A << endl;

            // 3. Add new particles to the cloud and redistribute the energy to them
            deltaQ += addParticle(cellI, R, intersectionP, incidentDir, nHat, false);
            deltaQ += addParticle(cellI, A, intersectionP, incidentDir, nHat, true);
            dQ_ -= deltaQ;

            // 4. Partially absorb the remaining energy
            absorb(cellI, min(1 - pathToInterface, td.alphaM(cellI)));

            // 5. Transmit if the energy is not exhausted
            if (dQ_ > SMALL)
            {
                DebugPout
                    << " +++ particle transmits " << dQ_ << " of energy further" << endl;

                isBeingAbsorbed_ = true;
                hitFace(s, cloud, td);
                continue;
            }

            // 6. Return particle to the intersection point and terminate its propagation
            trackToFace(intersectionP - position(), 1);
            terminate();
            break;
        }

    } while (!td.switchProcessor && stepFraction() < 1);

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


// ************************************************************************* //
