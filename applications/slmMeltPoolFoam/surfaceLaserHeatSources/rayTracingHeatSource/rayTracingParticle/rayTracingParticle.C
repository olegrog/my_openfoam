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
    td.switchProcessor = false;
    td.keepParticle = true;
    reset(); // to set stepFraction = 0

    // Function for translating energy from the particle to the heat source in the cell
    auto absorb = [&td, this](scalar cellI, scalar absorptionFraction)
    {
        if (absorptionFraction < SMALL) return;

        scalar deltaQ = min(dQ0_*absorptionFraction, dQ_);
        td.Q(cellI) += deltaQ;
        dQ_ -= deltaQ;

        DebugInfo
            << " --- absorbed energy = " << deltaQ
            << ", absorptionFraction = " << absorptionFraction
            << ", alphaM = " << td.alphaM(cellI) << endl;
    };

    // Function returning the unit normal vector directed to the gas using VoF subcell information
    auto getVoFNormal = [&td](scalar cellI)
    {
        vector nHat = -td.surf().normal()[cellI];
        nHat.normalise();
        DebugInfo << " --- nHat (VoF) = " << nHat << endl;
        return nHat;
    };

    // Function returning the unit normal vector directed to the gas using gradAlphaM
    auto getGradAlphaNormal = [&td](scalar cellI)
    {
        vector nHat = -td.gradAlphaM(cellI);
        nHat.normalise();
        DebugInfo << " --- nHat (gradAlphaM) = " << nHat << endl;
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
        DebugInfo << " --- pathToInterface = " << pathToInterface << endl;
        return pathToInterface;
    };

    // Function returning true if the entire traversed path goes inside the transparent medium
    auto pathGoesInsideGas = [](scalar cosTheta, scalar pathToInterface)
    {
        return
            (cosTheta < -SMALL && pathToInterface < 0) // outcoming && in the previous cell
         || (cosTheta > SMALL && pathToInterface > 1); // incoming && in the next cell
    };

    // surfCellTol is used in plicRDF.C and gradAlpha.C
    const scalar alphaTol = td.surf().modelDict().getOrDefault<scalar>("surfCellTol", 1e-8);
    const scalar maxTrackLength = mesh().bounds().mag();
    const scalar smallLength = SMALL*maxTrackLength;
    const vector s = p1_ - p0_;

    DebugInfo
        << "Particle #" << origId() << " is moving into direction " << s/mag(s)
        << " and transmit " << dQ0_ << " of the laser power" << endl;

    do
    {
        // Change the cell without moving
        hitFace(s, cloud, td);

        const point p0 = position();
        const label cellI = cell();

        // Move the particle through the cell
        trackToFace(s, 1 - stepFraction());

        const vector dsv = position() - p0;
        const scalar ds = mag(dsv);

        if (ds < smallLength)
        {
            DebugInfo << " --- ds = 0, stepFraction = " << stepFraction() << endl;
            continue;
        }

        const vector incidentDir = dsv/ds;

        // A. Do nothing if the cell consists of transparent medium only
        if (td.alphaM(cellI) < alphaTol)
        {
            if (isBeingAbsorbed_)
            {
                DebugInfo
                    << " +++ (1) particle to be absorbed enters the gas cell #"
                    << cellI << endl;
                isBeingAbsorbed_ = false;
            }
            continue;
        }

        // B. Replace scattering by full absorption for low-energy particles
        if (dQ_ < td.scattering().threshold() && dQ_ < td.alphaM(cellI))
        {
            DebugInfo << " +++ (2) low-energy particle" << endl;
            absorb(cellI, 1);
            if (dQ_ < SMALL)
            {
                stepFraction() = 1;
                break;
            }
            else
            {
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
                DebugInfo << " +++ (3) the first absorbing cell" << endl;
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
                        DebugInfo << " +++ (4) particle is no longer absorbed" << endl;
                        isBeingAbsorbed_ = false;
                    }
                    else
                    {
                        const scalar pathToInterfaceInMetal =
                            cosTheta < 0 ? pathToInterface : 1 - pathToInterface;
                        DebugInfo
                            << " +++ (5) absorbing cell with interface" << nl
                            << " --- pathToInterfaceInMetal = " << pathToInterfaceInMetal << endl;
                        absorb(cellI, min(td.alphaM(cellI), pathToInterfaceInMetal));
                    }
                }
                // Absorb proportionally to the metal fraction in the absence of subcell data
                else
                {
                    DebugInfo << " +++ (6) absorbing cell without interface" << endl;
                    absorb(cellI, td.alphaM(cellI));
                }
            }

            // Terminate propagation when all energy is exhausted
            if (dQ_ < SMALL)
            {
                stepFraction() = 1;
                DebugInfo << " +++ all energy has been absorbed" << endl;
                break;
            }
        }
        // D. Scatter if the particle is not being absorbed and hit the gas-metal interface
        else
        {
            scalar transmissivity = 0;
            point intersectionP = p0;
            vector nHat = getVoFNormal(cellI);
            scalar cosTheta = -incidentDir & nHat; // positive for incoming rays

            // 1. Determine the ray-interface intersection point and the related quantities
            if (mag(nHat) > SMALL)
            {
                const scalar pathToInterface = getPathToInterface(cellI, p0, dsv);

                // Do not scatter rays that has not reached the absorbing medium
                if (pathGoesInsideGas(cosTheta, pathToInterface))
                {
                    DebugInfo << " +++ (7) particle does not hit the interface" << endl;
                    continue;
                }
                // Scatter rays that enter the absorbing medium inside the cell
                else if (cosTheta > SMALL && pathToInterface > SMALL)
                {
                    DebugInfo << " +++ (8) particle hit the interface inside the cell" << endl;
                    intersectionP += pathToInterface*dsv;
                }
                // Scatter rays that enter the absorbing medium at the face
                else
                {
                    nHat = getGradAlphaNormal(cellI);

                    DebugInfo
                        << " +++ (9) particle penetrates into the metal at the face" << nl
                        << " --- alphaM = " << td.alphaM(cellI) << endl;
                }
            }
            else
            {
                nHat = getGradAlphaNormal(cellI);

                // Scatter rays that enter the absorbing medium at the face
                if (td.alphaM(cellI) > 1 - alphaTol)
                {
                    DebugInfo << " +++ (10) particle is inside the fully metal cell" << endl;
                }
                // Partially scatter rays that enter a cell that contains metal and no interface
                else
                {
                    DebugInfo
                        << " +++ (11) particle is inside the cell that contains metal" << nl
                        << " --- alphaM = " << td.alphaM(cellI) << endl;
                    transmissivity = 1 - td.alphaM(cellI);
                    // Assume that the scattering occurs at the middle of the cell path
                    intersectionP += dsv/2;
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
                DebugInfo
                    << " --- cosTheta = " << cosTheta << nl
                    << " +++ particle is outcoming" << endl;
                continue;
            }

            DebugInfo
                << "Particle #" << origId() << " is ready for scattering:" << nl
                << " --- stepFraction = " << stepFraction() << nl
                << " --- alphaM = " << td.alphaM(cellI) << nl
                << " --- cosTheta = " << cosTheta << nl
                << " --- incidentDir = " << incidentDir << nl
                << " --- p0 = " << p0 << nl
                << " --- intersectionP = " << intersectionP << nl
                << " --- p1 = " << position() << nl
                << " --- transmissivity = " << transmissivity << endl;

            // 2. Compute the reflectivity and absorptivity coefficients
            const scalar R = (1 - transmissivity)*td.scattering().R(cosTheta);
            const scalar A = (1 - transmissivity)*(1 - R);
            DebugInfo << " --- reflectivity = " << R << ", absorptivity = " << A << endl;

            // 3. Add the reflected particle to the cloud
            if (R*dQ_ > td.scattering().threshold())
            {
                const vector reflectionDir = td.scattering().reflection(incidentDir, nHat);
                const point reflectedP = intersectionP + reflectionDir*maxTrackLength;
                auto reflectedParticlePtr = autoPtr<rayTracingParticle>::New
                (
                    mesh(), intersectionP, reflectedP, R*dQ_, cellI, false
                );
                DebugInfo
                    << " --- particle #" << reflectedParticlePtr->origId() << " is reflected"
                    << " in direction " << reflectionDir << endl;
                cloud.addParticle(reflectedParticlePtr.release());
            }

            // 4. Add the refracted particle to the cloud
            if (A*dQ_ > td.scattering().threshold())
            {
                const vector refractionDir = td.scattering().refraction(incidentDir, nHat);
                const point refractedP = intersectionP + refractionDir*maxTrackLength;
                auto refractedParticlePtr = autoPtr<rayTracingParticle>::New
                (
                    mesh(), intersectionP, refractedP, A*dQ_, cellI, true
                );
                DebugInfo
                    << " --- particle #" << refractedParticlePtr->origId() << " is refracted"
                    << " in direction " << refractionDir << endl;
                cloud.addParticle(refractedParticlePtr.release());
            }

            // 5. Update the power contribution and determine the further life of the particle
            dQ_ *= transmissivity;
            // Do nothing if energy is not exhausted
            if (dQ_ > td.scattering().threshold())
            {
                DebugInfo
                    << " +++ particle transmits " << dQ_ << " of laser power further" << endl;
            }
            // Return particle to the intersection point and terminate its propagation
            else
            {
                absorb(cellI, 1);
                trackToFace(intersectionP - position(), 1);
                stepFraction() = 1;
                break;
            }
        }

    } while (td.keepParticle && !td.switchProcessor && stepFraction() < 1);

    return td.keepParticle;
}


void Foam::rayTracingParticle::hitProcessorPatch(Cloud<rayTracingParticle>&, trackingData& td)
{
    td.switchProcessor = true;

    DebugInfo << " +++ particle hit a processor patch" << endl;
}


void Foam::rayTracingParticle::hitWallPatch(Cloud<rayTracingParticle>&, trackingData& td)
{
    stepFraction() = 1;

    DebugInfo
        << " +++ particle hit a wall patch with " << dQ_
        << " of the laser power" << endl;
}


// ************************************************************************* //
