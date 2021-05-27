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
    const bool isTransmissive
)
:
    particle(mesh, position, cellI),
    p0_(position),
    p1_(targetPosition),
    dQ0_(dQ),
    dQ_(dQ),
    isTransmissive_(isTransmissive)
{}


Foam::rayTracingParticle::rayTracingParticle(const rayTracingParticle& p)
:
    particle(p),
    p0_(p.p0_),
    p1_(p.p1_),
    dQ0_(p.dQ0_),
    dQ_(p.dQ_),
    isTransmissive_(p.isTransmissive_)
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

    // surfCellTol is used in plicRDF.C and gradAlpha.C
    const scalar alphaTol = td.surf().modelDict().getOrDefault<scalar>("surfCellTol", 1e-8);
    const scalar maxTrackLength = mesh().bounds().mag();
    const scalar smallLength = SMALL*maxTrackLength;
    const vector s = p1_ - p0_;

    DebugInfo
        << "Particle #" << origId() << " is moving into direction " << s/mag(s) << endl;

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

        // If the particle have propagated through a transparent medium, then continue
        if (td.alphaM(cellI) < alphaTol) continue;

        // If the particle is transmissive and inside an absorbing cell, then absorb
        // proportionally to its traversed path in the cell
        if (isTransmissive_)
        {
            const scalar cellDelta = cbrt(mesh().cellVolumes()[cellI]);
            const scalar deltaQ = min(dQ0_*ds/cellDelta, dQ_);

            td.Q(cellI) += deltaQ;
            dQ_ -= deltaQ;

            DebugInfo
                << " --- absorbed energy = " << deltaQ
                << ", alphaM = " << td.alphaM(cellI) << endl;

            // Terminate propagation when all energy is exhausted
            if (dQ_ <= SMALL)
            {
                stepFraction() = 1;
                DebugInfo << " +++ all energy has been absorbed" << endl;
                break;
            }
        }
        // If the particle is not transmissive and hit the reconstructed interface, then scatter
        else
        {
            bool useGradAlpha = true;
            const vector incidentDir = dsv/ds;
            point intersectionP = p0;
            vector nHat = -td.surf().normal()[cellI]; // directed to the gas
            nHat.normalise();
            scalar cosTheta = -incidentDir & nHat;

            if (mag(nHat) > ROOTVSMALL)
            {
                DebugInfo << " --- normal = " << nHat << endl;

                // Do not scatter outcoming rays
                if (cosTheta < SMALL)
                {
                    DebugInfo
                        << " --- alphaM = " << td.alphaM(cellI) << nl
                        << " --- cosTheta = " << cosTheta << nl
                        << " +++ particle is outcoming" << endl;
                    continue;
                }

                // Use subcell VoF information to determine the ray-interface intersection point
                const vector C = td.surf().centre()[cellI];
                const plane interfacePlane(C, nHat, false);
                const scalar pathToInterface = interfacePlane.normalIntersect(p0, incidentDir);

                DebugInfo
                    << " --- centre = " << C << nl
                    << " --- pathToInterface = " << pathToInterface
                    << ", relative = " << pathToInterface/ds << endl;

                if (pathToInterface >= ds)
                {
                    DebugInfo << " +++ particle does not hit the interface inside the cell" << endl;
                    continue;
                }

                if (pathToInterface >= 0)
                {
                    intersectionP += pathToInterface*incidentDir;
                    useGradAlpha = false;
                }
            }
            else
            {
                // Value of gradAlphaM is not enough to identify the intersection point.
                // Therefore, we try to use it only when alphaM = 1 and alphaM = 0 or when
                // the intersection point is behind the cellI (pathToInterface < 0).
                const scalar alphaTol2 = min(1e-3, 1e2*alphaTol);

                // Skip scattering if this cell contains only small amount of metal
                if (td.alphaM(cellI) < alphaTol2)
                {
                    DebugInfo
                        << " +++ scattering is skipped since alphaM = "
                        << td.alphaM(cellI) << endl;
                    continue;
                }

                if (alphaTol2 < td.alphaM(cellI) && td.alphaM(cellI) < 1 - alphaTol2)
                {
                    // Sometimes this case corresponds to a droplet that is smaller than the cell
                    Warning
                        << "There is no reconstructed interface in cell #" << cellI
                        << ", where alphaM = " << td.alphaM(cellI) << endl;
                }

                // We use this crude approximation to ensure that the heat source is redistributed
                // according to the metal fraction
                intersectionP += (1 - td.alphaM(cellI))*dsv;
            }

            if (useGradAlpha)
            {
                DebugInfo << " +++ gradAlphaM is used to reconstruct nHat" << endl;

                nHat = -td.gradAlphaM(cellI);
                nHat.normalise();
                cosTheta = -incidentDir & nHat;
            }

            if (mag(nHat) < SMALL)
            {
                Warning
                    << "Particle #" << origId() << " in cell #" << cellI
                    << ": mag(nHat) = " << mag(nHat) << ", alphaM = " << td.alphaM(cellI) << endl;
                continue;
            }

            DebugInfo
                << "Particle #" << origId() << " is ready for scattering:" << nl
                << " --- stepFraction = " << stepFraction() << nl
                << " --- alphaM = " << td.alphaM(cellI) << nl
                << " --- dQ = " << dQ_ << nl
                << " --- cosTheta = " << cosTheta << nl
                << " --- normal = " << nHat << nl
                << " --- incidentDir = " << incidentDir << nl
                << " --- p0 = " << p0 << nl
                << " --- intersectionP = " << intersectionP << nl
                << " --- p1 = " << position() << endl;

            if (cosTheta < SMALL)
            {
                DebugInfo << " +++ particle is outcoming (gradAlphaM is used)" << endl;
                continue;
            }

            // Replace scattering by direct absorption for low-energy particles
            if (dQ_ < td.scattering().threshold() && dQ_ < td.alphaM(cellI))
            {
                td.Q(cellI) += dQ_;

                DebugInfo
                    << " --- absorbed energy " << dQ_ << ", alphaM = " << td.alphaM(cellI) << endl;

                stepFraction() = 1;
                break;
            }

            const scalar rho = td.scattering().rho(cosTheta);
            DebugInfo << " --- reflectivity = " << rho << endl;

            if (rho > SMALL)
            {
                // 1. Add the reflected particle to the cloud
                const vector reflectionDir = td.scattering().reflection(incidentDir, nHat);
                const point reflectedP = intersectionP + reflectionDir*maxTrackLength;
                auto reflectedParticlePtr = autoPtr<rayTracingParticle>::New
                (
                    mesh(), intersectionP, reflectedP, dQ_*rho, cellI, false
                );
                DebugInfo
                    << " --- particle #" << reflectedParticlePtr->origId() << " is reflected"
                    << " in direction " << reflectionDir << endl;
                cloud.addParticle(reflectedParticlePtr.release());
            }

            // 2. Add the refracted particle to the cloud
            const vector refractionDir = td.scattering().refraction(incidentDir, nHat);
            const point refractedP = intersectionP + refractionDir*maxTrackLength;
            auto refractedParticlePtr = autoPtr<rayTracingParticle>::New
            (
                mesh(), intersectionP, refractedP, dQ_*(1 - rho), cellI, true
            );
            DebugInfo
                << " --- particle #" << refractedParticlePtr->origId() << " is refracted"
                << " in direction " << refractionDir << endl;
            cloud.addParticle(refractedParticlePtr.release());

            // 3. Return particle to the intersection point and terminate its propagation
            trackToFace(intersectionP - position(), 1);
            stepFraction() = 1;
            break;
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

    DebugInfo << " +++ particle hit a wall patch" << endl;
}


// ************************************************************************* //
