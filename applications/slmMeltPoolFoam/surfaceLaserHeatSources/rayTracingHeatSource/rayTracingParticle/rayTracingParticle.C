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
    const scalar I,
    const label cellI,
    const scalar dA,
    const bool isTransmissive
)
:
    particle(mesh, position, cellI),
    p0_(position),
    p1_(targetPosition),
    I0_(I),
    I_(I),
    dA_(dA),
    isTransmissive_(isTransmissive)
{}


Foam::rayTracingParticle::rayTracingParticle
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
)
:
    particle(mesh, coordinates, celli, tetFacei, tetPti),
    p0_(position),
    p1_(targetPosition),
    I0_(I),
    I_(I),
    dA_(dA),
    isTransmissive_(isTransmissive)
{}


Foam::rayTracingParticle::rayTracingParticle(const rayTracingParticle& p)
:
    particle(p),
    p0_(p.p0_),
    p1_(p.p1_),
    I0_(p.I0_),
    I_(p.I_),
    dA_(p.dA_),
    isTransmissive_(p.isTransmissive_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::rayTracingParticle::move
(
    Cloud<rayTracingParticle>& spc,
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    do
    {
        // Cache old data of particle to use for reflected particle
        const point pos0 = position();
        const label cell1 = cell();

        scalar f = 1 - stepFraction();
        const vector s = p1_ - p0_;
        trackToAndHitFace(f*s, f, spc, td);

        const point p1 = position();
        vector dsv = p1 - pos0;
        scalar ds = mag(dsv);

        // NOTE:
        // Under the new barocentric tracking alghorithm the newly
        // inserted particles are tracked to the nearest cell centre first,
        // then, given the direction, to a face. In both occasions the first call
        // to trackToAndHitFace returns ds = 0. In this case we do an extra
        // call to trackToAndHitFace to start the tracking.
        // This is a temporary fix until the tracking can handle it.
        if (ds < SMALL*mesh().bounds().mag())
        {
            trackToAndHitFace(f*s, f, spc, td);
            dsv = p1 - position();
            ds = mag(dsv);
        }

        bool isReflectiveCell = td.reflectingCells()[cell1];
        DebugInfo<< " === " << p1[2] << ", isReflectiveCell_ = " << isReflectiveCell << endl;

        if (isReflectiveCell && !isTransmissive_)
        {
            // Create a new reflected particle when the particles is not
            // transmissive and larger than an absolute I

            scalar rho(0);

            if (I_ > 0.01*I0_ && ds > 0)
            {
                vector pDir = dsv/ds;

                cellPointWeight cpw(mesh(), position(), cell1, face());
                vector nHat = -td.gradAlphaMInterp().interpolate(cpw);

                nHat /= mag(nHat) + ROOTSMALL;
                scalar cosTheta(-pDir & nHat);

                DebugInfo<< " --- Ready for reflection: cosTheta = " << cosTheta << ' '
                    << pDir << ' ' << nHat << endl;

                // Only new incoming rays
                if (cosTheta > SMALL)
                {
                    vector newDir = td.reflection().R(pDir, nHat);
                    DebugInfo<< " --- newDir = " << newDir << endl;

                    // reflectivity
                    rho = min(max(td.reflection().rho(cosTheta), 0.0), 0.98);

                    scalar deltaM = cbrt(mesh().cellVolumes()[cell1]);

                    const point insertP(position() + newDir*0.01*deltaM);
                    label cellI = mesh().findCell(insertP);

                    if (cellI > -1)
                    {
                        auto pPtr = autoPtr<rayTracingParticle>::New
                        (
                            mesh(),
                            insertP,
                            insertP + newDir*mesh().bounds().mag(),
                            I_*rho,
                            cellI,
                            dA_,
                            false
                        );

                        // Add to cloud
                        spc.addParticle(pPtr.release());
                    }
                }
            }

            // Change isTransmissive of the particle
            isTransmissive_ = true;
            I0_ = I0_*(1 - rho);
            I_ = I_*(1 - rho);
            DebugInfo<< "Reflection: " << p1[2] << ", rho = " << rho << endl;
        }

        if (isReflectiveCell)
        {
            scalar a = td.AInterp().interpolate(pos0, cell1);

            scalar deltaI = I0_*ds*a;
            td.Q(cell1) += deltaI*dA_;
            I_ -= deltaI;

            DebugInfo<< "Position: " << p1[2] << ", I0_ = " << I0_ << ", I_ = " << I_ << ", deltaI = " << deltaI
                << ", ds*a = " << ds*a << ", a = " << a << ", stepFraction = " << stepFraction() << endl;

            if (I_ <= 0.01*I0_)
            {
                stepFraction() = 1.0;
                DebugInfo<< "Stop!" << endl;
                break;
            }
        }

    } while (td.keepParticle && !td.switchProcessor && stepFraction() < 1);

    return td.keepParticle;
}


void Foam::rayTracingParticle::hitProcessorPatch
(
    Cloud<rayTracingParticle>&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::rayTracingParticle::hitWallPatch
(
    Cloud<rayTracingParticle>&,
    trackingData& td
)
{
    td.keepParticle = false;
}


bool Foam::rayTracingParticle::hitPatch
(
    Cloud<rayTracingParticle>&,
    trackingData& td
)
{
    return false;
}


// ************************************************************************* //
