/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019-2021 Oleg Rogozin
-------------------------------------------------------------------------------
License
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

Application
    solidificationFoam

Description
    Solver for the multicomponent phase-field model with linearised phase
    diagrams, convection, and tracking of the grain boundaries.
    The directional solidification problem is assumed to be solved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "zeroGradientFvPatchField.H"

#include "multicomponentAlloy/multicomponentAlloy.H"

using constant::mathematical::pi;
using constant::mathematical::twoPi;

#include "grainGenerationFunctions.H"

tmp<volVectorField> calcNormal(const volVectorField& vField)
{
    dimensionedVector smallVector("small", vField.dimensions(), vector(0, SMALL, 0));
    return (vField + smallVector)/mag(vField + smallVector);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // --- Additional variables

    volScalarField theta = 0*phase;
    volScalarField theta0 = 0*phase;
    const tensor rot(0, -1, 0, 1, 0, 0, 0, 0, 1);

    // Derived quantities
    dimensionedScalar pullingSpeed = coolingRate/tempGradient;
    const label nGrains = crystallographicAngles.size();

    // Coordinates-related constants and variables
    const volVectorField coord = mesh.C();
    const boundBox& bounds = mesh.bounds();
    const dimensionedScalar xmax("xmax", dimLength, bounds.max().x());
    const dimensionedScalar xmin("xmin", dimLength, bounds.min().x());
    const dimensionedScalar ymax("ymax", dimLength, bounds.max().y());
    const dimensionedScalar ymin("ymin", dimLength, bounds.min().y());
    const dimensionedScalar height("height", ymax - ymin);
    const dimensionedScalar width("width", xmax - xmin);
    const dimensionedVector center("center", dimLength, (bounds.max() + bounds.min())/2);
    const dimensionedScalar frontPosition = ymin + frontPositionRel*height;
    const dimensionedScalar initialWidth = interfaceWidth/interfaceNarrowing;
    dimensionedScalar tipSpeed = pullingSpeed;
    dimensionedScalar tipPosition = frontPosition;
    dimensionedScalar tipPositionPrev = tipPosition;
    label tipCell(0), tipCellPrev(0);
    dimensionedScalar tipTimeInterval("tipTimeInterval", dimTime, 0);

    #include "initialConditions.H"
    #include "printProblemParameters.H"
    #include "doubleWellPotentialFunctions.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "phaseCourantNo.H"
        #include "setDeltaT.H"

        scalar minDeltaT =
            runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", 0);
        if (runTime.deltaTValue() < minDeltaT)
        {
            Info<< "Time step becomes too small!" << endl;
            runTime.writeAndEnd();
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "tipPosition.H"

        // --- Calculate temperature

        T = alloy.liquidus() - undercooling
            + tempGradient*(coord.component(vector::Y) - frontPosition)
            - coolingRate*runTime;

        tau = a1*a2*pow3(interfaceWidth)/alloy.interfaceEnergy()*alloy.sumRestrictionFactors(T);

        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    if (correctPhi)
                    {
                        #include "correctPhi.H"
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            label nCorrPhase(readLabel(pimple.dict().lookup("nPhaseCorrectors")));
            for (label corrPhase = 1; corrPhase <= nCorrPhase; ++corrPhase)
            {
                #include "phaseEqn.H"
            }

            scalar phaseTolerance = mesh.solverDict(phase.name()).get<scalar>("tolerance");
            if (gMin(phase) < -phaseTolerance || gMax(phase) > 1 + phaseTolerance)
            {
                FatalError
                    << "Phase is out of bounds."
                    << abort(FatalError);
            }

            phase.clip(0, 1);

            if (!pimple.frozenFlow())
            {
                #include "UEqn.H"

                while (pimple.correct())
                {
                    #include "pEqn.H"
                }
            }

            #include "CEqn.H"
        }

        #include "grainEqn.H"

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
