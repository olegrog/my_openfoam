/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                  interFoam | Copyright (C) 2011-2017 OpenFOAM Foundation
                isoAdvector | Copyright (C) 2016 DHI
              Modified work | Copyright (C) 2018 Johan Roenby
               interIsoFoam | Copyright (C) 2019-2020 DLR
            slmMeltPoolFoam | Copyright (C) 2019-2023 Oleg Rogozin
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
    slmMeltPoolFoam

Description
    Solver for thermo-fluid-dynamic model of selective laser melting (SLM)
    based on interIsoFoam.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "dynamicRefineFvMesh.H"
#include "IOmanip.H"

#include "incompressibleGasMetalMixture.H"
#include "surfaceLaserHeatSource.H"
#include "movingReferenceFrame.H"

// For debug
auto pp = [](const volScalarField& f)
{
    Info<< f.name() << ": " << min(f).value() << " " << max(f).value() << endl;
    // The following methods works only for internal field but globally
    // f.writeMinMax(Info);
    // Info<< f.name() << ": " << gMin(f) << " " << gMax(f) << endl;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "temperatureChange.H"
        #include "setDeltaT.H"

        ++runTime;
        bool reduceTimeStep = false;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Calculate time-dependent quantities
        const dimensionedScalar totalEnthalpy = fvc::domainIntegrate(rho*h);
        // NB: SMALL is too small to be used in the denominator
        const volScalarField divUInMetal("divUInMetal", fvc::div(phi)/(alpha1 + ROOTSMALL));

        // --- Alpha-enthalpy-pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.corr() == pimple.nCorrPIMPLE())
            {
                reduceTimeStep = true;
            }

            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                if (isA<dynamicRefineFvMesh>(mesh))
                {
                    advector.surf().reconstruct();
                }

                const scalar meshUpdateStartTime = runTime.elapsedCpuTime();

                mesh.update();

                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    if (isA<dynamicRefineFvMesh>(mesh))
                    {
                        advector.surf().mapAlphaField();
                        alpha2 = 1.0 - alpha1;
                        alpha2.correctBoundaryConditions();
                        rho == alpha1*rho1 + alpha2*rho2;
                        rho.correctBoundaryConditions();
                        rho.oldTime() = rho;
                        alpha2.oldTime() = alpha2;

                        // We have to update the alpha-dependent fields
                        mixture.correct();
                        mixture.correctThermo();

                        // Need for ddt in correctPassiveFields() and mixture.divPhi()
                        const std::vector<std::reference_wrapper<const volScalarField>> fields
                        {
                            std::cref(mixture.liquidFractionInMetal()),
                            std::cref(liquidFraction),
                            std::cref(T)
                        };
                        for (auto& field : fields)
                        {
                            const_cast<volScalarField&>(field.get()).oldTime() = field;
                        }
                    }

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        // Use velocity divergence from the old mesh
                        divU = alpha1*divUInMetal;

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        // NB: the following line is uncommented in interIsoFoam
                        // mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }

                meshUpdateTime = runTime.elapsedCpuTime() - meshUpdateStartTime;

                // Moving reference frame
                if (MRF->moving())
                {
                    MRF->correct();
                    phi -= MRF->phiRel();
                }

                // Update position of the laser beam
                laserHeatSource->update();
            }

            // --- Advect alpha field
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();  // update mixture.gradAlpha needed for laserHeatSource
            laserHeatSource->correct();

            // --- Momentum predictor
            #include "UEqn.H"

            // --- Enthalpy corrector loop
            label nCorrEnthalpy(readLabel(pimple.dict().lookup("nEnthalpyCorrectors")));
            for (label corrEnthalpy = 1; corrEnthalpy <= nCorrEnthalpy; ++corrEnthalpy)
            {
                #include "hEqn.H"
            }

            if (pimple.frozenFlow())
            {
                continue;
            }

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        #include "updatePassiveFields.H"
        #include "effectiveAbsorptivity.H"

        runTime.write();
        #include "timeConsumption.H"
        runTime.printExecutionTime(Info);

        if (reduceTimeStep)
        {
            Info<< "Halve the time step since all PIMPLE iterations were used" << endl;
            runTime.setDeltaT(runTime.deltaTValue()/2, false);
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
