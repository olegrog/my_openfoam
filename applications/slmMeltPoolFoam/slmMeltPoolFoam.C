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
            slmMeltPoolFoam | Copyright (C) 2019-2020 Oleg Rogozin
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
    Solver for thermo-fluid-dynamic model of SLM based on interIsoFoam.

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

#include "gasMetalMixture.H"

volScalarField surfaceGaussian
(
    const volVectorField& x,
    dimensionedVector x0,
    dimensionedScalar radius
)
{
    using constant::mathematical::pi;
    volVectorField r = x - x0;
    r.replace(2, 0);
    return 2*exp(-2*magSqr(r/radius))/pi/sqr(radius);
}


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

    using constant::mathematical::pi;
    using constant::physicoChemical::R;     // universal gas constant
    using constant::physicoChemical::sigma; // Stefan--Boltzmann constant

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        const dimensionedVector laserCoordinate
        (
            "laserCoordinate", coordStart + laserVelocity*runTime
        );
        const dimensionedScalar totalEnthalpy = fvc::domainIntegrate(rho*h);
        // vector beamDirection(0, 0, -1);
        // mag(fvc::grad(alpha2) & beamDirection) can be used instead of mag(fvc::grad(alpha2))

        // --- Enthalpy--pressure--velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }


                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            // Enthalpy-independent quantities
            laserHeatSource = (runTime < timeStop)*absorptivity*laserPower
                *surfaceGaussian(mesh.C(), laserCoordinate, laserRadius);

            label nCorrEnthalpy(readLabel(pimple.dict().lookup("nEnthalpyCorrectors")));
            for (label corrEnthalpy = 1; corrEnthalpy <= nCorrEnthalpy; ++corrEnthalpy)
            {
                #include "hEqn.H"
            }

            rhok = rho*(1.0 - beta*(T - ambientTemperature));
            vaporPressure = ambientPressure*exp(molarMass*Hvapour/R*(1/Tboiling - 1/T));

            // For debug: check that heatConduction is almost equal to heatConduction2
            heatConvection = fvc::div(rhoPhi, h);
            heatConduction = fvc::laplacian(k, T);
            heatConduction2 = fvc::laplacian(k/Cp*(1 - Hfusion*liquidFraction.dEnthalpy()), h)
                + fvc::laplacian(k/Cp*hPrimeGasFraction, alpha1);

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

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

        // Update passive quantities
        mixture.finalCorrect();

        Info<< "Real energy input = " << (fvc::domainIntegrate(rho*h) - totalEnthalpy).value()
            << ", laser input = "
            << fvc::domainIntegrate(laserHeatSource*mag(fvc::grad(alpha2))*runTime.deltaT()).value()
            << ", theoretical value = " << (absorptivity*laserPower*runTime.deltaT()).value()
            << endl;

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
