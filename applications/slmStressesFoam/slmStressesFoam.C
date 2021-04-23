/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
            slmStressesFoam | Copyright (C) 2020-2021 Oleg Rogozin
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
    slmStressesFoam

Description
    Solver for thermo-mechanical model of selective laser melting (SLM)
    based on solidDisplacementFoam.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicRefineFvMesh.H"
#include "pimpleControl.H"

#include "quiescentGasMetalMixture.H"
#include "volumetricLaserHeatSource.H"

volScalarField surfaceGaussian
(
    const volVectorField& x,
    dimensionedVector x0,
    dimensionedScalar radius
)
{
    using constant::mathematical::pi;
    volVectorField r = x - x0;
    r.replace(vector::Z, 0);
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
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    #include "laserCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "laserCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Calculate time-dependent quantities
        const dimensionedScalar totalEnthalpy = fvc::domainIntegrate(rho*h);

        // --- Alpha-enthalpy-pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter())
            {
                mesh.update();

                if (mesh.changing())
                {
                    mixture.correct();
                }

                laserHeatSource->correct();
            }

            #include "hEqn.H"
        }

        wasMelted = Foam::max(wasMelted, liquidFraction);

        #include "effectiveAbsorptivity.H"

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
