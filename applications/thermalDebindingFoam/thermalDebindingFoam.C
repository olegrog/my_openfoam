/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020 Oleg Rogozin
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
    thermalDebindingFoam

Description
    Solver for the thermal debinding model of ceramic--polymer mixture.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"
#include "zeroGradientFvPatchField.H"

#include "multicomponentPolymer/multicomponentPolymer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run() && gMax(T) < maxTemperature)
    {
        #include "readTimeControls.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Calculate the temperature
        T += heatingRate*runTime.deltaT();

        // --- Calculate the burning of polymer components
        forAllIters(polymer.components(), iter)
        {
            polymerComponent& y = iter();

            fvScalarMatrix yEqn
            (
                fvm::ddt(y) + fvm::Sp(y.degradationRate(T), y)
            );

            yEqn.solve(mesh.solverDict("massFraction"));
        }

        // --- Solve the transport equation for monomer
        while (pimple.loop())
        {
            D1 = polymer.diffusion(rho, T);
            p = polymer.pressure(rho, T);
            permeability = permeability0
                *pow(polymer.poresFraction()/polymer.totalVolumeFraction(), particleSizeExponent)
                *pow(polymer.poresFraction(), 3)/sqr(1 - polymer.poresFraction());
            D2 == permeability/mu/polymer.poresFraction()*p;

            fvScalarMatrix rhoEqn
            (
                fvm::ddt(rho) == fvm::laplacian(D1, rho)
            );

            if (conservativeForm)
            {
                D2.correctBoundaryConditions();
                D2.boundaryFieldRef() /= 2; // rho = 0 is assumed at the boundary
                rhoEqn -= fvm::laplacian(D2, rho);
            }
            else
            {
                rhoEqn -= D2*fvm::laplacian(rho) + (fvc::grad(D2) & fvc::grad(rho));
            }

            forAllIters(polymer.components(), iter)
            {
                rhoEqn += polymer.initialVolumeFraction()*polymer.rho()*fvc::ddt(iter());
            }

            rhoEqn.solve();
        }

        porosity = polymer.poresFraction();
        Info<< "max(D1) = " << gMax(D1) << ", max(D2) = " << gMax(D2) << endl;

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
