/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020-2020 Oleg Rogozin
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
#include "simpleControl.H"

#include "multicomponentPolymer/multicomponentPolymer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop() && gMax(T) < maxTemperature)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Calculate the temperature

        T += heatingRate*runTime.deltaT();

        // --- Calculate the burning of polymer components

        forAllIter(PtrDictionary<polymerComponent>, polymer.components(), iter)
        {
            polymerComponent& y = iter();

            fvScalarMatrix yEqn
            (
                fvm::ddt(y) + fvm::Sp(y.degradationRate(T), y)
            );

            yEqn.solve(mesh.solverDict("massFraction"));
        }

        // --- Solve the transport equation for monomer

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix rhoEqn
            (
                fvm::ddt(rho)
             ==
                fvm::laplacian(polymer.diffusion(rho, T), rho)
                // + convection term
            );

            forAllIter(PtrDictionary<polymerComponent>, polymer.components(), iter)
            {
                rhoEqn += polymer.initialVolumeFraction()*polymer.rho()*fvc::ddt(iter());
            }

            rhoEqn.solve();
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
