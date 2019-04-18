/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.2.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
-------------------------------------------------------------------------------
Application
    solidificationFoam

Description
    Solver for the phase-field equations for the solidification of binary alloys.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"

#include "multicomponentAlloy/multicomponentAlloy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    pimpleControl pimple(mesh);

    #include "createFields.H"

    const boundBox& bounds = mesh.bounds();
    const dimensionedScalar ymax("ymax", dimLength, bounds.max().y());
    const dimensionedScalar ymin("ymin", dimLength, bounds.min().y());
    phase = atan((coordY - (ymax + ymin) / 2) / interfaceWidth) / mathematicalConstant::pi + 0.5;

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        /** Temperature calculation */

        /** Phase field calculation */

        while (pimple.correct()) {
            volVectorField normal = fvc::grad(phase) / mag(fvc::grad(phase));
            volScalarField theta = 0 * phase;
            forAll(theta, cellI) {
                theta[cellI] = normal[cellI].x() / normal[cellI].y();
            }
            volScalarField a_s = 1 + epsilon4 * Foam::cos(4 * theta);

            volVectorField A = 0 * fvc::grad(phase);
            forAll(A, cellI) {
                A[cellI].x() = - 4 * epsilon4 * Foam::sin(4 * theta[cellI]) * Foam::cos(theta[cellI]);
                A[cellI].y() =   4 * epsilon4 * Foam::sin(4 * theta[cellI]) * Foam::sin(theta[cellI]);
            }

            tmp<volScalarField> f_prime = 2 * phase * (1 - 2 * sqr(phase));
            tmp<volScalarField> g_prime = 30 * sqr(phase) * sqr(1 - phase);
            fvScalarMatrix phaseEqn(
                tau * sqr(a_s) * fvm::ddt(phase)
                == sqr(interfaceWidth) * fvm::laplacian(sqr(a_s), phase)
                + sqr(interfaceWidth) * fvc::div(A) - f_prime
                + a1 * interfaceWidth / sigma * g_prime * alloy.chemicalDrivingForce(phase, T)
                    / molarVolume
            );
            phaseEqn.solve();
        }

        /** Concentrations calculation */

        volVectorField normal = fvc::grad(phase) / mag(fvc::grad(phase));
        volScalarField h = alloy.partition(phase);
        forAllIter(PtrDictionary<alloyComponent>, alloy.components(), iter)
        {
            alloyComponent& C = iter();

            surfaceScalarField phi = fvc::snGrad(phase) * mesh.magSf()
                * fvc::interpolate(C.diffusion(phase) / sqr(h)) * alloy.partition_prime();
            fvScalarMatrix CEqn(
                fvm::ddt(C) == fvm::laplacian(C.diffusion(phase) / h, C) - fvm::div(phi, C)
                - fvc::laplacian(C.diffusion(phase), C.equilibrium(phase, T) / h)
                + fvc::div(interfaceWidth / sqrt(2) * C.deltaA() * normal * fvc::ddt(phase))
            );
            CEqn.solve(mesh.solutionDict().solver("concentration"));
        }

        runTime.write();
        Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
