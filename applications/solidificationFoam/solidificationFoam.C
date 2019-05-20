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
    const dimensionedScalar xmax("xmax", dimLength, bounds.max().x());
    const dimensionedScalar xmin("xmin", dimLength, bounds.min().x());
    const dimensionedScalar ymax("ymax", dimLength, bounds.max().y());
    const dimensionedScalar ymin("ymin", dimLength, bounds.min().y());
    const dimensionedScalar front = (ymin + ymax*frontPosition) / (1 + frontPosition);
    phase = Foam::atan((front - coordY) / 2 / interfaceWidth) / mathematicalConstant::pi + 0.5;
    phase += Foam::exp(-(sqr(coordX - 0.25*(xmax + xmin)) + sqr(coordY - front)) / sqr(xmax - xmin) * 200)/10;
    phase += Foam::exp(-(sqr(coordX - 0.50*(xmax + xmin)) + sqr(coordY - front)) / sqr(xmax - xmin) * 200)/10;
    phase += Foam::exp(-(sqr(coordX - 0.75*(xmax + xmin)) + sqr(coordY - front)) / sqr(xmax - xmin) * 200)/10;

    phase.write();
    assert(Foam::max(phase).value() <= 1 && Foam::min(phase).value() >= 0);

    T = alloy.liquidus() - undercooling + tempGradient * (coordY - ymin/3 - 2*ymax/3);
    T.write();

    Info<< "Dimensionless parameters:" << endl
        << " -- minimal undercooling = "
        << alloy.undercooling(Foam::min(T)).value() << endl
        << " -- maximum undercooling = "
        << alloy.undercooling(Foam::max(T)).value() << endl
        << " -- interface width / capillary length = "
        << (interfaceWidth / alloy.capillaryLength()).value() << endl
        << " -- mesh step / interface width = "
        << (1./ interfaceWidth / Foam::max(mesh.surfaceInterpolation::deltaCoeffs())).value() << endl;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        /** Temperature calculation */

        T = alloy.liquidus() - undercooling
            + tempGradient * (coordY - ymin/3 - 2*ymax/3)
            - coolingRate * runTime;

        while (pimple.loop()) {
            /** Phase field calculation */

            volVectorField gradPhase = fvc::grad(phase);
            volVectorField normal = (gradPhase + smallVector) / mag(gradPhase + smallVector);
            volScalarField theta = 0 * phase;
            forAll(theta, cellI) {
                theta[cellI] = Foam::atan2(normal[cellI].x(), normal[cellI].y());
            }
            volScalarField a_s = 1 + epsilon4 * Foam::cos(4 * theta);

            volVectorField A = 4 * epsilon4 * a_s * mag(gradPhase) * vector::one;
            forAll(A, cellI) {
                A[cellI].x() *= - Foam::sin(4 * theta[cellI]) * Foam::cos(theta[cellI]);
                A[cellI].y() *= + Foam::sin(4 * theta[cellI]) * Foam::sin(theta[cellI]);
            }

            tmp<volScalarField> f_prime = 2 * phase * (1 - phase) * (1 - 2 * phase);
            tmp<volScalarField> g_prime = 30 * sqr(phase) * sqr(1 - phase);
            fvScalarMatrix phaseEqn(
                tau * sqr(a_s) * fvm::ddt(phase)
                == sqr(interfaceWidth) * fvm::laplacian(sqr(a_s), phase)
                + sqr(interfaceWidth) * fvc::div(A) - f_prime
                + a1 * interfaceWidth / alloy.interfaceEnergy()
                    * g_prime * alloy.chemicalDrivingForce(phase, T)
            );
            phaseEqn.solve();

            /** Concentrations calculation */

            //normal = fvc::grad(phase) / mag(fvc::grad(phase));
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
