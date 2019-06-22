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

tmp<volScalarField> fPrime(const volScalarField& phase) {
    return 2 * phase * (1 - phase) * (1 - 2 * phase);
}

tmp<volScalarField> gPrime(const volScalarField& phase) {
    return 30 * sqr(phase) * sqr(1 - phase);
}

tmp<volScalarField> addSeed(
    const volScalarField& coordX,
    const volScalarField& coordY,
    const dimensionedScalar& centerX,
    const dimensionedScalar& centerY,
    const dimensionedScalar& radius)
{
    return 2*Foam::exp(-sqr((mag(coordX - centerX) + mag(coordY - centerY)) / radius));
}

tmp<volVectorField> calcNormal(const volVectorField& gradPhase) {
    dimensionedVector smallVector("small", dimless / dimLength, vector(0, 1e-20, 0));
    return (gradPhase + smallVector) / mag(gradPhase + smallVector);
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    pimpleControl pimple(mesh);

    #include "createFields.H"

    /** Initial conditions */

    const boundBox& bounds = mesh.bounds();
    const dimensionedScalar xmax("xmax", dimLength, bounds.max().x());
    const dimensionedScalar xmin("xmin", dimLength, bounds.min().x());
    const dimensionedScalar ymax("ymax", dimLength, bounds.max().y());
    const dimensionedScalar ymin("ymin", dimLength, bounds.min().y());
    const dimensionedScalar front = ymin + (ymax - ymin) * frontPosition;
    const dimensionedScalar initialWidth = interfaceWidth / interfaceNarrowing;

    phase = Foam::atan(pow3((front - coordY) / initialWidth)) / mathematicalConstant::pi + 0.5;
    const dimensionedScalar radius = (xmax - xmin) / nSeeds / seedNarrowing;
    for (int i = 0; i < nSeeds; i++) {
        phase += addSeed(coordX, coordY, xmin + (i+.5)/nSeeds*(xmax - xmin), front, radius);
    }

    //phase += addSeed(coordX, coordY, (xmax + xmin)/2, (ymax + ymin)/2, radius);

    phase = min(max(phase, scalar(0)), scalar(1));
    phase.write();

    T = alloy.liquidus() - undercooling + tempGradient * (coordY - ymin/3 - 2*ymax/3);
    T.write();

    forAllIter(PtrDictionary<alloyComponent>, alloy.components(), iter) {
        alloyComponent& C = iter();
        C == C.equilibrium(phase, T);
        C.write();
    }

    dimensionedScalar tip_velocity = coolingRate / tempGradient;

    Info<< "Dimensionless parameters:" << endl
        << " -- minimal undercooling = "
        << alloy.undercooling(Foam::min(T)).value() << endl
        << " -- maximum undercooling = "
        << alloy.undercooling(Foam::max(T)).value() << endl
        << " -- interface width / capillary length = "
        << (interfaceWidth / alloy.capillaryLength()).value() << endl
        << " -- mesh step / interface width = "
        << (1./ interfaceWidth / Foam::max(mesh.surfaceInterpolation::deltaCoeffs())).value() << endl
        << " -- interface stability parameter = "
        << (alloy.diffusionL() * alloy.capillaryLength() / tip_velocity / sqr(interfaceWidth)).value() << nl << endl;

    Info<< "Dimensioned parameters:" << endl
        << " -- tip velocity (m/s) = "
        << tip_velocity.value() << endl;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run()) {
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        /** Temperature calculation */

        T = alloy.liquidus() - undercooling
            + tempGradient * (coordY - ymin/3 - 2*ymax/3)
            - coolingRate * runTime;

        /** Adaptive time step */

#       include "createTimeControls.H" // read adjustTimeStep, maxCo, maxDeltaT from controlDict
        scalar minDeltaT =
            runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", 0);
        if (adjustTimeStep) {
            // Estimate the interface velocity via rhs
            volScalarField rhsPhase = (sqr(interfaceWidth) * fvc::laplacian(phase) - fPrime(phase)
                + a1 * interfaceWidth / alloy.interfaceEnergy()
                    * gPrime(phase) * alloy.chemicalDrivingForce(phase, T)) / tau;
            Info<< "Max interface change = "
                << max(mag(rhsPhase) * runTime.deltaT()).value() << endl;
            Co = mag(rhsPhase) * runTime.deltaT();
            // Estimate evolution of concentration via rhs
            forAllIter(PtrDictionary<alloyComponent>, alloy.components(), iter) {
                alloyComponent& C = iter();
                volScalarField h = alloy.partition(phase);
                volVectorField normal = calcNormal(fvc::grad(phase));
                volScalarField rhsC = (
                    fvc::laplacian(C.diffusion(phase), (C - C.equilibrium(phase, T)) / h)
                    + fvc::div(interfaceWidth / Foam::sqrt(2.) * C.deltaA() * normal * rhsPhase)
                ) / C.equilibrium(phase, T);
                Info<< "Max relative change of " << C.name() << " = "
                    << max(mag(rhsC)* runTime.deltaT()).value() << endl;
                Co = max(Co, mag(rhsC) * runTime.deltaT());
            }

            // Reset the timestep to maintain a constant maximum courant Number.
            // Reduction of time-step is immediate, but increase is damped to avoid
            // unstable oscillations.
            scalar maxDeltaTFact = maxCo/(max(Co).value() + SMALL);
            scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

            runTime.setDeltaT(min(
                deltaTFact*runTime.deltaT().value(),
                maxDeltaT
            ));
            if (deltaTFact*runTime.deltaT().value() < minDeltaT) {
                Info<< "Time step becomes too small!" << endl;
                runTime.writeAndEnd();
            }

            Info<< "Time = " << runTime.timeName()
                << " deltaT = " <<  runTime.deltaT().value()
                << " deltaTFact = " << deltaTFact
                << endl;
        }
        while (pimple.loop()) {
            /** Phase field calculation */

            volVectorField gradPhase = fvc::grad(phase);
            volVectorField normal = calcNormal(gradPhase);
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

            fvScalarMatrix phaseEqn(
                tau * sqr(a_s) * fvm::ddt(phase)
                == sqr(interfaceWidth) * fvm::laplacian(sqr(a_s), phase)
                + sqr(interfaceWidth) * fvc::div(A) - fPrime(phase)
                + a1 * interfaceWidth / alloy.interfaceEnergy()
                    * gPrime(phase) * alloy.chemicalDrivingForce(phase, T)
            );
            phaseEqn.solve();

            /** Concentrations calculation */

            // Recompute normal
            normal = calcNormal(fvc::grad(phase));

            volScalarField h = alloy.partition(phase);
            forAllIter(PtrDictionary<alloyComponent>, alloy.components(), iter) {
                alloyComponent& C = iter();

                surfaceScalarField phi = fvc::snGrad(phase) * mesh.magSf()
                    * fvc::interpolate(C.diffusion(phase) / sqr(h)) * alloy.partitionPrime();
                fvScalarMatrix CEqn(
                    fvm::ddt(C) == fvm::laplacian(C.diffusion(phase) / h, C) - fvm::div(phi, C)
                    - fvc::laplacian(C.diffusion(phase), C.equilibrium(phase, T) / h)
                    + fvc::div(interfaceWidth / Foam::sqrt(2.) * C.deltaA() * normal * fvc::ddt(phase))
                );
                CEqn.solve(mesh.solutionDict().solver("concentration"));
            }
        }

        runTime.write();
        if (!adjustTimeStep) {
            Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
        }
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
