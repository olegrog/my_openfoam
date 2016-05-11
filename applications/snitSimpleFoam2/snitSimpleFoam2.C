/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
-------------------------------------------------------------------------------
Application
    snitSimpleFoam2

Description
    Solves fluid-dynamic-type equations derived from the asymptotic Sone-Bobylev
    analysis of Boltzmann equation for hard-sphere gas with additional boundary
    conditions.

    For details look at Sone, Aoki, Takata, Sugimoto, Bobylev in Phys. Fluids 8
    (1996).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "fixedFluxPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        curvature.correctBoundaryConditions();

        /** temperature equation */
        fvScalarMatrix TEqn(
            fvm::laplacian(0.5 * gamma2 * sqrt(T), T) == fvc::div(U)
        );
        TEqn.solve();

        /** SIMPLE algorithm for solving pressure-velocity equations */

        /** velocity predictor */
        tmp<fvVectorMatrix> UEqn(
              fvm::div(phi, U)
            - fvm::laplacian(0.5 * gamma1 * sqrt(T), U) ==
              0.5 * gamma1 * fvc::div(sqrt(T) * dev2(::T(fvc::grad(U))))
            + gamma7 / T * (sqr(fvc::grad(T)) & (U / gamma2 / sqrt(T)))
        );
        UEqn().relax();
        solve(UEqn() ==
            fvc::reconstruct((
                - fvc::snGrad(p2)
                - 0.25 * gamma7 * fvc::interpolate(magSqr(fvc::grad(T)) / T) * fvc::snGrad(T)
            ) * mesh.magSf())
        );

        /** pressure corrector */
        volScalarField rAU(1./UEqn().A());
        surfaceScalarField rAUbyT("rhorAUf", fvc::interpolate(rAU / T));
        volVectorField HbyA("HbyA", U);
        HbyA = rAU * UEqn().H();
        UEqn.clear();

        surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(HbyA / T) & mesh.Sf());
        adjustPhi(phiHbyA, U, p2);
        surfaceScalarField phif(
            - 0.25 * gamma7 * rAUbyT
            * fvc::interpolate(magSqr(fvc::grad(T)) / T)
            * fvc::snGrad(T) * mesh.magSf()
        );
        phiHbyA += phif;

        // Update the fixedFluxPressure BCs to ensure flux consistency
        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
            p2.boundaryField(),
            phiHbyA.boundaryField()
            / (mesh.magSf().boundaryField() * rAUbyT.boundaryField())
        );

        while (simple.correctNonOrthogonal()) {
            fvScalarMatrix pEqn(
                fvm::laplacian(rAUbyT, p2) == fvc::div(phiHbyA)
            );
            pEqn.setReference(pRefCell, getRefCellValue(p2, pRefCell));
            pEqn.solve();
            if (simple.finalNonOrthogonalIter()) {
                // Calculate the conservative fluxes
                phi = phiHbyA - pEqn.flux();
                p2.relax();
                /** velocity corrector */
                U = HbyA + rAU * fvc::reconstruct((phif - pEqn.flux()) / rAUbyT);
                U.correctBoundaryConditions();
            }
        }
        #include "continuityErrs.H"

        // Pressure is defined up to a constant factor,
        // we adjust it to maintain the initial domainIntegrate
        p2 += (initialPressure - fvc::domainIntegrate(p2)) / totalVolume;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

