/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
-------------------------------------------------------------------------------
Application
    snitSimpleFoam

Description
    Solves fluid-dynamic-type equations derived from the asymptotic Sone-Bobylev
    analysis of Boltzmann equation for hard-sphere gas.

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

        /** temperature equation */
        fvScalarMatrix TEqn(
            fvm::laplacian(0.5 * gamma2 * sqrt(T0), T0) == fvc::div(phi, T0)
        );
        TEqn.solve();

        /** SIMPLE algorithm for solving pressure-velocity equations */

        /** velocity predictor */
        tmp<fvVectorMatrix> UEqn(
              fvm::div(phi, U1)
            - fvm::laplacian(0.5 * gamma1 * sqrt(T0), U1) ==
              0.5 * gamma1 * fvc::div(sqrt(T0) * dev2(::T(fvc::grad(U1))))
        );
        UEqn().relax();
        solve(UEqn() ==
            fvc::reconstruct((
                - fvc::snGrad(p2)
                + gamma7 * fvc::snGrad(T0) * fvc::interpolate(
                    (U1 / gamma2 / sqrt(T0) - 0.25*fvc::grad(T0)) & fvc::grad(T0)/T0
                )
            ) * mesh.magSf())
        );

        /** pressure corrector */
        volScalarField rAU(1./UEqn().A());
        surfaceScalarField rAUbyT("rhorAUf", fvc::interpolate(rAU / T0));
        volVectorField HbyA("HbyA", U1);
        HbyA = rAU * UEqn().H();
        UEqn.clear();

        surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(HbyA / T0) & mesh.Sf());
        adjustPhi(phiHbyA, U1, p2);
        surfaceScalarField phif(
            gamma7 * rAUbyT * fvc::snGrad(T0) * mesh.magSf() * fvc::interpolate(
                (U1 / gamma2 / sqrt(T0) - 0.25*fvc::grad(T0)) & fvc::grad(T0)/T0
            )
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
                U1 = HbyA + rAU * fvc::reconstruct((phif - pEqn.flux()) / rAUbyT);
                U1.correctBoundaryConditions();
            }
        }
        #include "continuityErrs.H"

        // Pressure is defined up to a constant factor,
        // we adjust it to maintain the initial domainIntegrate
        p2 += (initialPressure - fvc::domainIntegrate(p2)) / totalVolume;
        p0 = totalVolume / fvc::domainIntegrate(1/T0);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

