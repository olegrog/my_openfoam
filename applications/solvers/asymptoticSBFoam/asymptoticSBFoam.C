/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
-------------------------------------------------------------------------------
Application
    asymptoticSBFoam

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
            fvm::laplacian(0.5 * gamma2 * sqrt(T), T) == fvm::div(phi, T)
        );
        TEqn.solve();

        /** SIMPLE algorithm for solving pressure-velocity equations */

        /** velocity predictor */
        tmp<fvVectorMatrix> UEqn(
              fvm::div(phi, U)
            - fvm::laplacian(0.5 * gamma1 * sqrt(T), U) ==
              0.5 * gamma1 * fvc::div(sqrt(T) * dev2(::T(fvc::grad(U))))
            + gamma7 / T * (sqr(fvc::grad(T)) & (U / gamma2 / sqrt(T) - 0.25 * fvc::grad(T)))
        );
        UEqn().relax();
        solve(UEqn() == -0.5 * fvc::grad(p));

        /** pressure corrector */
        volScalarField rAU(1./UEqn().A());
        volVectorField HbyA("HbyA", U);
        HbyA = rAU * UEqn().H();
        UEqn.clear();

        surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(HbyA / T) & mesh.Sf());
        adjustPhi(phiHbyA, U, p);

        // Update the fixedFluxPressure BCs to ensure flux consistency
        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
            p.boundaryField(),
            phiHbyA.boundaryField() * T.boundaryField()
            / (mesh.magSf().boundaryField() * rAU.boundaryField())
        );

        while (simple.correctNonOrthogonal()) {
            fvScalarMatrix pEqn(
                fvm::laplacian(0.5 * rAU / T, p) == fvc::div(phiHbyA)
            );
            pEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));
            pEqn.solve();
            if (simple.finalNonOrthogonalIter()) {
                phi = phiHbyA - pEqn.flux();
            }
        }
        #include "continuityErrs.H"
        p.relax();

        /** velocity corrector */
        U = HbyA - 0.5 * rAU * fvc::grad(p);
        U.correctBoundaryConditions();

        runTime.write();
    }

    Info<< "End\n" << endl;
    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

