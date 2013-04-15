/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        volVectorField DT = fvc::grad(T);
        volTensorField DU = fvc::grad(U);
        p.storePrevIter();

        /** temperature equation */
        fvScalarMatrix TEqn(
            fvm::laplacian((g2/2)*sqrt(T), T) == tr(DU)
        );
        TEqn.solve();
        DT = fvc::grad(T);

        /** using SIMPLE algorithm for solving pressure-velocity equations */
        
        /** velocity predictor */
        fvVectorMatrix UEqn(
              fvm::div(phi, U)
            - fvm::laplacian(g1*sqrt(T)/2, U)
            + fvc::laplacian(g1*sqrt(T)/2, U)
            - g1 * fvc::div(sqrt(T)*dev(symm(DU)))
            - g7/T * (sqr(DT) & (U/g2/sqrt(T) - DT/4))
        );
        UEqn.relax();
        solve(UEqn == -fvc::grad(p));
        p.boundaryField().updateCoeffs();

        /** pressure corrector */
        volScalarField rAU = 1./UEqn.A();
        U = rAU*UEqn.H();
        phi = fvc::interpolate(U/T) & mesh.Sf();
        adjustPhi(phi, U, p);
        fvScalarMatrix pEqn(
            fvm::laplacian(rAU/T, p) == fvc::div(phi)
        );
        pEqn.setReference(0, 0);
        pEqn.solve();
        phi -= pEqn.flux();
        p.relax();

        /** velocity corrector */
        U -= rAU*fvc::grad(p);
        U.correctBoundaryConditions();

        runTime.write();
    }

    Info<< "End\n" << endl;
    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

