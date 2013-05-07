/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.2.0                                 |
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
#include "fvIOoptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    simpleControl simple(mesh);

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        volVectorField DT = fvc::grad(T);
        volTensorField DU = fvc::grad(U);

        /** temperature equation */
        fvScalarMatrix TEqn(
            fvm::laplacian((g2/2)*sqrt(T), T) == tr(DU)
        );
        TEqn.solve();
        DT = fvc::grad(T);

        /** SIMPLE algorithm for solving pressure-velocity equations */
        
        /** velocity predictor */
        volScalarField invT(1/T);
        tmp<fvVectorMatrix> UEqn(
              fvm::div(phi, U)
            - fvm::laplacian(g1*sqrt(T)/2, U)
            + fvc::laplacian(g1*sqrt(T)/2, U)
            - g1 * fvc::div(sqrt(T)*dev(symm(DU)))
            - g7/T * (sqr(DT) & (U/g2/sqrt(T) - DT/4))
	    == fvOptions(invT, U)
        );
        UEqn().relax();
	fvOptions.constrain(UEqn());
        solve(UEqn() == -fvc::grad(p));
	fvOptions.correct(U);

        /** pressure corrector */
        volScalarField rAU(1./UEqn().A());
	volVectorField HbyA("HbyA", U);
        HbyA = rAU*UEqn().H();
	UEqn.clear();

	surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(HbyA/T) & mesh.Sf());
        adjustPhi(phiHbyA, U, p);
	fvOptions.relativeFlux(phiHbyA);

	while (simple.correctNonOrthogonal()) {
            fvScalarMatrix pEqn(
                fvm::laplacian(rAU/T, p) == fvc::div(phiHbyA)
            );
            pEqn.setReference(0, 0);
            pEqn.solve();
            if (simple.finalNonOrthogonalIter())
		phi = phiHbyA - pEqn.flux();
	}
	#include "continuityErrs.H"
        p.relax();
	
        /** velocity corrector */
        U = HbyA - rAU*fvc::grad(p);
        U.correctBoundaryConditions();
	fvOptions.correct(U);

        runTime.write();
    }

    Info<< "End\n" << endl;
    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

