volScalarField rAU("rAU", liquidFraction/UlEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UlEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

adjustPhi(phiHbyA, U, p);

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0/max(1.0/rAU - UlEqn.H1()/liquidFraction, 0.1/rAU);
    phiHbyA += liquidFractionf*fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= liquidFraction*(rAU - rAtU())*fvc::grad(p);
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU());

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU()*liquidFraction, p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        // Calculate the conservative fluxes
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

Ul = HbyA - rAtU*liquidFraction*fvc::grad(p);
Ul.correctBoundaryConditions();
U = liquidFraction*Ul;

