if (correctPhi)
{
    rAU.ref() = 1.0/UEqn.A();
}
else
{
    rAU = 1.0/UEqn.A();
}

surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
  + fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf)
);

if (p_rgh.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p_rgh);
    fvc::makeAbsolute(phiHbyA, U);
}

surfaceScalarField phig
(
    (
        (
            mixture.surfaceTensionForce()
          + recoilPressure
        )*2*fvc::interpolate(rho/(rho1 + rho2), "interpolate(rho)")*rAUf
      - ghf*fvc::snGrad(rhok)*rAUf
    )*mesh.magSf()
);

phiHbyA += phig;

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, U, phiHbyA, rAUf);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        fvm::laplacian(rAUf, p_rgh) + mixture.divPhi() == fvc::div(phiHbyA)
    );

    p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

    p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - p_rghEqn.flux();

        p_rgh.relax();

        U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
        U.correctBoundaryConditions();
    }
}

const volScalarField& divU = mixture.divPhi();
#include "continuityErrs.H"

// Correct Uf if the mesh is moving
fvc::correctUf(Uf, U, phi);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

if (MRF->moving())
{
    phi -= MRF->phiRel();
}

p == p_rgh + rhok*gh;

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rhok*gh;
}

if (!correctPhi)
{
    rAU.clear();
}
