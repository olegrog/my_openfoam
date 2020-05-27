const volScalarField liquidFraction("liquidFraction", 1 - phase + SMALL);
const surfaceScalarField liquidFractionf(fvc::interpolate(liquidFraction));
const volScalarField& nu = laminarTransport.nu();

fvVectorMatrix UEqn
(
    fvm::ddt(U)
  + fvm::div(phi/liquidFractionf, U)
 ==
    fvm::laplacian(nu, U)
  - fvm::Sp(nu*hStar/sqr(interfaceWidth)*phase*(1-phase), U)
);

UEqn.relax();

if (pimple.momentumPredictor())
{
    solve(UEqn == -liquidFraction*fvc::grad(p));
}
