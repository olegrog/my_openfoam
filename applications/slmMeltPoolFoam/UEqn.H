// gradT = fvc::grad(T) is a less smooth alternative
gradT = (1 - Hfusion*liquidFraction.dEnthalpy())/Cp*fvc::grad(h)
    + hPrimeGasFraction/Cp*fvc::grad(alpha1);

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  + turbulence->divDevRhoReff(rho, U)
 ==
  - fvm::Sp(mixture.solidPhaseDamping(), U)
  + mixture.marangoniForce(gradT)*2*rho/(rho1 + rho2)
);

UEqn.relax();

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                (
                    mixture.surfaceTensionForce()
                  - (3 - evaporationCoeff)/4*fvc::interpolate(vaporPressure)*fvc::snGrad(alpha2)
                )*2*fvc::interpolate(rho)/(rho1 + rho2)
              - ghf*fvc::snGrad(rhok)
              - fvc::snGrad(p_rgh)
            )*mesh.magSf()
        )
    );
}
