// Update temperature-dependent fields used in the momentum equation
rhok = rho*(1.0 - beta*(T - ambientTemperature));

surfaceScalarField recoilPressure
(
    (3 - evaporationCoeff)/4
   *fvc::interpolate(mixture.vapourPressure(ambientPressure), "interpolate(vapourPressure)")
   *fvc::snGrad(alpha1)
);

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  - fvm::Sp(contErr, U)
  + turbulence->divDevRhoReff(rho, U)
  + fvm::Sp(mixture.solidPhaseDamping(), U)
 ==
    mixture.marangoniForce()*2*rho/(rho1 + rho2)
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
                  + recoilPressure
                )*2*fvc::interpolate(rho/(rho1 + rho2), "interpolate(rho)")
              - ghf*fvc::snGrad(rhok)
              - fvc::snGrad(p_rgh)
            )*mesh.magSf()
        )
    );
}
