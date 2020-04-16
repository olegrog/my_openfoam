    MRF.correctBoundaryVelocity(U);

    volVectorField normal("normal", fvc::reconstruct(mixture.nHatf()));
    volTensorField I_nn("I_nn", tensor::one - normal*normal);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
      - fvm::Sp(permeability * sqr(1 - liquidFraction) / pow(liquidFraction + SMALL, 3) * alpha1, U)
      + dsigma_dT * mag(fvc::grad(alpha2)) * (fvc::grad(T) & I_nn)
      - (3 - evaporationCoeff)/4 * vaporPressure * fvc::grad(alpha2)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
