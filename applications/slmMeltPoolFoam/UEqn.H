    MRF.correctBoundaryVelocity(U);

    volVectorField normal("normal", fvc::reconstruct(mixture.nHatf()));
    volTensorField I_nn("I_nn", tensor::one - normal*normal);
    const dimensionedScalar zeroT(dimDensity/dimTemperature);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
      - fvm::Sp(permeability * sqr(1 - meltIndicator) / pow(meltIndicator + SMALL, 3) * alpha1, U)
      + dsigma_dT * (fvc::grad(T) & I_nn)
            * magModifiedGradAlpha(T, meltIndicator, alpha2, rho1, rho1, zeroT, zeroT, rho2)
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
                  - (3-evaporationCoeff)/4 * fvc::interpolate(vaporPressure) * fvc::snGrad(alpha2)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
