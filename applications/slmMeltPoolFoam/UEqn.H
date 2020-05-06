    MRF.correctBoundaryVelocity(U);

    const dimensionedScalar zeroT(dimDensity/dimTemperature);

    const volVectorField normal(fvc::reconstruct(mixture.nHatf()));
    const volTensorField I_nn(tensor::one - normal*normal);
    const volScalarField darcyDamping = mushyCoeff * alpha1
        * sqr(1 - liquidFraction.inMetal()) / pow(liquidFraction.inMetal() + SMALL, 2);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
      - fvm::Sp(darcyDamping, U)
      + dsigma_dT * (fvc::grad(T) & I_nn)
            * magModifiedGradAlpha(T, liquidFraction, alpha2, rho1, rho1, zeroT, zeroT, rho2)
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
