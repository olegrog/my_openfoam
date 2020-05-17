    const dimensionedScalar zeroT(dimDensity/dimTemperature);

    const volVectorField normal(fvc::reconstruct(mixture.nHatf()));
    const volTensorField I_nn(tensor::one - sqr(normal));
    const volScalarField darcyDamping = mushyCoeff * alpha1
        * sqr(1 - liquidFraction.inMetal()) / (sqr(liquidFraction.inMetal()) + SMALL);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
      - fvm::Sp(darcyDamping, U)
      + dsigma_dT * (gradT & I_nn) * 2*rho / (rho1 + rho2)
            * mag(fvc::reconstruct(fvc::snGrad(alpha2) * mesh.magSf()))
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
                      - (3-evaporationCoeff)/4 * fvc::interpolate(vaporPressure) * fvc::snGrad(alpha2)
                    )*2*fvc::interpolate(rho)/(rho1 + rho2)
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );
    }
