    MRF.correctBoundaryVelocity(U);

    volVectorField normal("normal", fvc::reconstruct(mixture.nHatf()));
    dimensionedScalar evapor = enthalpyBoiling * molarMass / constant::physicoChemical::R;

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
      + fvm::Sp(permeability * sqr(1-liquidFraction) / pow(liquidFraction + SMALL, 3) * alpha1, U)
     ==
        fvOptions(rho, U)
      + mag(fvc::grad(alpha1)) * marangoni * (fvc::grad(T) - normal * (normal & fvc::grad(T)))
//      - anisimov * ambientPressure * evapor * fvc::snGrad(T) *
//        fvc::interpolate(exp(evapor * alpha1 * (1./T_boiling - 1./T)) / sqr(T))
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
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
