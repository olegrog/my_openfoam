pp(T);

const surfaceScalarField kByCp = mixture.kf()/mixture.Cpf();

// const surfaceScalarField liquidFractionPrimeEnthalpyf =
//     fvc::interpolate(mixture.liquidFractionPrimeEnthalpy());

fvScalarMatrix hEqn
(
    fvm::ddt(rho, h)
 ==
    fvm::laplacian(kByCp, h)
//  - fvm::laplacian(Hfus*kByCp*liquidFractionPrimeEnthalpyf, h)
//  + fvc::laplacian(Hfus*kByCp*liquidFractionPrimeEnthalpyf, h)
  - fvc::laplacian(Hfus*kByCp, liquidFraction)
  - fvc::laplacian(kByCp*mixture.HsPrimeAlphaGf(), alpha2)
  + (*laserHeatSource)()
);
hEqn.solve();
pp(h);

mixture.correctThermo();
