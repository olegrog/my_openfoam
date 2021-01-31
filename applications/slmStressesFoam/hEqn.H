pp(T);

const surfaceScalarField kByCp = mixture.kf()/mixture.Cpf();
const volScalarField redistribution =
    mixture.interfacialHeatSourceRedistribution(rho, rho1, rho2);

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
  + redistribution*absorptivity*laserHeatSource
);
hEqn.solve();
pp(h);

mixture.correctThermo();
