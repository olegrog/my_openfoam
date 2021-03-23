radiativeCooling = emissivity*sigma*(pow4(T) - pow4(ambientTemperature));
evaporativeCooling =
    evaporationCoeff*mixture.vapourPressure(ambientPressure)*Hvap*Foam::sqrt(molarMass/2/pi/R/T);

const surfaceScalarField kByCp = mixture.kf()/mixture.Cpf();
const volScalarField redistribution =
    mixture.interfacialHeatSourceRedistribution(rho, rho1, rho2);
const volScalarField evaporativeCoolingPrime = evaporativeCooling*molarMass*Hvap/R/Cp/sqr(T);
// const surfaceScalarField liquidFractionPrimeEnthalpyf =
//     fvc::interpolate(mixture.liquidFractionPrimeEnthalpy());

energySource =
(
  - fvc::laplacian(Hfus*kByCp, liquidFraction)
  - fvc::laplacian(kByCp*mixture.HsPrimeAlphaGf(), alpha2)
  + redistribution*
    (
        laserHeatSource.value()
      - mag(gradAlpha1)*(radiativeCooling + evaporativeCooling)
    )
)/rho;

pp(T);

fvScalarMatrix hEqn
(
    fvm::ddt(rho, h)
  + fvm::div(rhoPhi, h)
 ==
    fvm::laplacian(kByCp, h)
//  - fvm::laplacian(Hfus*kByCp*liquidFractionPrimeEnthalpyf, h)
//  + fvc::laplacian(Hfus*kByCp*liquidFractionPrimeEnthalpyf, h)
  + rho*energySource
  + redistribution*mag(gradAlpha1)
   *(evaporativeCoolingPrime*h - fvm::Sp(evaporativeCoolingPrime, h))
  //+ fvc::div(rho*U & turbulence->devReff())
);
hEqn.solve();
pp(h);

if (gMin(h) < 0)
{
    runTime.writeAndEnd();
}

mixture.correctThermo();
