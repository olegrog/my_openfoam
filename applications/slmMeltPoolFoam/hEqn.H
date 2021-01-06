radiativeCooling = emissivity*sigma*(pow4(T) - pow4(ambientTemperature));
evaporativeCooling =
    evaporationCoeff*mixture.vapourPressure(ambientPressure)
   *Hvapour*Foam::sqrt(molarMass/2/pi/R/T);
const volScalarField evaporativeCoolingPrime =
    evaporativeCooling*molarMass*Hvapour/R/Cp/sqr(T)*sign(T - Tboiling);

const volScalarField Cp1 = thermo.Cp(T, liquidFraction, geometricUniformField<scalar>(0));
const volScalarField Cp2 = thermo.Cp(T, liquidFraction, geometricUniformField<scalar>(1));

absorptivity =
    2*rho*Cp/(rho1*Cp1 + rho2*Cp2)
   *max(dimensionedScalar(inv(dimLength)), gradAlpha1 & beamDirection)
   *(1 - (1 - absorptivity0)*exp(min(Zero, mesh.C().component(vector::Z)/2/laserRadius)));

pp(T);

fvScalarMatrix hEqn
(
    fvm::ddt(rho, h)
  + fvm::div(rhoPhi, h)
 ==
    fvm::laplacian(k*mixture.TPrimeEnthalpy(), h)
  + fvc::laplacian(k*mixture.TPrimeMetalFraction(), alpha1)
  + absorptivity*laserHeatSource
  + 2*rho*Cp/(rho1*Cp1 + rho2*Cp2)*mag(gradAlpha1)
   *(
      - radiativeCooling
      - evaporativeCooling
      - fvm::Sp(evaporativeCoolingPrime, h) + evaporativeCoolingPrime*h
    )
  //+ fvc::div(rho*U & turbulence->devReff())
);
hEqn.solve();
pp(h);

mixture.correctThermo();
