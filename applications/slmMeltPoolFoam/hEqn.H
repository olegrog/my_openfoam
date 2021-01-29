radiativeCooling = emissivity*sigma*(pow4(T) - pow4(ambientTemperature));
evaporativeCooling =
    evaporationCoeff*mixture.vapourPressure(ambientPressure)
   *Hvapour*Foam::sqrt(molarMass/2/pi/R/T);

const volScalarField evaporativeCoolingPrime =
    evaporativeCooling*molarMass*Hvapour/R/Cp/sqr(T)*sign(T - Tboiling);
const volScalarField rdf = mixture.interfacialHeatSourceRedistribution(rho, rho1, rho2);

absorptivity =
   rdf*max(dimensionedScalar(inv(dimLength)), gradAlpha1 & beamDirection)
   *(1 - (1 - absorptivity0)*exp(min(Zero, mesh.C().component(vector::Z)/2/laserRadius)));

pp(T);

fvScalarMatrix hEqn
(
    fvm::ddt(rho, h)
  + fvm::div(rhoPhi, h)
 ==
    fvm::laplacian(mixture.kf()*mixture.TPrimeEnthalpyf(), h)
  + fvc::laplacian(mixture.kf()*mixture.TPrimeMetalFractionf(), alpha1)
  + absorptivity*laserHeatSource
  + rdf*mag(gradAlpha1)
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
