{
    radiativeCooling = emissivity*sigma*(pow4(T) - pow4(ambientTemperature));
    vaporPressure = ambientPressure*exp(molarMass*Hvapour/R*(1/Tboiling - 1/T));
    evaporativeCooling = evaporationCoeff*vaporPressure*Hvapour*Foam::sqrt(molarMass/2/pi/R/T);
    const volScalarField evaporativeCoolingPrime
    (
        evaporativeCooling*molarMass*Hvapour/R/Cp/sqr(T)*sign(T - Tboiling)
    );

    const volScalarField Cp1 = thermo.Cp(T, liquidFraction(), geometricUniformField<scalar>(0));
    const volScalarField Cp2 = thermo.Cp(T, liquidFraction(), geometricUniformField<scalar>(1));

    pp(liquidFraction());
    pp(T);

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h)
      + fvm::div(rhoPhi, h)
     ==
        fvm::laplacian(k/Cp*(1 - Hfusion*liquidFraction.dEnthalpy()), h)
      + fvc::laplacian(k/Cp*hPrimeGasFraction, alpha1)
      + 2*rho*Cp/(rho1*Cp1 + rho2*Cp2)*mag(fvc::grad(alpha1))
       *(
            laserHeatSource
          - radiativeCooling
          - evaporativeCooling
          - evaporativeCoolingPrime*fvm::SuSp(sign(T - Tboiling), h)
          + evaporativeCoolingPrime*h
        )
      + fvc::div(rho*U & turbulence->devReff())
    );
    hEqn.solve();
    pp(h);

    mixture.correct();
}
