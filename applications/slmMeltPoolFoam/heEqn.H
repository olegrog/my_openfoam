{
    // mag(fvc::grad(alpha2) & beamDirection) can be used instead of mag(fvc::grad(alpha2))
    const dimensionedScalar totalEnthalpy = fvc::domainIntegrate(rho*he);

    heatCapacity = threePhaseParameter(T, meltIndicator, alpha2, Cp_sol, Cp_liq, dCp_sol, dCp_liq, Cp_gas);
    thermalConductivity = threePhaseParameter(T, meltIndicator, alpha2, k_sol, k_liq, dk_sol, dk_liq, k_gas);
    diffusivityf = fvc::interpolate(thermalConductivity/heatCapacity);
    laserHeatSource = (runTime < timeStop)
        * absorptivity * laserPower * surfaceGaussian(mesh.C(), laserCoordinate, laserRadius);
    radiativeCooling = emissivity * constant::physicoChemical::sigma
        * (pow4(T) - pow4(T_ambient));
    evaporativeCooling = evaporationCoeff * vaporPressure * enthalpyBoiling
        * Foam::sqrt(molarMass/2/constant::mathematical::pi/constant::physicoChemical::R/T);

    fvScalarMatrix heEqn
    (
        fvm::ddt(rho, he)
      + fvm::div(rhoPhi, he)
      - fvm::laplacian(thermalConductivity / heatCapacity, he)
      + enthalpyFusion
            * fvm::laplacian(alpha1 * thermalConductivity / heatCapacity * dMeltIndicator, he)
     ==
        (
            laserHeatSource
          - radiativeCooling
          - evaporativeCooling
        ) * magModifiedGradAlpha(T, meltIndicator, alpha2, rho1*Cp_sol, rho1*Cp_liq,
                rho1*dCp_sol, rho1*dCp_liq, rho2*Cp_gas)
      + fvc::div(rho * U & turbulence->devReff())
    );
    heEqn.solve();

    pp(he);

    Info<< "Real energy input = " << (fvc::domainIntegrate(rho*he) - totalEnthalpy).value()
        << ", laser input = " << fvc::domainIntegrate(laserHeatSource * runTime.deltaT()).value()
        << ", theoretical value = " << (absorptivity * laserPower * runTime.deltaT()).value()
        << endl;

    // Update other variables
    calcMeltIndicator(meltIndicator, he, he_melting, enthalpyFusion);
    calcMeltIndicator(dMeltIndicator, he, he_melting, enthalpyFusion, true);
    wasMelted = max(wasMelted, alpha1*meltIndicator);
    calcTemperature(T, he, meltIndicator, alpha2,
        Cp_sol, dCp_sol, Cp_liq, dCp_liq, Cp_gas, T_melting, enthalpyFusion);
    rhok = rho*(1.0 - beta*(T - T_ambient));
    vaporPressure = ambientPressure * exp(molarMass * enthalpyBoiling
        / constant::physicoChemical::R * (1/T_boiling - 1/T));

    pp(meltIndicator);
    pp(T);
}
