{
    // mag(fvc::grad(alpha2) & beamDirection) can be used instead of mag(fvc::grad(alpha2))
    const dimensionedScalar totalEnthalpy = fvc::domainIntegrate(rho*he);

    // We have to update all alpha-dependent fields
    liquidFraction.update();
    calcTemperature(T, he, liquidFraction, alpha2,
        Cp_sol, dCp_sol, Cp_liq, dCp_liq, Cp_gas, T_melting, enthalpyFusion);
    calcEnthalpy(heDAlpha, T, liquidFraction, alpha2,
        Cp_sol, dCp_sol, Cp_liq, dCp_liq, Cp_gas, T_melting, enthalpyFusion, true);
    heatCapacity = threePhaseParameter(T, liquidFraction, alpha2,
        Cp_sol, Cp_liq, dCp_sol, dCp_liq, Cp_gas);
    thermalConductivity = threePhaseParameter(T, liquidFraction, alpha2,
        k_sol, k_liq, dk_sol, dk_liq, k_gas);
    radiativeCooling = emissivity * constant::physicoChemical::sigma
        * (pow4(T) - pow4(T_ambient));
    vaporPressure = ambientPressure * exp(molarMass * enthalpyBoiling
        / constant::physicoChemical::R * (1/T_boiling - 1/T));
    evaporativeCooling = evaporationCoeff * vaporPressure * enthalpyBoiling
        * Foam::sqrt(molarMass/2/constant::mathematical::pi/constant::physicoChemical::R/T);

    pp(liquidFraction);
    pp(T);

    fvScalarMatrix heEqn
    (
        fvm::ddt(rho, he)
      + fvm::div(rhoPhi, he)
      - fvm::laplacian(thermalConductivity / heatCapacity
            * (1 - enthalpyFusion * liquidFraction.dEnthalpy()), he)
      + fvc::laplacian(thermalConductivity / heatCapacity
            * (heDAlpha + enthalpyFusion * liquidFraction.dAlpha()), alpha2)
     ==
        (
            laserHeatSource
          - radiativeCooling
          - evaporativeCooling
        ) * magModifiedGradAlpha(T, liquidFraction, alpha2, rho1*Cp_sol, rho1*Cp_liq,
                rho1*dCp_sol, rho1*dCp_liq, rho2*Cp_gas)
      + fvc::div(rho * U & turbulence->devReff())
    );
    heEqn.solve();
    pp(he);

    Info<< "Real energy input = " << (fvc::domainIntegrate(rho*he) - totalEnthalpy).value()
        << ", laser input = " << fvc::domainIntegrate(laserHeatSource * runTime.deltaT()).value()
        << ", theoretical value = " << (absorptivity * laserPower * runTime.deltaT()).value()
        << endl;

    if (pimple.corrPISO() == pimple.nCorrPISO())
    {
        // Update alpha-dependent quantities that are used in the momentum equation
        liquidFraction.update();
        calcTemperature(T, he, liquidFraction, alpha2,
            Cp_sol, dCp_sol, Cp_liq, dCp_liq, Cp_gas, T_melting, enthalpyFusion);
        pp(liquidFraction);
        pp(T);
        rhok = rho*(1.0 - beta*(T - T_ambient));
        vaporPressure = ambientPressure * exp(molarMass * enthalpyBoiling
            / constant::physicoChemical::R * (1/T_boiling - 1/T));

        // For debug: check that heatConduction is almost equal to heatConduction2
        heatConvection = fvc::div(rhoPhi, he);
        heatConduction = fvc::laplacian(thermalConductivity, T);
        heatConduction2 = fvc::laplacian(thermalConductivity / heatCapacity
                * (1 - enthalpyFusion * liquidFraction.dEnthalpy()), he)
            - fvc::laplacian(thermalConductivity / heatCapacity
                * (heDAlpha + enthalpyFusion * liquidFraction.dAlpha()), alpha2);
    }
}
