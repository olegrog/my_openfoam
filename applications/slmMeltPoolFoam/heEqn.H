{
    // We have to update all alpha-dependent fields
    liquidFraction.correct();
    calcTemperature(T, he, liquidFraction, alpha2,
        Cp_sol, dCp_sol, Cp_liq, dCp_liq, Cp_gas, T_melting, enthalpyFusion);
    calcEnthalpy(hePrimeAlpha, T, liquidFraction, alpha2,
        Cp_sol, dCp_sol, Cp_liq, dCp_liq, Cp_gas, T_melting, enthalpyFusion, true);
    heatCapacity = threePhaseParameter(T, liquidFraction.sharp(), alpha2,
        Cp_sol, Cp_liq, dCp_sol, dCp_liq, Cp_gas);
    thermalConductivity = threePhaseParameter(T, liquidFraction(), alpha2,
        k_sol, k_liq, dk_sol, dk_liq, k_gas);
    radiativeCooling = emissivity * constant::physicoChemical::sigma
        * (pow4(T) - pow4(T_ambient));
    vaporPressure = ambientPressure * exp(molarMass * enthalpyBoiling
        / constant::physicoChemical::R * (1/T_boiling - 1/T));
    evaporativeCooling = evaporationCoeff * vaporPressure * enthalpyBoiling
        * Foam::sqrt(molarMass/2/constant::mathematical::pi/constant::physicoChemical::R/T);

    volScalarField Cp1 = threePhaseParameter(T, liquidFraction.sharp(), zero,
        Cp_sol, Cp_liq, dCp_sol, dCp_liq, Cp_gas);
    volScalarField Cp2 = threePhaseParameter(T, liquidFraction.sharp(), one,
        Cp_sol, Cp_liq, dCp_sol, dCp_liq, Cp_gas);

    pp(liquidFraction());
    pp(T);

    fvScalarMatrix heEqn
    (
        fvm::ddt(rho, he)
      + fvm::div(rhoPhi, he)
      - fvm::laplacian(thermalConductivity / heatCapacity
            * (1 - enthalpyFusion * liquidFraction.dEnthalpy()), he)
      + fvc::laplacian(thermalConductivity / heatCapacity * hePrimeAlpha, alpha2)
     ==
        (
            laserHeatSource
          - radiativeCooling
          - evaporativeCooling
        ) * 2*rho*heatCapacity / (rho1*Cp1 + rho2*Cp2)
            * mag(fvc::reconstructMag(fvc::snGrad(alpha2) * mesh.magSf()))
      + fvc::div(rho * U & turbulence->devReff())
    );
    heEqn.solve();
    pp(he);

    if (pimple.corrPISO() == pimple.nCorrPISO())
    {
        // Update alpha-dependent quantities that are used in the momentum equation
        liquidFraction.correct();
        calcTemperature(T, he, liquidFraction, alpha2,
            Cp_sol, dCp_sol, Cp_liq, dCp_liq, Cp_gas, T_melting, enthalpyFusion);
        pp(liquidFraction());
        pp(T);
        rhok = rho*(1.0 - beta*(T - T_ambient));
        vaporPressure = ambientPressure * exp(molarMass * enthalpyBoiling
            / constant::physicoChemical::R * (1/T_boiling - 1/T));

        // For debug: check that heatConduction is almost equal to heatConduction2
        heatConvection = fvc::div(rhoPhi, he);
        heatConduction = fvc::laplacian(thermalConductivity, T);
        heatConduction2 = fvc::laplacian(thermalConductivity / heatCapacity
                * (1 - enthalpyFusion * liquidFraction.dEnthalpy()), he)
            - fvc::laplacian(thermalConductivity / heatCapacity * hePrimeAlpha, alpha2);
    }
}
