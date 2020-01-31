{
    // mag(fvc::grad(alpha2) & beamDirection) can be used instead of mag(fvc::grad(alpha2))
    // vector beamDirection(0, 0, -1);
    dimensionedVector laserCoordinate("laserCoordinate", coordStart + laserVelocity * mesh.time());
    dimensionedScalar totalEnthalpy = fvc::domainIntegrate(rho*he);

    heatCapacity = threePhaseParameter(T, liquidFraction, alpha2, Cp_sol, Cp_liq, dCp_sol, dCp_liq, Cp_gas);
    thermalConductivity = threePhaseParameter(T, liquidFraction, alpha2, k_sol, k_liq, dk_sol, dk_liq, k_gas);
    laserHeatSource = (runTime < timeStop) * mag(fvc::grad(alpha2))
        * absorptivity * laserPower * surfaceGaussian(mesh.C(), laserCoordinate, laserRadius);

    fvScalarMatrix heEqn
    (
        fvm::ddt(rho, he)
        + fvm::div(rhoPhi, he)
        - fvm::laplacian(thermalConductivity / heatCapacity, he)
        + enthalpyFusion * fvc::laplacian(thermalConductivity / heatCapacity, liquidFraction)
        - fvc::div(rho * U & turbulence->devReff())
        - laserHeatSource
    );
    heEqn.solve();

    Info<< "Real energy input = " << (fvc::domainIntegrate(rho*he) - totalEnthalpy).value()
        << ", laser input = " << fvc::domainIntegrate(laserHeatSource * runTime.deltaT()).value()
        << ", theoretical value = " << (absorptivity * laserPower * runTime.deltaT()).value()
        << endl;

    // Update other variables
    calcLiquidFraction(liquidFraction, he, he_sol, he_liq);
    calcLiquidFractionDer(dLiquidFraction, he, he_sol, he_liq);
    wasMelted = max(wasMelted, liquidFraction);
    calcTemperature(T, he, liquidFraction,
        Cp_sol, Cp_liq, dCp_sol, dCp_liq, T_solidus, T_liquidus, enthalpyFusion);
}
