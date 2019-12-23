{
    volScalarField Cp = fourParameterModel(liquidFraction, T, Cp_sol, Cp_liq, dCp_sol, dCp_liq, T0);
    volScalarField k = fourParameterModel(liquidFraction, T, k_sol, k_liq, dk_sol, dk_liq, T0) * alpha1;
    diffusivity = k / Cp;

    dimensionedVector laserCoordinate("laserCoordinate", coordStart + laserVelocity * mesh.time());
    tmp<volScalarField> laserHeatSource = mag(fvc::grad(alpha1))
        * absorptivity * laserPower * gaussian(mesh.C(), laserCoordinate, laserRadius);
    tmp<volScalarField> fusionTerm = fvc::laplacian(diffusivity * enthalpyFusion, liquidFraction);

    fvScalarMatrix heEqn
    (
        fvm::ddt(rho, he)
        + fvm::div(rhoPhi, he)
        - fvm::laplacian(diffusivity, he)
        - fvc::div(rho * U & turbulence->devReff())
        - (runTime < timeStop) * laserHeatSource()
        + fusionTerm()
    );
    heEqn.solve();

    liquidFraction = tanhSmooth(he, (he_solidus + he_liquidus)/2, (he_liquidus - he_solidus)/2);
    wasMelted = max(wasMelted, liquidFraction);
    T = temperatureCalc(he, liquidFraction, Cp_sol, Cp_liq, dCp_sol, dCp_liq, T_solidus, T_liquidus, enthalpyFusion);
}


