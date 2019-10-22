Info<< "Reading thermal properties\n" << endl;

IOdictionary thermalProperties
(
    IOobject
    (
        "thermalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    )
);

dimensionedScalar T0(thermalProperties.lookup("T0"));
dimensionedScalar k_sol(thermalProperties.lookup("k_sol"));
dimensionedScalar dk_sol(thermalProperties.lookup("dk_sol"));
dimensionedScalar k_liq(thermalProperties.lookup("k_liq"));
dimensionedScalar dk_liq(thermalProperties.lookup("dk_liq"));
dimensionedScalar Cp_sol(thermalProperties.lookup("Cp_sol"));
dimensionedScalar dCp_sol(thermalProperties.lookup("dCp_sol"));
dimensionedScalar Cp_liq(thermalProperties.lookup("Cp_liq"));
dimensionedScalar dCp_liq(thermalProperties.lookup("dCp_liq"));
dimensionedScalar T_solidus(thermalProperties.lookup("T_solidus"));
dimensionedScalar T_liquidus(thermalProperties.lookup("T_liquidus"));
dimensionedScalar T_boiling(thermalProperties.lookup("T_boiling"));
dimensionedScalar enthalpyFusion(thermalProperties.lookup("enthalpyFusion"));
dimensionedScalar enthalpyEvaporation(thermalProperties.lookup("enthalpyEvaporation"));
