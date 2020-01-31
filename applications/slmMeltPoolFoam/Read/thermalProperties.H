Info<< "Reading thermal properties\n" << endl;

IOdictionary thermalProperties
(
    IOobject
    (
        "thermalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar k_sol("k_sol", thermalProperties);
dimensionedScalar dk_sol("dk_sol", thermalProperties);
dimensionedScalar k_liq("k_liq", thermalProperties);
dimensionedScalar dk_liq("dk_liq", thermalProperties);
dimensionedScalar k_gas("k_gas", thermalProperties);

dimensionedScalar Cp_sol("Cp_sol", thermalProperties);
dimensionedScalar dCp_sol("dCp_sol", thermalProperties);
dimensionedScalar Cp_liq("Cp_liq", thermalProperties);
dimensionedScalar dCp_liq("dCp_liq", thermalProperties);
dimensionedScalar Cp_gas("Cp_gas", thermalProperties);

dimensionedScalar T_solidus("T_solidus", thermalProperties);
dimensionedScalar T_liquidus("T_liquidus", thermalProperties);
dimensionedScalar T_boiling("T_boiling", thermalProperties);

dimensionedScalar enthalpyFusion("enthalpyFusion", thermalProperties);
dimensionedScalar enthalpyBoiling("enthalpyBoiling", thermalProperties);

dimensionedScalar T_1 = (T_liquidus + T_solidus)/2;
dimensionedScalar dT_1 = (T_liquidus - T_solidus)/2;
dimensionedScalar he_sol = Cp_sol*T_solidus + dCp_sol*sqr(T_solidus)/2;
dimensionedScalar he_1 = Cp_sol*T_1 + dCp_sol*sqr(T_1)/2 + enthalpyFusion/2;
dimensionedScalar he_liq = he_1 + Cp_liq*dT_1 + dCp_liq*sqr(dT_1)/2 + enthalpyFusion/2;
dimensionedScalar zero(0), one(1);

Info<< "Solidus enthalpy = " << he_sol.value() << endl
    << "Liquidus enthalpy = " << he_liq.value() << endl
    << "Solidus heat capacity = " << threePhaseParameter(
            T_solidus, zero, zero, Cp_sol, Cp_liq, dCp_sol, dCp_liq, Cp_gas).value() << endl
    << "Liquidus heat capacity = " << threePhaseParameter(
            T_liquidus, one, zero, Cp_sol, Cp_liq, dCp_sol, dCp_liq, Cp_gas).value() << endl
    << "Solidus heat conductivity = " << threePhaseParameter(
            T_solidus, zero, zero, k_sol, k_liq, dk_sol, dk_liq, k_gas).value() << endl
    << "Liquidus heat conductivity = " << threePhaseParameter(
            T_liquidus, one, zero, k_sol, k_liq, dk_sol, dk_liq, k_gas).value() << endl;
