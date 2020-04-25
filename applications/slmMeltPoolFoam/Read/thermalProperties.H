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

dimensionedScalar T_ambient("T_ambient", thermalProperties);
dimensionedScalar T_melting("T_melting", thermalProperties);
dimensionedScalar T_boiling("T_boiling", thermalProperties);

dimensionedScalar enthalpyFusion("enthalpyFusion", thermalProperties);
dimensionedScalar enthalpyBoiling("enthalpyBoiling", thermalProperties);

auto he_S = [&](dimensionedScalar temp) {
    return Cp_sol*temp + dCp_sol*sqr(temp)/2;
};
auto he_L = [&](dimensionedScalar temp) {
    return Cp_liq*temp + dCp_liq*sqr(temp)/2;
};
auto he_G = [&](dimensionedScalar temp) {
    return Cp_gas*temp;
};

dimensionedScalar he_ambient = he_S(T_ambient);
dimensionedScalar he_melting = he_S(T_melting) + enthalpyFusion/2;
dimensionedScalar he_preboiling = he_melting + enthalpyFusion/2 + he_L(T_boiling) - he_L(T_melting);
dimensionedScalar zero(0), one(1);

Info<< "Ambient enthalpy = " << he_ambient.value() << endl
    << "Melting enthalpy = " << he_melting.value() << endl
    << "Gas enthalpy at the melting temperature = " << he_G(T_melting).value() << endl
    << "Pre-boiling enthalpy = " << he_preboiling.value() << endl
    << "Solidus heat capacity = " << threePhaseParameter(
            T_melting, zero, zero, Cp_sol, Cp_liq, dCp_sol, dCp_liq, Cp_gas).value() << endl
    << "Liquidus heat capacity = " << threePhaseParameter(
            T_melting, one, zero, Cp_sol, Cp_liq, dCp_sol, dCp_liq, Cp_gas).value() << endl
    << "Solidus heat conductivity = " << threePhaseParameter(
            T_melting, zero, zero, k_sol, k_liq, dk_sol, dk_liq, k_gas).value() << endl
    << "Liquidus heat conductivity = " << threePhaseParameter(
            T_melting, one, zero, k_sol, k_liq, dk_sol, dk_liq, k_gas).value() << endl
    << endl;
