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

dimensionedScalar T0("T0", thermalProperties);
dimensionedScalar k_sol("k_sol", thermalProperties);
dimensionedScalar dk_sol("dk_sol", thermalProperties);
dimensionedScalar k_liq("k_liq", thermalProperties);
dimensionedScalar dk_liq("dk_liq", thermalProperties);
dimensionedScalar Cp_sol("Cp_sol", thermalProperties);
dimensionedScalar dCp_sol("dCp_sol", thermalProperties);
dimensionedScalar Cp_liq("Cp_liq", thermalProperties);
dimensionedScalar dCp_liq("dCp_liq", thermalProperties);
dimensionedScalar T_solidus("T_solidus", thermalProperties);
dimensionedScalar T_liquidus("T_liquidus", thermalProperties);
dimensionedScalar T_boiling("T_boiling", thermalProperties);
dimensionedScalar enthalpyFusion("enthalpyFusion", thermalProperties);
dimensionedScalar enthalpyBoiling("enthalpyBoiling", thermalProperties);
