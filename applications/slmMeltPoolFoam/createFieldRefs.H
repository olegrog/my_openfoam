// References

volScalarField& h = mixture.h();
const volScalarField& Cp = mixture.Cp();
const volScalarField& liquidFraction = mixture.liquidFraction();

const gasMetalThermo& thermo = mixture.thermo();
const dimensionedScalar& molarMass = thermo.metalM();
const dimensionedScalar& Hfus = thermo.Hfusion();
const dimensionedScalar& Hvap = thermo.Hvapour();
const dimensionedScalar& Tboiling = thermo.Tboiling();

// Nondimensional constants

scalar emissivity(problemProperties.get<scalar>("emissivity"));
scalar evaporationCoeff(problemProperties.get<scalar>("evaporationCoeff"));
