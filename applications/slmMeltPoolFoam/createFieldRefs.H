// References

volScalarField& h = mixture.h();
const volScalarField& Cp = mixture.Cp();
const volScalarField& k = mixture.k();
const volScalarField& liquidFraction = mixture.liquidFraction();

const dimensionedScalar& molarMass = thermo.metalM();
const dimensionedScalar& Hvapour = thermo.Hvapour();
const dimensionedScalar& Tboiling = thermo.Tboiling();

// Nondimensional constants

scalar emissivity(problemProperties.get<scalar>("emissivity"));
scalar evaporationCoeff(problemProperties.get<scalar>("evaporationCoeff"));
