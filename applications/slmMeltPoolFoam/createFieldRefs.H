// References

volScalarField& h = mixture.h();
const volScalarField& Cp = mixture.Cp();
const volScalarField& liquidFraction = mixture.liquidFraction();

const gasMetalThermo& thermo = mixture.thermo();
const dimensionedScalar& molarMass = thermo.metalM();
const dimensionedScalar& Hfus = thermo.Hfusion();
const dimensionedScalar& Hvap = thermo.Hvapour();

// Nondimensional constants

const scalar emissivity(problemProperties.get<scalar>("emissivity"));
const scalar evaporationCoeff(problemProperties.get<scalar>("evaporationCoeff"));

// Variables
scalar meshUpdateTime = 0;
