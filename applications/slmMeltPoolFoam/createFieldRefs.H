// Define references

volScalarField& h = mixture.h();
volScalarField& alpha1 = mixture.alpha1();
volScalarField& alpha2 = mixture.alpha2();
const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
const volScalarField& Cp = mixture.Cp();
const volScalarField& k = mixture.k();
const volScalarField& liquidFraction = mixture.liquidFraction();
const volScalarField& vapourFraction = mixture.vapourFraction();

// Define nondimensional constants

scalar emissivity(problemProperties.get<scalar>("emissivity"));
scalar evaporationCoeff(problemProperties.get<scalar>("evaporationCoeff"));
