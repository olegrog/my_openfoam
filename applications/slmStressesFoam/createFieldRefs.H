volScalarField& h = mixture.h();
const volScalarField& T = mixture.T();
const volScalarField& liquidFraction = mixture.liquidFraction();

const gasMetalThermo& thermo = mixture.thermo();
const dimensionedScalar& Hfus = thermo.Hfusion();

