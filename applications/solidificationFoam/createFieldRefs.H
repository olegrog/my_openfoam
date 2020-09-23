// Numerical constants
const auto a1(phaseFieldProperties.lookupOrDefault<scalar>("a1", Foam::sqrt(2.)/6));
const auto a2(phaseFieldProperties.lookupOrDefault<scalar>("a2", 2.35));
const auto a3(phaseFieldProperties.lookupOrDefault<scalar>("a3", 1));
const auto a4(phaseFieldProperties.lookupOrDefault<scalar>("a4", 10));
const auto aAT(phaseFieldProperties.lookupOrDefault<scalar>("aAT", 1/Foam::sqrt(2.)));
const auto aIF(phaseFieldProperties.lookupOrDefault<scalar>("aIF", 2.757));

// Conditions for stable and accurate simulations
const auto maxWbyD(phaseFieldProperties.get<scalar>("maxWbyD"));
const auto maxHbyW(phaseFieldProperties.get<scalar>("maxHbyW"));
const auto minLbyD(phaseFieldProperties.get<scalar>("minLbyD"));
const auto minSigmaW(phaseFieldProperties.get<scalar>("minSigmaW"));

// Model parameters
const auto nAnisotropyModes(alloy.get<label>("nAnisotropyModes"));
const auto capillaryAnisotropy(alloy.get<scalar>("capillaryAnisotropy"));
const auto kineticAnisotropy(alloy.get<scalar>("kineticAnisotropy"));

// Initial conditions
const auto frontPositionRel(problemProperties.get<scalar>("frontPosition"));
const auto interfaceNarrowing(problemProperties.get<scalar>("interfaceNarrowing"));
const auto seedNarrowing(problemProperties.get<scalar>("seedNarrowing"));
const auto nSeeds(problemProperties.get<label>("nSeeds"));
const auto addRandomSeeds(problemProperties.get<bool>("addRandomSeeds"));

