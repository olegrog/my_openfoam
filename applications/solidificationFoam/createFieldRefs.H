scalar a1(phaseFieldProperties.lookupOrDefault<scalar>("a1", Foam::sqrt(2.)/6));
scalar a2(phaseFieldProperties.lookupOrDefault<scalar>("a2", 2.35));
scalar a3(phaseFieldProperties.lookupOrDefault<scalar>("a3", 1));
scalar a4(phaseFieldProperties.lookupOrDefault<scalar>("a4", 10));
scalar epsilon4(phaseFieldProperties.lookupOrDefault<scalar>("epsilon4", 0));
scalar hStar(phaseFieldProperties.lookupOrDefault<scalar>("hStar", 2.757));

scalar frontPositionRel(problemProperties.lookupOrDefault<scalar>("frontPosition", 0));
scalar interfaceNarrowing(problemProperties.lookupOrDefault<scalar>("interfaceNarrowing", 1));
scalar seedNarrowing(problemProperties.lookupOrDefault<scalar>("seedNarrowing", 1));
label nSeeds(problemProperties.lookupOrDefault<label>("nSeeds", 0));
bool addRandomSeeds(problemProperties.lookupOrDefault<bool>("addRandomSeeds", false));

