{
    // phase + grain
    phase = Foam::atan(pow3((frontPosition - coord.component(vector::Y))/initialWidth))/pi + .5;
    const dimensionedScalar radius = width/nSeeds/seedNarrowing;
    for (int i = 0; i < nSeeds; i++)
    {
        dimensionedVector position = center;
        position.replace(vector::X, xmin + (i+.5)/nSeeds*width);
        position.replace(vector::Y, frontPosition);
        phase = max(phase, generateSeed(coord, position, radius));
    }
    phase.clip(0, 1);

    const volScalarField grainNumber
    (
        sign((coord.component(vector::X) - center.component(vector::X))/width)
    );
    addGrain(grain, phase*grainNumber, 1, nGrains);

    if (addRandomSeeds) {
        dimensionedVector position = center;
        position.replace(vector::Y, ymin/3 + 2*ymax/3);
        volScalarField seed = generateSeed(coord, position, radius/2);
        phase += seed;
        addGrain(grain, seed, 0, nGrains);
    }

    // number of grain
    calcNGrain(nGrain, grain, nGrains);

    // temperature
    T = alloy.liquidus() - undercooling
        + tempGradient*(coord.component(vector::Y) - frontPosition);

    // concentrations
    forAllIters(alloy.components(), iter)
    {
        alloyComponent& C = iter();
        const volScalarField temp = alloy.solidus()*phase + alloy.liquidus()*(1 - phase);
        C == C.equilibrium(phase, temp);
    }

}
