if (adjustTimeStep)
{
    scalar maxDeltaTFact =
        min
        (
            min(maxCo/(CoNum + SMALL), maxAlphaCo/(alphaCoNum + SMALL)),
            maxDeltaTemp/(maxDeltaTempNum + SMALL)
        );

    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;

    scalar minDeltaT(runTime.controlDict().getOrDefault("minDeltaT", 1e-20));

    if (runTime.deltaTValue() < minDeltaT)
    {
        runTime.writeNow();
        returnReduceOr(true);   // MPI barrier
        FatalError
            << "The time step becomes too low!"
            << exit(FatalError);
    }
}

