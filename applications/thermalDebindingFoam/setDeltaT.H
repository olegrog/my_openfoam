if (adjustTimeStep)
{
    // 1. Calculate the Courant number

    surfaceScalarField phi = fvc::interpolate(D1 + D2)*fvc::snGrad(rho)*mesh.magSf();
    dimensionedScalar rhoP = polymer.rho();
    scalarField sumPhi((fvc::surfaceSum(mag(phi))/(rho + rhoP*SMALL))().primitiveField());
    scalar CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();
    scalar meanCoNum = 0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();

    Info<< "Courant Number mean: " << meanCoNum  << " max: " << CoNum << endl;

    // 2. Calculate the maximum change of polymer

    scalar maxDeltaY(runTime.controlDict().get<scalar>("maxDeltaY"));
    scalar totalRHS = 0;
    forAllIters(polymer.components(), iter)
    {
        const auto& y = iter();
        volScalarField RHS = y*y.degradationRate(T);
        totalRHS += gMax(RHS);
    }
    scalar maxDeltaYNum = runTime.deltaTValue()*totalRHS;

    // 3. Set the time step

    scalar maxDeltaTFact =
        min
        (
            maxCo/(CoNum + SMALL),
            maxDeltaY/(maxDeltaYNum + SMALL)
        );
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.5);
    runTime.setDeltaT(min(deltaTFact*runTime.deltaTValue(), maxDeltaT));

    Info<< "deltaT = " << runTime.deltaTValue() << endl;
}

