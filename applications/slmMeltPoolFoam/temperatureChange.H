scalar maxDeltaTemp
(
    runTime.controlDict().get<scalar>("maxDeltaTemp")
);

scalar maxDeltaTempNum = gMax(mag(energySource/Cp)())*runTime.deltaTValue();
scalar meanDeltaTempNum = gSum((energySource/Cp)())/mesh.nCells()*runTime.deltaTValue();

Info<< "Temperature change mean: " << meanDeltaTempNum
    << " max: " << maxDeltaTempNum << endl;

