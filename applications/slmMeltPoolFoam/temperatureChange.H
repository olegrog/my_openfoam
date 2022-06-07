scalar maxDeltaTemp
(
    runTime.controlDict().get<scalar>("maxDeltaTemp")
);

const scalarField TSource = energySource.primitiveField()/Cp.primitiveField();
const scalar maxDeltaTempNum = gMax(mag(TSource))*runTime.deltaTValue();
const scalar meanDeltaTempNum = gSum(TSource)/mesh.nCells()*runTime.deltaTValue();

Info<< "Temperature change mean: " << meanDeltaTempNum
    << " max: " << maxDeltaTempNum << endl;

