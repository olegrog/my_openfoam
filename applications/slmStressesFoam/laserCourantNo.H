const scalar minMeshStep = 1/gMax(mesh.surfaceInterpolation::deltaCoeffs());
const scalar CoNum = mag(laserHeatSource->velocity()).value()*runTime.deltaTValue()/minMeshStep;

Info<< "Laser Courant number = " << CoNum << endl;
