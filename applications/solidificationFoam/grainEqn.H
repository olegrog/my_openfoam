// Calculate grain number and crystallographic orientation

volVectorField gNormal = calcNormal(grain);
volVectorField gTangent = rot & gNormal;

fvVectorMatrix grainEqn
(
    tau*fvm::ddt(grain)
  + a3*phase*(gNormal*(mag(grain) - 1) + gTangent*sin(2*pi*nGrain))
 ==
    a4*pow3(interfaceWidth)*fvm::laplacian(mag(fvc::grad(phase)), grain)
);
grainEqn.solve();

calcNGrain(nGrain, grain, nGrains);
forAll(theta0, cellI)
{
    theta0[cellI] =
        crystallographicAngles.lookupOrDefault(name(lround(nGrain[cellI])), 0)*pi/180;
}


