const volVectorField nHat = calcNormal(fvc::grad(phase));
const volScalarField k = alloy.partition(phase, T);

for(alloyComponent& C : alloy.components())
{
    const volScalarField D = C.diffusion(phase, k);
    const surfaceScalarField phi2 = fvc::snGrad(k)*fvc::interpolate(D/sqr(k))*mesh.magSf();

    fvScalarMatrix CEqn
    (
        fvm::ddt(C)
      + fvm::div(phi/fvc::interpolate(k), C)
     ==
        fvm::laplacian(D/k, C)
      - fvm::div(phi2, C)
      - fvc::laplacian(D, C.equilibrium(phase, T)/k)
      + fvc::div(aAT*interfaceWidth*C.delta(T)*nHat*fvc::ddt(phase))
    );
    CEqn.solve(mesh.solverDict("concentration"));

    if (gMin(C) < SMALL)
    {
        Warning << "min(" << C.name() << ") = " << gMin(C) << endl;
    }
}

