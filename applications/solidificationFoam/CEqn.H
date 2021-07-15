const volVectorField nHat = calcNormal(fvc::grad(phase));
const volScalarField k = alloy.partition(phase, T);

for (alloyComponent& C : alloy.components())
{
    const volScalarField DbyK = C.diffusion(phase, k)/k;
    const volScalarField Ceq = C.equilibrium(phase, T);

    fvScalarMatrix CEqn
    (
        fvm::ddt(C)
      + fvm::div(phi/fvc::interpolate(k), C)
     ==
        fvm::laplacian(DbyK, C)
      - fvc::laplacian(DbyK, Ceq)
      - fvc::laplacian(DbyK*(C - Ceq)/k*alloy.partitionPrimePhase(T), phase)
      - fvc::laplacian(DbyK*(C - Ceq)/k*alloy.partitionPrimeT(phase, T), T)
      + fvc::div(aAT*interfaceWidth*C.delta(T)*nHat*fvc::ddt(phase))
    );
    CEqn.solve(mesh.solverDict("concentration"));

    if (gMin(C) < SMALL)
    {
        Warning << "min(" << C.name() << ") = " << gMin(C) << endl;
    }
}

