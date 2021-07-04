const volVectorField nHat = calcNormal(fvc::grad(phase));
const volScalarField k = alloy.partition(phase, T);

for (alloyComponent& C : alloy.components())
{
    const volScalarField D = C.diffusion(phase, k);
    const volScalarField Ceq = C.equilibrium(phase, T);

    fvScalarMatrix CEqn
    (
        fvm::ddt(k, C.byK())
      + fvm::div(phi, C.byK())
     ==
        fvm::laplacian(D, C.byK())
      - fvc::laplacian(D, Ceq/k)
      + fvc::div(aAT*interfaceWidth*C.delta(T)*nHat*fvc::ddt(phase))
    );
    CEqn.solve(mesh.solverDict("concentration"));

    C.volScalarField::operator=(k*C.byK());

    if (gMin(C) < SMALL)
    {
        Warning << "min(" << C.name() << ") = " << gMin(C) << endl;
    }
}

