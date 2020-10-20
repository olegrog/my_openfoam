const volVectorField normal = calcNormal(fvc::grad(phase));
const volScalarField partition = alloy.partition(phase);

forAllIters(alloy.components(), iter)
{
    alloyComponent& C = iter();

    const surfaceScalarField phi2 = fvc::snGrad(phase)*mesh.magSf()
        *fvc::interpolate(C.diffusion(phase)/sqr(partition))*alloy.partitionPrime();

    fvScalarMatrix CEqn
    (
        fvm::ddt(C)
      + fvm::div(phi/fvc::interpolate(alloy.partition(phase)), C)
     ==
        fvm::laplacian(C.diffusion(phase)/partition, C)
      - fvm::div(phi2, C)
      - fvc::laplacian(C.diffusion(phase), C.equilibrium(phase, T)/partition)
      + fvc::div(aAT*interfaceWidth*C.deltaA()*normal*fvc::ddt(phase))
    );
    CEqn.solve(mesh.solverDict("concentration"));
}

