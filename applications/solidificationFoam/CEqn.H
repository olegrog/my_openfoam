const volVectorField normal = calcNormal(fvc::grad(phase));
const volScalarField h = alloy.partition(phase);

forAllIter(PtrDictionary<alloyComponent>, alloy.components(), iter)
{
    alloyComponent& C = iter();

    const surfaceScalarField phi2 = fvc::snGrad(phase)*mesh.magSf()
        *fvc::interpolate(C.diffusion(phase)/sqr(h))*alloy.partitionPrime();

    fvScalarMatrix CEqn
    (
        fvm::ddt(C)
      + fvm::div(phi, C)
     ==
        fvm::laplacian(C.diffusion(phase)/h, C)
      - fvm::div(phi2, C)
      - fvc::laplacian(C.diffusion(phase), C.equilibrium(phase, T)/h)
      + fvc::div(interfaceWidth/Foam::sqrt(2.)*C.deltaA()*normal*fvc::ddt(phase))
    );
    CEqn.solve(mesh.solverDict("concentration"));
}

