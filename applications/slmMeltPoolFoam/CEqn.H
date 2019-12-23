    volScalarField phase = 1 - liquidFraction;
    volScalarField h = pAlloy->partition(phase);
    forAllIter(PtrDictionary<alloyComponent>, pAlloy->components(), iter) {
        alloyComponent& C = iter();
        surfaceScalarField phi = fvc::snGrad(phase) * mesh.magSf()
            * fvc::interpolate(C.diffusion(phase) / sqr(h)) * pAlloy->partitionPrime();
        fvScalarMatrix CEqn(
            fvm::ddt(C) == fvm::laplacian(C.diffusion(phase) / h, C) - fvm::div(phi, C)
            - fvc::laplacian(C.diffusion(phase), C.equilibrium(phase, T) / h)
        );
        CEqn.solve();
    }
