//TODO(optimization): calculate anisotropy and antitrapping terms only near the interface
//TODO(optimization): replace cos(n*atan) with polynomial

if (capillaryAnisotropy > 0 || kineticAnisotropy > 0)
{
    const volVectorField gradPhase = fvc::grad(phase);
    const volVectorField nHat = calcNormal(gradPhase);
    forAll(theta, cellI)
    {
        theta[cellI] = Foam::atan2(nHat[cellI].x(), nHat[cellI].y()) - theta0[cellI];
    }
    const volScalarField aS = 1 + capillaryAnisotropy*Foam::cos(nAnisotropyModes*theta);
    const volScalarField aTau = 1 + kineticAnisotropy*Foam::cos(nAnisotropyModes*theta);


    volVectorField A = nAnisotropyModes*capillaryAnisotropy*aS*mag(gradPhase)*vector::one;
    forAll(A, cellI)
    {
        A[cellI].x() *= - Foam::sin(nAnisotropyModes*theta[cellI])*Foam::cos(theta[cellI]);
        A[cellI].y() *= + Foam::sin(nAnisotropyModes*theta[cellI])*Foam::sin(theta[cellI]);
    }

    fvScalarMatrix phaseEqn
    (
        tau*sqr(aTau)*fvm::ddt(phase)
     ==
        sqr(interfaceWidth)*fvm::laplacian(sqr(aS), phase)
      + sqr(interfaceWidth)*fvc::div(A)
      - fPrime(phase)
      - fPrime2(phase)*(fvm::Sp(pTraits<scalar>::one, phase) - phase)
      + a1*interfaceWidth/alloy.interfaceEnergy()*gPrime(phase)
            *alloy.chemicalDrivingForce(phase, T)
    );
    phaseEqn.solve();
}
else
{
    fvScalarMatrix phaseEqn
    (
        tau*fvm::ddt(phase)
     ==
        sqr(interfaceWidth)*fvm::laplacian(phase)
      - fPrime(phase)
      - fPrime2(phase)*(fvm::Sp(pTraits<scalar>::one, phase) - phase)
      + a1*interfaceWidth/alloy.interfaceEnergy()*gPrime(phase)
            *alloy.chemicalDrivingForce(phase, T)
    );
    phaseEqn.solve();
}

Info<< "Solid fraction: min = " << gMin(phase)
    << " max = 1 + " << gMax(phase) - 1 << endl;
