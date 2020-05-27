const dimensionedScalar zero(0), one(1);

if (epsilon4 > 0)
{
    volVectorField gradPhase = fvc::grad(phase);
    volVectorField normal = calcNormal(gradPhase);
    forAll(theta, cellI)
    {
        theta[cellI] = Foam::atan2(normal[cellI].x(), normal[cellI].y()) - theta0[cellI];
    }
    volScalarField a_s = 1 + epsilon4*Foam::cos(4*theta);

    volVectorField A = 4*epsilon4*a_s*mag(gradPhase)*vector::one;
    forAll(A, cellI)
    {
        A[cellI].x() *= - Foam::sin(4*theta[cellI])*Foam::cos(theta[cellI]);
        A[cellI].y() *= + Foam::sin(4*theta[cellI])*Foam::sin(theta[cellI]);
    }

    fvScalarMatrix phaseEqn
    (
        tau*sqr(a_s)*fvm::ddt(phase)
     ==
        sqr(interfaceWidth)*fvm::laplacian(sqr(a_s), phase)
      + sqr(interfaceWidth)*fvc::div(A)
      - fPrime(phase)
      - fPrime2(phase)*(fvm::Sp(one, phase) - phase)
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
      - fPrime2(phase)*(fvm::Sp(one, phase) - phase)
      + a1*interfaceWidth/alloy.interfaceEnergy()*gPrime(phase)
            *alloy.chemicalDrivingForce(phase, T)
    );
    phaseEqn.solve();
}

