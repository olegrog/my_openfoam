wasMelted = Foam::max(wasMelted, liquidFraction);

mixture.correctPassiveFields();

if (writeAllFields)
{
    heatConvection = fvc::div(rhoPhi, h);

    // For debug: check that heatConduction is almost equal to heatConduction2
    heatConduction = fvc::laplacian(mixture.kappaf(), T);

    const surfaceScalarField kappaByCp = mixture.kappaf()/mixture.Cpf();

    heatConduction2 =
        fvc::laplacian(kappaByCp, h)
      - fvc::laplacian(Hfus*kappaByCp, liquidFraction)
      + fvc::laplacian(kappaByCp*mixture.HsPrimeAlphaGf(), alpha1);
}


