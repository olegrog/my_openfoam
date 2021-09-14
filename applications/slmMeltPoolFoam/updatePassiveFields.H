wasMelted = Foam::max(wasMelted, liquidFraction);

mixture.correctPassiveFields();

if (writeAllFields)
{
    heatConvection = fvc::div(rhoPhi, h);

    // For debug: check that heatConduction is almost equal to heatConduction2
    heatConduction = fvc::laplacian(mixture.kf(), T);

    const surfaceScalarField kByCp = mixture.kf()/mixture.Cpf();

    heatConduction2 =
        fvc::laplacian(kByCp, h)
      - fvc::laplacian(Hfus*kByCp, liquidFraction)
      - fvc::laplacian(kByCp*mixture.HsPrimeAlphaGf(), alpha2);
}


