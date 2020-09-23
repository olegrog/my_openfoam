{
    const scalar minMeshStep = 1/gMax(mesh.surfaceInterpolation::deltaCoeffs());
    const dimensionedScalar Tmin("Tmin", dimTemperature, gMin(T));
    const dimensionedScalar Tmax("Tmax", dimTemperature, gMax(T));
    const dimensionedScalar seedSpacing = (xmax - xmin)/nSeeds;
    const dimensionedScalar diffusionLength = alloy.diffusionL()/pullingSpeed;

    // Steady-state estimations from Kundin et al. 2019
    const scalar anisotropyStrength = max(capillaryAnisotropy, kineticAnisotropy);
    const scalar tipSpeed = (pullingSpeed*alloy.capillaryLength()/alloy.diffusionL()).value();
    const scalar tipUndercooling = Foam::pow(pow4(pi)*tipSpeed/820/anisotropyStrength, 1./5.4);
    dimensionedScalar tipRadius
    (
        2*pi*Foam::sqrt(alloy.interfaceEnergy()/pullingSpeed/alloy.sumRestrictionFactors())
    );
    if (anisotropyStrength > VSMALL) tipRadius *= Foam::pow(50*anisotropyStrength, -0.225);
    const dimensionedScalar cellSpacing =
        Foam::sqrt(3*(1 - tipUndercooling)*alloy.deltaTemp()/tempGradient*tipRadius);

    const scalar sigmaW =
    (
        alloy.diffusionL()*alloy.capillaryLength()/pullingSpeed/sqr(interfaceWidth)
    ).value();
    const scalar WbyD = (interfaceWidth/alloy.capillaryLength()).value();
    const scalar HbyW = (minMeshStep/interfaceWidth).value();
    const scalar LbyD = ((xmax - xmin)/diffusionLength).value();

    // Check the dimensions
    seedSpacing + cellSpacing + diffusionLength + tipRadius + interfaceWidth
        + alloy.capillaryLength() + pullingSpeed*tau;

    Info<< nl << "Dimensioned parameters:" << endl
        << " -- computational domain (m) = " << bounds.max() - bounds.min() << endl
        << " -- estimated cell spacing (m) = " << cellSpacing.value() << endl
        << " -- seed spacing (m) = " << seedSpacing.value() << endl
        << " -- estimated tip radius (m) = " << tipRadius.value() << endl
        << " -- diffusion length (m) = " << diffusionLength.value() << endl
        << " -- interface width (m) = " << interfaceWidth.value() << endl
        << " -- mesh step (m) = " << minMeshStep << endl
        << " -- capillary length (m) = " << alloy.capillaryLength().value() << endl
        << " -- pulling speed (m/s) = " << pullingSpeed.value() << endl
        << " -- characteristic speed (m/s) = "
        << (alloy.diffusionL()/alloy.capillaryLength()).value() << endl
        << " -- relaxation time (s) = "  << tau.value() << endl
        << nl;

    Info<< "Dimensionless parameters:" << endl
        << " -- minimal undercooling = " << alloy.undercooling(Tmin).value() << endl
        << " -- maximum undercooling = " << alloy.undercooling(Tmax).value() << endl
        << " -- estimated tip undercooling = " << tipUndercooling << endl
        << " -- tip speed = " << tipSpeed << endl
        << " -- Reynolds number = " << max(mag(U)/laminarTransport.nu()*seedSpacing).value() << endl
        << nl;

    Info<< "Accuracy and stability parameters:" << endl
        << " -- interface width / capillary length = " << WbyD << " < " << maxWbyD << endl
        << " -- mesh step / interface width = " << HbyW << " < " << maxHbyW << endl
        << " -- computational width / diffusion length = " << LbyD << " > " << minLbyD << endl
        << " -- interface stability parameter = " << sigmaW << " > " << minSigmaW << endl
        << nl;

    if (WbyD > maxWbyD)
    {
        FatalError
            << "(interface width / capillary length) is too high:" << endl
            << WbyD << " > " << maxWbyD
            << abort(FatalError);
    }

    if (HbyW > maxHbyW)
    {
        FatalError
            << "(mesh step / interface width) is too low:" << endl
            << HbyW << " > " << maxHbyW
            << abort(FatalError);
    }

    if (LbyD < minLbyD)
    {
        FatalError
            << "(computational width / diffusion length) is too low:" << endl
            << LbyD << " < " << minLbyD
            << abort(FatalError);
    }

    if (sigmaW < minSigmaW)
    {
        FatalError
            << "Interface stability parameter is too low:" << endl
            << sigmaW << " < " << minSigmaW
            << abort(FatalError);
    }

}
