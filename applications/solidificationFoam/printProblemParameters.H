{
    const scalar minMeshStep = 1/gMax(mesh.surfaceInterpolation::deltaCoeffs());
    const dimensionedScalar Tmin("Tmin", dimTemperature, gMin(T));
    const dimensionedScalar Tmax("Tmax", dimTemperature, gMax(T));
    const dimensionedScalar PDAS = (xmax - xmin)/nSeeds;
    const dimensionedScalar diffusionLength = alloy.diffusionL()/pullingVelocity;

    const scalar sigmaW =
    (
        alloy.diffusionL()*alloy.capillaryLength()/pullingVelocity/sqr(interfaceWidth)
    ).value();
    const scalar WbyD = (interfaceWidth/alloy.capillaryLength()).value();
    const scalar HbyW = (minMeshStep/interfaceWidth).value();
    const scalar LbyD = ((xmax - xmin)/diffusionLength).value();

    Info<< nl << "Dimensionless parameters:" << endl
        << " -- minimal undercooling = " << alloy.undercooling(Tmin).value() << endl
        << " -- maximum undercooling = " << alloy.undercooling(Tmax).value() << endl
        << " -- interface width / capillary length = " << WbyD << " < " << maxWbyD << endl
        << " -- mesh step / interface width = " << HbyW << " < " << maxHbyW << endl
        << " -- computational width / diffusion length = " << LbyD << " > " << minLbyD << endl
        << " -- interface stability parameter = " << sigmaW << " > " << minSigmaW << endl
        << " -- pulling velocity = "
        << (pullingVelocity*alloy.capillaryLength()/alloy.diffusionL()).value() << endl
        << " -- Reynolds number = " << max(mag(U)/laminarTransport.nu()*PDAS).value() << endl
        << nl;

    Info<< "Dimensioned parameters:" << endl
        << " -- computational domain (m) = " << bounds.max() - bounds.min() << endl
        << " -- generated PDAS (m) = " << PDAS.value() << endl
        << " -- diffusion length (m) = " << diffusionLength.value() << endl
        << " -- interface width (m) = " << interfaceWidth.value() << endl
        << " -- mesh step (m) = " << minMeshStep << endl
        << " -- capillary length (m) = " << alloy.capillaryLength().value() << endl
        << " -- pulling velocity (m/s) = " << pullingVelocity.value() << endl
        << " -- characteristic velocity (m/s) = "
        << (alloy.diffusionL()/alloy.capillaryLength()).value() << endl
        << " -- relaxation time (s) = "  << tau.value() << endl
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
