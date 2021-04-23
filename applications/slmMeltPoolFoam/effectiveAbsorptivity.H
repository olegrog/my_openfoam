const dimensionedScalar effectiveLaserPower = fvc::domainIntegrate((*laserHeatSource)());

Info<< "Real energy input = " << (fvc::domainIntegrate(rho*h) - totalEnthalpy).value()
    << ", laser input = " << (effectiveLaserPower*runTime.deltaT()).value()
    << ", effective absorptivity = " << (effectiveLaserPower/laserHeatSource->power()).value()
    << endl;
