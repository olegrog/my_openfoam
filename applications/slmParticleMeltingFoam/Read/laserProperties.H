Info<< "Reading laser properties\n" << endl;

IOdictionary laserProperties
(
    IOobject
    (
        "laserProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar laserPower(laserProperties.lookup("laserPower"));
dimensionedScalar absorptivity(laserProperties.lookup("absorptivity"));
dimensionedVector laserVelocity(laserProperties.lookup("laserVelocity"));
dimensionedScalar laserRadius(laserProperties.lookup("laserRadius"));
dimensionedVector coordStart(laserProperties.lookup("coordStart"));
dimensionedScalar timeStop(laserProperties.lookup("timeStop"));
