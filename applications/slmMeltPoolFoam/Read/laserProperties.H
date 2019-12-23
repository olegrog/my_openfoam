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

dimensionedScalar laserPower("laserPower", laserProperties);
dimensionedScalar absorptivity("absorptivity", laserProperties);
dimensionedVector laserVelocity("laserVelocity", laserProperties);
dimensionedScalar laserRadius("laserRadius", laserProperties);
dimensionedVector coordStart("coordStart", laserProperties);
dimensionedScalar timeStop("timeStop", laserProperties);
