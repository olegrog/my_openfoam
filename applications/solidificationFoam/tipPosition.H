tipTimeInterval += runTime.deltaT();

// Find the current tip position and the corresponding cell
tipPosition = ymin;
forAll(phase, cellI)
{
    if (phase[cellI] > .5)
    {
        dimensionedScalar tipPositionNew =
            dimensionedScalar("y", dimLength, coord[cellI].y())
          + pi*(phase[cellI] - .5)*interfaceWidth;

        if (tipPositionNew > tipPosition)
        {
            tipPosition = tipPositionNew;
            tipCell = cellI;
        }
    }
}

if (tipCell != tipCellPrev)
{
    dimensionedScalar tipSpeed = (tipPosition - tipPositionPrev)/tipTimeInterval;
    Info<< "Tip position(m) = " << (tipPosition).value()
        << " relative = " << ((tipPosition - ymin)/height).value() << endl
        << "Tip speed(m/s) = " << tipSpeed.value()
        << " dimless = "
        << (tipSpeed*alloy.capillaryLength()/alloy.diffusionL()).value() << endl
        << "Tip undercooling(K) = " << alloy.liquidus().value() - T[tipCell]
        << " relative = "
        << alloy.undercooling(dimensionedScalar("T", dimTemperature, T[tipCell])).value()
        << nl << endl;

    tipPositionPrev = tipPosition;
    tipCellPrev = tipCell;
    tipTimeInterval *= 0;
}

