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
    dimensionedScalar tipTemp("tipTemp", T.dimensions(), T[tipCell]);

    Info<< "Tip position(m) = " << (tipPosition).value()
        << " relative = " << ((tipPosition - ymin)/height).value() << endl
        << "Tip speed(m/s) = " << tipSpeed.value()
        << " dimless = "
        << (tipSpeed*alloy.capillaryLength(tipTemp)/alloy.diffusionL(tipTemp)).value() << endl
        << "Tip undercooling(K) = " << (alloy.liquidus() - tipTemp).value()
        << " relative = " << alloy.undercooling(tipTemp).value()
        << nl << endl;

    tipPositionPrev = tipPosition;
    tipCellPrev = tipCell;
    tipTimeInterval *= 0;
}

