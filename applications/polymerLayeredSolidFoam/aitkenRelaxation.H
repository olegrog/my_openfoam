if (iCorr == 0)
{
    aitkenInitialRes = gMax(mag(D.primitiveField()));
}

// aitken acceleration
aitkenDelta.storePrevIter();

// update delta
aitkenDelta = (D - D.prevIter()) / aitkenInitialRes;

// update relaxation factor
if (iCorr > 0)
{
    const vectorField b = aitkenDelta.internalField() - aitkenDelta.prevIter().internalField();

    scalar sumMagB = gSum(magSqr(b));
    if(sumMagB < SMALL)
    {
        Warning << "Aitken under-relaxation: denominator less then SMALL" << endl;
        sumMagB += SMALL;
    }

    aitkenTheta = -aitkenTheta*gSum(aitkenDelta.prevIter().internalField() & b)/sumMagB;
    Info<< "aitkenTheta = " << aitkenTheta << endl;
}

// correction to the latest U
D += aitkenTheta*aitkenDelta*aitkenInitialRes;
