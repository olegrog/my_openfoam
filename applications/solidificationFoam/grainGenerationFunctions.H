tmp<volScalarField> generateSeed
(
    const volVectorField& coord,
    const dimensionedVector& center,
    const dimensionedScalar& radius
)
{
    return min
    (
        2*Foam::exp(-magSqr((coord - center)/radius)),
        scalar(1)
    );
}


void addGrain(volVectorField& grain, const volScalarField& phase, label nGrain, label nGrains)
{
    forAll(grain, cellI)
    {
        scalar argument = 2*pi*sign(phase[cellI])*nGrain/nGrains;
        scalar magnitude = fabs(phase[cellI]);
        grain[cellI].x() += magnitude*Foam::cos(argument);
        grain[cellI].y() += magnitude*Foam::sin(argument);
    }
}


void calcNGrain(volScalarField& nGrain, const volVectorField& grain, label nGrains)
{
    forAll(grain, cellI)
    {
        nGrain[cellI] = Foam::atan2(grain[cellI].y(), grain[cellI].x())/2/pi*nGrains;
    }
}


