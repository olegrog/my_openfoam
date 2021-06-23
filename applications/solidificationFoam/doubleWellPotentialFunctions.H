auto fPrime = [](const volScalarField& phase)
{
    return 2*phase*(1 - phase)*(1 - 2*phase);
};

auto fPrime2 = [](const volScalarField& phase)
{
    return 3*sqr(1 - 2*phase) - 1;
};

auto gPrime = [](const volScalarField& phase)
{
    return 30*sqr(phase)*sqr(1 - phase);
};

