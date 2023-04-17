CorrectPhi
(
    U,
    phi,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    mixture.divPhi(),
    pimple
);

#include "continuityErrs.H"
