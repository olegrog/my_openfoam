// Update alpha1
#include "alphaSuSp.H"
advector.advect(Sp, fvc::div(phi)());

#include "rhofs.H"
rhoPhi = (rho1f - rho2f)*advector.alphaPhi() + rho2f*phi;

alpha2 = 1.0 - alpha1;

Info<< "Phase-1 volume fraction = "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") - 1 = " << max(alpha1).value() - 1
    << endl;
