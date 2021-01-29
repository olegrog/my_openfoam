const volScalarField Cp1 = thermo.Cp(T, liquidFraction, geometricUniformField<scalar>(0));
const volScalarField Cp2 = thermo.Cp(T, liquidFraction, geometricUniformField<scalar>(1));

pp(T);

fvScalarMatrix hEqn
(
    fvm::ddt(rho, h)
 ==
    fvm::laplacian(k*mixture.TPrimeEnthalpy(), h)
  + fvc::laplacian(k*mixture.TPrimeMetalFraction(), alpha1)
  + 2*rho*Cp/(rho1*Cp1 + rho2*Cp2)
   *(
        absorptivity*laserHeatSource
    )
);
hEqn.solve();
pp(h);

mixture.correctThermo();
