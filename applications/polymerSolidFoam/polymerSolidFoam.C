/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
      solidDisplacementFoam | Copyright (C) 2011-2016 OpenFOAM Foundation
           polymerSolidFoam | Copyright (C) 2020 Oleg Rogozin
-------------------------------------------------------------------------------
License
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    polymerSolidFoam

Description
    Segregated finite-volume solver of elastic small-strain deformation of a
    solid body with chemical stresses.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating displacement field\n" << endl;

    symmTensor Ixy = symmTensor::I; Ixy.zz() = 0;

    while (simple.loop())
    {
        Info<< "Iteration: " << runTime.value() << nl << endl;

        fvVectorMatrix DEqn
        (
            fvm::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
          + divSigmaExp
        );
        DEqn.solve();

        volTensorField gradD(fvc::grad(D));
        sigma = mu*twoSymm(gradD) + lambda*I*tr(gradD) - threeK*epsilonChemicalMax
            *(polymerization*I + interlayerCohesion*(nLayer - 1)*Ixy);
        divSigmaExp = fvc::div(sigma - (2*mu + lambda)*gradD, "div(sigma)");

        if (runTime.writeTime())
        {
            epsilon = symm(gradD);
            sigmaEq = sqrt(1.5*magSqr(dev(sigma)));
            Info<< "Max sigmaEq = " << max(sigmaEq).value() << endl;
            runTime.write();
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
