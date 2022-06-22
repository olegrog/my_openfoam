/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
      solidDisplacementFoam | Copyright (C) 2011-2016 OpenFOAM Foundation
           polymerSolidFoam | Copyright (C) 2020 Oleg Rogozin
    polymerLayeredSolidFoam | Copyright (C) 2022 Oleg Rogozin
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
    polymerLayeredSolidFoam

Description
    Segregated finite-volume solver of elastic small-strain deformation of a
    solid body with chemical stresses and layering.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    const label nLayers = round(gMax(nLayer));
    runTime.setEndTime(2*nLayers*runTime.deltaT());

    while (runTime.loop())
    {
        activeCells = pos(runTime.timeIndex() + 2 - 2*nLayer);
        E = activeCells*Efull + (1 - activeCells)*E0;

        mu = E/(2.0*(1.0 + nu));
        lambda = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));
        threeK = E/(1.0 - 2.0*nu);

        // Set free BC at all the pathes before the final stage
        if (runTime.timeIndex() == 2*nLayers)
        {
            auto& DBf = D.boundaryFieldRef();
            forAll(DBf, patchi)
            {
                if (DBf[patchi].type() != "zeroTraction")
                {
                    Info<< "Set zeroTraction BC at patch " << DBf[patchi].patch().name()
                        << " instead of " << DBf[patchi].type() << endl;

                    DBf.set
                    (
                        patchi,
                        fvPatchField<vector>::New
                        (
                            "zeroTraction",
                            mesh.boundary()[patchi],
                            DBf[patchi].internalField()
                        )
                    );
                }
                else
                {
                    refCast<fixedGradientFvPatchVectorField>(DBf[patchi]).gradient() = Zero;
                }
            }
            D = 0*D;
            sigma = 2*mu*epsilonRes + lambda*tr(epsilonRes)*I
                - threeK*epsilonChemicalMax*polymerization*I;
            divSigmaExp = fvc::div(sigma, "div(sigma)");
        }

        Info<< "Iteration: " << runTime.value() << endl;
        Info<< "Calculating displacement field" << endl;

        label iCorr = 0;
        scalar initialResidual = 0;

        do
        {
            Info<< nl << "Correction: " << iCorr << endl;
            fvVectorMatrix DEqn
            (
                fvm::laplacian(2*mu + lambda, D, "laplacian(DD,D)") + divSigmaExp
            );
            initialResidual = DEqn.solve().max().initialResidual();

            gradD = fvc::grad(D);
            sigma = mu*(twoSymm(gradD) + 2*epsilonRes) + lambda*(tr(gradD) + tr(epsilonRes))*I
                - threeK*epsilonChemicalMax*polymerization*activeCells*I;
            divSigmaExp = fvc::div(sigma - (2*mu + lambda)*gradD, "div(sigma)");

        } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);

        epsilon = symm(gradD);
        sigmaEq = sqrt(1.5*magSqr(dev(sigma)));
        Info<< "Max sigmaEq = " << gMax(sigmaEq) << endl;

        if (runTime.timeIndex() % 2 == 1)
        {
            epsilonRes += epsilon*(1 - activeCells);
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
