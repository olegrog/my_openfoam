/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    wallShearStress1

Description
    Calculates and reports wall shear stress for all patches, for the
    specified times.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();
        #include "createFields.H"

        const surfaceVectorField shearStress1(
            twoSymm(dev(fvc::interpolate(
                -g1 * sqrt(T) * fvc::grad(U)
            ))) & mesh.Sf()
        );

        // the second derivative of T is asymmetric due to interpolation
        const surfaceVectorField shearStress2(
            symm(dev(fvc::interpolate(
                g3 * T * fvc::grad(fvc::grad(T))
            ))) & mesh.Sf()
        );

        const surfaceVectorField shearStress3(
            dev(fvc::interpolate(
                g7 * fvc::grad(T) * fvc::grad(T)
            )) & mesh.Sf()
        );

        const surfaceVectorField shearStress(
            shearStress1 + shearStress2 + shearStress3
        );

        Info<< "\nThe force acting on" << endl;
        forAll(shearStress.boundaryField(), patchi)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchi]))
            {
                Info<< mesh.boundary()[patchi].name()
                    << " "
                    << gSum(shearStress.boundaryField()[patchi])
                    << " consisted of "
                    << gSum(shearStress1.boundaryField()[patchi])
                    << gSum(shearStress2.boundaryField()[patchi])
                    << gSum(shearStress3.boundaryField()[patchi])
                    << endl;
            }
        }
        Info<< endl;

        dimensionedVector zeroShearStress(
            "wallShearStress",
            shearStress.dimensions(),
            vector::zero
        );

        volVectorField wallShearStress1
        (
            IOobject(
                "wallShearStress1",
                runTime.timeName(),
                mesh
            ),
            mesh,
            zeroShearStress
        );
        volVectorField wallShearStress2
        (
            IOobject(
                "wallShearStress2",
                runTime.timeName(),
                mesh
            ),
            mesh,
            zeroShearStress
        );
        volVectorField wallShearStress3
        (
            IOobject(
                "wallShearStress3",
                runTime.timeName(),
                mesh
            ),
            mesh,
            zeroShearStress
        );
        volVectorField wallShearStress
        (
            IOobject(
                "wallShearStress",
                runTime.timeName(),
                mesh
            ),
            mesh,
            zeroShearStress
        );

        forAll(wallShearStress1.boundaryField(), patchi)
        {
            wallShearStress1.boundaryField()[patchi] =
                - shearStress1.boundaryField()[patchi] / mesh.magSf().boundaryField()[patchi];
        }
        wallShearStress1.write();
        forAll(wallShearStress2.boundaryField(), patchi)
        {
            wallShearStress2.boundaryField()[patchi] =
                - shearStress2.boundaryField()[patchi] / mesh.magSf().boundaryField()[patchi];
        }
        wallShearStress2.write();
        forAll(wallShearStress3.boundaryField(), patchi)
        {
            wallShearStress3.boundaryField()[patchi] =
                - shearStress3.boundaryField()[patchi] / mesh.magSf().boundaryField()[patchi];
        }
        wallShearStress3.write();
        forAll(wallShearStress.boundaryField(), patchi)
        {
            wallShearStress.boundaryField()[patchi] =
                - shearStress.boundaryField()[patchi] / mesh.magSf().boundaryField()[patchi];
        }
        wallShearStress.write();
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
