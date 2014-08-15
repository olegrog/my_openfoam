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

        const surfaceVectorField shearStress
        (
            symm(dev(fvc::interpolate
            (
                - g1 * 2 * sqrt(T) * fvc::grad(U)
                + g7 * fvc::grad(T) * fvc::grad(T)
                + g3 * T * fvc::grad(fvc::grad(T))
            ))) & mesh.Sf()
        );

        const surfaceVectorField::GeometricBoundaryField& patchShearStress =
            shearStress.boundaryField();

        Info<< "\nThe force acting on" << endl;
        forAll(shearStress.boundaryField(), patchi)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchi]))
            {
                Info<< mesh.boundary()[patchi].name()
                    << " "
                    << gSum(patchShearStress[patchi])
                    << endl;
            }
        }
        Info<< endl;

        volVectorField wallShearStress
        (
            IOobject
            (
                "wallShearStress",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedVector(
                "wallShearStress",
                shearStress.dimensions(),
                vector::zero
            )
        );

        forAll(wallShearStress.boundaryField(), patchi)
        {
            wallShearStress.boundaryField()[patchi] =
                - patchShearStress[patchi] / mesh.magSf().boundaryField()[patchi];
        }
        wallShearStress.write();
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
