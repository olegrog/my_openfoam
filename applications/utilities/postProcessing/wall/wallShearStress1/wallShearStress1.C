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

void writeWallField(
    const std::string name,
    const surfaceVectorField& field,
    const Foam::fvMesh& mesh,
    const Foam::Time& runTime
)
{
    dimensionedVector zeroField(
        "zeroForce",
        field.dimensions(),
        vector::zero
    );
    volVectorField wallField(
        IOobject("wall" + name, runTime.timeName(), mesh),
        mesh,
        zeroField
    );
    forAll(wallField.boundaryField(), patchi) {
        wallField.boundaryField()[patchi] =
            - field.boundaryField()[patchi] / mesh.magSf().boundaryField()[patchi];
    }
    wallField.write();
}

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

        const surfaceVectorField pressure0(
            fvc::interpolate(p) * mesh.Sf()
        );

        const surfaceVectorField pressure3(
            fvc::interpolate(
                -2. / 3 * g3 * fvc::laplacian(T, T)
            ) * mesh.Sf()
        );

        const surfaceVectorField pressure7(
            fvc::interpolate(
                g7 / 6 * magSqr(fvc::grad(T))
            ) * mesh.Sf()
        );

        const surfaceVectorField shearStress1(
            twoSymm(dev(fvc::interpolate(
                -g1 * sqrt(T) * fvc::grad(U)
            ))) & mesh.Sf()
        );

        // the second derivative of T is asymmetric due to interpolation
        const surfaceVectorField shearStress3(
            symm(dev(fvc::interpolate(
                g3 * T * fvc::grad(fvc::grad(T))
            ))) & mesh.Sf()
        );

        const surfaceVectorField shearStress7(
            dev(fvc::interpolate(
                g7 * fvc::grad(T) * fvc::grad(T)
            )) & mesh.Sf()
        );

        const surfaceVectorField force(
            pressure0 + pressure3 + pressure7 +
            shearStress1 + shearStress3 + shearStress7
        );

        Info<< "\nThe force acting on" << endl;
        forAll(force.boundaryField(), patchi)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchi]))
            {
                Info<< mesh.boundary()[patchi].name()
                    << " " << gSum(force.boundaryField()[patchi])
                    << " consisted of shear stresses:\n"
                    << "gamma_1 " << gSum(shearStress1.boundaryField()[patchi])
                    << ", gamma_3 " << gSum(shearStress3.boundaryField()[patchi])
                    << ", gamma_7 " << gSum(shearStress7.boundaryField()[patchi])
                    << " and hydrostatic pressures:\n"
                    << "p^dag " << gSum(pressure0.boundaryField()[patchi])
                    << ", gamma_3 " << gSum(pressure3.boundaryField()[patchi])
                    << ", gamma_7 " << gSum(pressure7.boundaryField()[patchi])
                    << endl;
            }
        }
        Info<< endl;

        writeWallField("Pressure0", pressure0, mesh, runTime);
        writeWallField("Pressure3", pressure3, mesh, runTime);
        writeWallField("Pressure7", pressure7, mesh, runTime);
        writeWallField("ShearStress1", shearStress1, mesh, runTime);
        writeWallField("ShearStress3", shearStress3, mesh, runTime);
        writeWallField("ShearStress7", shearStress7, mesh, runTime);
        writeWallField("Force", force, mesh, runTime);
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
