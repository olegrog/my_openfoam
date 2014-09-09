/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    wallSnitForces

Description
    Calculates and reports the forces acting on the wall for all patches and
    the specific times.

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
            field.boundaryField()[patchi] / mesh.magSf().boundaryField()[patchi];
    }
    wallField.write();
}

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();
        #include "createFields.H"

        // NB: mesh.Sf() is a normal, pointed out of the gas on the boundary
        const surfaceVectorField pressure0(
            fvc::interpolate(2 * p) * mesh.Sf()
        );
        // we avoid laplacian scheme for benefit of more precise gradient one
        const surfaceVectorField pressure3(
            fvc::interpolate(
                //-2 * g3 / 3 * fvc::laplacian(T, T)
                /**/
                - g3 / 3 * (
                    magSqr(fvc::grad(T)) + 4 * (U & fvc::grad(T)) / (g2 * sqrt(T))
                )
                /**/
            ) * mesh.Sf()
        );
        const surfaceVectorField pressure7(
            fvc::interpolate(
                g7 / 6 * magSqr(fvc::grad(T))
            ) * mesh.Sf()
        );
        const surfaceVectorField shearStress1(
            //- g1 * (mesh.Sf() & fvc::interpolate(sqrt(T) * fvc::grad(U)))
            //- g1 * fvc::interpolate(sqrt(T)) * fvc::snGrad(U) * mesh.magSf()
            /**/
            twoSymm(dev(fvc::interpolate(
                -g1 * sqrt(T) * fvc::grad(U)
            ))) & mesh.Sf()
            /**/
        );
        const surfaceVectorField shearStress3(
            //g3 * fvc::snGrad(T) * (fvc::interpolate(T * curvature) - fvc::snGrad(T) / 3) * mesh.Sf()
            /**/
            // the second derivative of T can be asymmetric due to interpolation
            symm(dev(fvc::interpolate(
                g3 * T * fvc::grad(fvc::grad(T))
            ))) & mesh.Sf()
            /**/
        );
        const surfaceVectorField shearStress7(
            //2 * g7 / 3 * sqr(fvc::snGrad(T)) * mesh.Sf()
            /**/
            dev(fvc::interpolate(
                g7 * sqr(fvc::grad(T))
            )) & mesh.Sf()
            /**/
        );
        const surfaceVectorField force(
            pressure0 + pressure3 + pressure7 +
            shearStress1 + shearStress3 + shearStress7
        );
        const surfaceVectorField moment(
            mesh.Cf() ^ force
        );

        Info<< "\nThe force acting on" << endl;
        forAll(force.boundaryField(), patchi)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchi]))
            {
                Info<< mesh.boundary()[patchi].name()
                    << " " << gSum(force.boundaryField()[patchi])
                    << " consisted of\n\tshear stresses: "
                    << "gamma_1 " << gSum(shearStress1.boundaryField()[patchi])
                    << ", gamma_3 " << gSum(shearStress3.boundaryField()[patchi])
                    << ", gamma_7 " << gSum(shearStress7.boundaryField()[patchi])
                    << ",\n\tand pressures: "
                    << "p^dag " << gSum(pressure0.boundaryField()[patchi])
                    << ", gamma_3 " << gSum(pressure3.boundaryField()[patchi])
                    << ", gamma_7 " << gSum(pressure7.boundaryField()[patchi])
                    << ".\n\tmoment of force relative to origin: "
                    << gSum(moment.boundaryField()[patchi])
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
        writeWallField("Moment", moment, mesh, runTime);

        const volVectorField DP(fvc::grad(p));
        const volVectorField DT(fvc::grad(T));
        DP.write();
        DT.write();
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
