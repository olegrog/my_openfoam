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
        wallField.boundaryFieldRef()[patchi] =
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

        const surfaceVectorField hydrostatic(
            2 * fvc::interpolate(p2/p0) * mesh.Sf()
        );
        const surfaceVectorField viscous(
            g1*fvc::interpolate(pow(T0, s)/p0) * (
                - (s + 1./3) * fvc::interpolate(fvc::div(U1)) * mesh.Sf()
                - fvc::snGrad(U1) * mesh.magSf()
            )
        );
        const surfaceVectorField thermal(
            - g7 * fvc::interpolate(pow(T0, 2*s-1)*fvc::grad(T0)/p0) * fvc::snGrad(T0) * mesh.magSf()
            + g7/2 * magSqr(fvc::interpolate(fvc::grad(T0)))*fvc::interpolate(pow(T0, 2*s-1)/p0) * mesh.Sf()
        );
        const surfaceVectorField force(
            hydrostatic + viscous + thermal
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
                    << " consisted of:\n"
                    << "\thydrostatic " << gSum(hydrostatic.boundaryField()[patchi])
                    << ", viscous " << gSum(viscous.boundaryField()[patchi])
                    << ", thermal " << gSum(thermal.boundaryField()[patchi])
                    << ",\n\tmoment of force relative to origin: "
                    << gSum(moment.boundaryField()[patchi])
                    << "." << endl;
            }
        }
        Info<< endl;

        writeWallField("HydrostaticForce", hydrostatic, mesh, runTime);
        writeWallField("ViscousForce", viscous, mesh, runTime);
        writeWallField("ThermalForce", thermal, mesh, runTime);
        writeWallField("Force", force, mesh, runTime);
        writeWallField("Moment", moment, mesh, runTime);

        const volVectorField DP(fvc::grad(p2));
        const volVectorField DT(fvc::grad(T0));
        DP.write();
        DT.write();
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
