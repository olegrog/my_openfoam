/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    calcPressure

Description
    Calculates pressure (p) from temperature (T) and density (rho)
    Calculates Mach number (Ma) from velocity (U)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();
        IOobjectList objects(mesh, runTime.timeName());

        volScalarField T(*objects.lookup("T"), mesh);
        volScalarField rho(*objects.lookup("rho"), mesh);
        volVectorField U(*objects.lookup("U"), mesh);

        volScalarField p("p", rho*T);
        // 2/gamma = 1.2 for monatomic gas
        volScalarField Ma("Ma", mag(U)*std::sqrt(1.2));

        p.write();
        Ma.write();
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
