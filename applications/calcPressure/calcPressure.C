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
        
        if (!objects.lookup("Ma")) {
            Info<< "Writing Ma" << endl;
            volVectorField U(*objects.lookup("U"), mesh);
            volScalarField Ma("Ma", mag(U)*std::sqrt(1.2));     // 2/gamma = 1.2 for monatomic gas
            Ma.write();
        }

        volScalarField T(*objects.lookup("T"), mesh);
        if (!objects.lookup("p")) {
            Info<< "Writing p" << endl;
            volScalarField rho(*objects.lookup("rho"), mesh);
            volScalarField p("p", rho*T);
            p.write();
        }
        if (!objects.lookup("rho")) {
            Info<< "Writing rho" << endl;
            volScalarField p(*objects.lookup("p"), mesh);
            volScalarField rho("rho", p/T);
            rho.write();
        }

    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
