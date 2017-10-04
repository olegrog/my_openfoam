/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    dsmcCalcPressure

Description
    Calculates pressure (p) from temperature (overallT) and density (rhoMMean)
    Calculates Mach number (Ma) from velocity (UMean)

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
        
        if (objects.lookup("overallT")) {
            volScalarField T(*objects.lookup("overallT"), mesh); 
            if (!objects.lookup("p")) {
                Info<< "Writing p" << endl;
            } else {
                Info<< "Rewriting p" << endl;
            }
            volScalarField rho(*objects.lookup("rhoMMean"), mesh);
            volScalarField p("p", rho*T);
            p.write();

            forAll(T, celli) {
                if (T[celli] <= 1e-10) {
                    T[celli] = 1;
                }
            }
            volScalarField boundaryT(*objects.lookup("boundaryT"), mesh);
            T.boundaryFieldRef() = boundaryT.boundaryField();

            if (!objects.lookup("Ma")) {
                Info<< "Writing Ma" << endl;
                volVectorField U(*objects.lookup("UMean"), mesh);
                volScalarField Ma("Ma", mag(U)*sqrt(1.2/T));     // 2/gamma = 1.2 for monatomic gas
                Ma.write();
            }
        }

    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
