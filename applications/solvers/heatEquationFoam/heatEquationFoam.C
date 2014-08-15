/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.2.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
-------------------------------------------------------------------------------
Application
    heatEquationFoam

Description
    Solves a heat equation for hard-sphere gas.
    The thermal conductivity is proportional to the square root of the
    temperature.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    simpleControl simple(mesh);

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        solve (
            fvm::laplacian (sqrt(T), T)
        );

    	runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
