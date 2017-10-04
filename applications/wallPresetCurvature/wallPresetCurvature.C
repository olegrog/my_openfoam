/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    wallPresetCurvature

Description
    Calculates and presets the curvature of wall patches for the specific
    times.

\*---------------------------------------------------------------------------*/

#include <cstdlib>
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkOption(Foam::argList args, word option)
{
    if (!args.optionFound(option)) {
        FatalErrorIn("int main()")
            << "The " << option << " has not been specified."
            << abort(FatalError);
    }
}

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    timeSelector::addOptions();

    argList::addBoolOption(
        "planar",
        "consider boundary as planar"
    );
    argList::addBoolOption(
        "cylinder",
        "consider boundary as out of cylinder"
    );
    argList::addBoolOption(
        "sphere",
        "consider boundary as out of sphere"
    );
    argList::addBoolOption(
        "elliptic",
        "consider boundary as out of elliptic cylinder"
    );

    argList::addOption(
        "patch",
        "word",
        "patch name"
    );
    argList::addOption(
        "radius",
        "scalar",
        "specific radius of cylinder or sphere"
    );
    argList::addOption(
        "center",
        "vector",
        "specific center (for elliptic)"
    );
    argList::addOption(
        "major",
        "scalar",
        "major semiaxis (for elliptic)"
    );
    argList::addOption(
        "minor",
        "scalar",
        "minor semiaxis (for elliptic)"
    );
    argList::addOption(
        "phi",
        "scalar",
        "angle of rotation (part of Pi) (for elliptic)"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    volScalarField curvature
    (
        IOobject
        (
            "curvature",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        mesh,
        dimensionedScalar("zero", dimless/dimLength, 0)
    );

    checkOption(args, "patch");
    label patch = mesh.boundary().findPatchID(args["patch"]);
    if (patch == -1) {
        FatalErrorIn("int main()")
            << "The patch has not been found."
            << abort(FatalError);
    } else {
        Info<< "Updating patch: " << args["patch"] << endl;
    }

    if (args.optionFound("cylinder")) {
        checkOption(args, "radius");
        scalar value = 1./strtod(args["radius"].c_str(), NULL);
        Info<< "Preset cylindrical curvature = " << value << endl;
        curvature.boundaryFieldRef()[patch] = value;
    }
    if (args.optionFound("sphere")) {
        checkOption(args, "radius");
        scalar value = 2./strtod(args["radius"].c_str(), NULL);
        Info<< "Preset spherical curvature = " << value << endl;
        curvature.boundaryFieldRef()[patch] = value;
    }
    if (args.optionFound("elliptic")) {
        checkOption(args, "major");
        checkOption(args, "minor");
        checkOption(args, "phi");
        scalar major = strtod(args["major"].c_str(), NULL);
        scalar minor = strtod(args["minor"].c_str(), NULL);
        scalar phi = strtod(args["phi"].c_str(), NULL);
        forAll(curvature.boundaryFieldRef()[patch], celli) {
            scalar x = mesh.boundary()[patch].Cf()[celli].x();
            scalar y = mesh.boundary()[patch].Cf()[celli].y();
            scalar pi = constant::mathematical::pi;
            scalar x_ = x * Foam::sin(phi*pi) + y * Foam::cos(phi*pi);
            scalar y_ = y * Foam::sin(phi*pi) - x * Foam::cos(phi*pi);
            curvature.boundaryFieldRef()[patch][celli] =
                sign(major) * pow4(major*minor) /
                pow3(Foam::sqrt(pow4(minor)*x_*x_ + pow4(major)*y_*y_));
        }
        Info<< "Preset elliptical curvature"
            << ": min = " << min(curvature.boundaryFieldRef()[patch])
            << ", max = " << max(curvature.boundaryFieldRef()[patch])
            << endl;
    }

    Info<< endl;
    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        curvature.write();
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
