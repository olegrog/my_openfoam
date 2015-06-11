/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    correctBC

Description
    Execute correctBoundaryConditions() for all fields.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void correctFields
(
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    IOobjectList fieldObjects(objects.lookupClass(GeoField::typeName));
    forAllIter(IOobjectList, fieldObjects, iter)
    {
        if (selectedFields.empty() || selectedFields.found(iter()->name()))
        {
            GeoField field(*iter(), mesh);
            field.correctBoundaryConditions();
            field.write();
        }
    }
}

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::addOption
    (
        "fields",
        "wordList",
        "only convert the specified fields - e.g. '(p T U)'"
    );
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    HashSet<word> selectedFields;
    args.optionReadIfPresent("fields", selectedFields);

    instantList timeDirs = timeSelector::select0(runTime, args);
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();
        
        IOobjectList objects(mesh, runTime.timeName());

        correctFields<volScalarField>(mesh, objects, selectedFields);
        correctFields<volVectorField>(mesh, objects, selectedFields);
        correctFields<volSphericalTensorField>(mesh, objects, selectedFields);
        correctFields<volSymmTensorField>(mesh, objects, selectedFields);
        correctFields<volTensorField>(mesh, objects, selectedFields);

    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
