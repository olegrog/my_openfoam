/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    emptyFoam

Description
    Solver that do nothing except functionObjects.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
class readFields
:
    public IOobjectList
{
    PtrList<GeoField> geoFieldList_;

public:
    readFields
    (
        const fvMesh& mesh
    )
    :
        IOobjectList(mesh, mesh.time().timeName())
    {
        IOobjectList fieldObjects(lookupClass(GeoField::typeName));
        forAllIter(IOobjectList, fieldObjects, iter)
        {
            Info<< "Read field " << iter()->name() << endl;
            geoFieldList_.append(
                new GeoField
                (
                    IOobject
                    (
                        iter()->name(),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
    }
};

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        readFields<volScalarField> scalarFields(mesh);
        readFields<volVectorField> vectorFields(mesh);
        readFields<volSphericalTensorField> sphericalTensorFields(mesh);
        readFields<volSymmTensorField> symmTensorFields(mesh);
        readFields<volTensorField> tensorFields(mesh);
        Info<< nl;

        runTime.functionObjects().execute();
    }
    
    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
