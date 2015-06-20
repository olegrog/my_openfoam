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
class createFields
:
    public IOobjectList
{
    PtrList<GeoField> geoFieldList_;

public:
    createFields
    (
        const fvMesh& mesh
    )
    :
        IOobjectList(mesh, mesh.time().timeName())
    {
        IOobjectList fieldObjects(lookupClass(GeoField::typeName));
        forAllIter(IOobjectList, fieldObjects, iter)
        {
            Info<< "Create field " << iter()->name() << endl;
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
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    createFields<volScalarField> scalarFields(mesh);
    createFields<volVectorField> vectorFields(mesh);
    createFields<volSphericalTensorField> sphericalTensorFields(mesh);
    createFields<volSymmTensorField> symmTensorFields(mesh);
    createFields<volTensorField> tensorFields(mesh);
    Info<< nl;

    while(runTime.loop()) {
        Info<< "Time = " << runTime.timeName() << nl << endl;
    }
    
    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
