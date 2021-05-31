/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019-2020 Oleg Rogozin
-------------------------------------------------------------------------------
License
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    generatePowderBed

Description
    Set initial conditions for alpha field, which represent the powder bed on
    a substrate.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cutCellIso.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using constant::mathematical::pi;

scalarField generateBall
(
    const vectorField& points,
    const vector center,
    const scalar radius
)
{
    return mag(points - center) - radius;
}


tmp<volScalarField> VolumeOfFluid(const fvMesh& mesh, scalarField& f)
{
    cutCellIso cutCell(mesh, f);
    auto tres = volScalarField::New("result", mesh, dimless);
    auto& res = tres.ref();

    forAll(res, cellI)
    {
        label cellStatus = cutCell.calcSubCell(cellI, Zero);

        if (cellStatus == -1)
        {
            res[cellI] = 1;
        }
        else if (cellStatus == 1)
        {
            res[cellI] = 0;
        }
        else if (cellStatus == 0)
        {
            if (mag(cutCell.faceArea()) != 0)
            {
                res[cellI] = max(min(cutCell.VolumeOfFluid(), 1), 0);
            }
        }
    }

    return tres;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::addArgument("field");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    word alphaName = args.get<word>(1);

    Info<< "Reading field " << alphaName << endl;
    volScalarField alpha
    (
        IOobject
        (
            alphaName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        ),
        mesh
    );

    const word dictName("powderBedProperties");
    Info<< "Reading " << dictName << endl;
    IOdictionary dict
    (
        IOobject
        (
            dictName,
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Dimensionless parameters of the powder distribution
    const label seed(dict.lookupOrDefault<label>("seed", 0));
    const scalar amplitudeRadius(dict.lookupOrDefault<scalar>("amplitudeRadius", 0));
    const scalar amplitudePosition(dict.lookupOrDefault<scalar>("amplitudePosition", 0));
    const scalar latticeStep(dict.lookupOrDefault<scalar>("latticeStep", 2));

    // Dimensioned parameters
    const dimensionedScalar ballRadius("ballRadius", dict);
    const dimensionedScalar substratePosition("substratePosition", dict);

    // Auxiliary constants
    const vector substrateNormal(0, 0, 1);
    const scalar domainVolume = gSum(mesh.V());
    const boundBox& bounds = mesh.bounds();

    scalar volumeFraction = 0;  // for theoretical prediction

    // --- Generate substrate
    {
        scalarField f = substratePosition.value() - (mesh.points() & substrateNormal);
        alpha = VolumeOfFluid(mesh, f);
        volumeFraction = alpha.weightedAverage(mesh.Vsc()).value();
    }

    // --- Generate balls
    Random random(seed);
    for (int j = -2; j <= 2; j++)
    {
        for (int i = -2; i < 15-2; i++)
        {
            const scalar R = ballRadius.value()*(1 + amplitudeRadius*random.position(-1, 1));
            const scalar X = (amplitudePosition*random.position(-1, 1) + latticeStep*i)*R;
            const scalar Y = (amplitudePosition*random.position(-1, 1) + latticeStep*j)*R;
            const scalar Z = substratePosition.value() + R;

            scalarField f = -generateBall(mesh.points(), vector(X, Y, Z), R);
            alpha += VolumeOfFluid(mesh, f);

            // Evaluate the volume of ball (full or cut)
            // TODO(olegrog): improve accuracy by adding spherical caps
            if
            (
                bounds.min().x() < X && X < bounds.max().x()
             && bounds.min().y() < Y && Y < bounds.max().y()
            )
            {
                volumeFraction += 4./3*pi*pow(R, 3)/domainVolume;
            }
        }
    }

    // --- Save the result
    alpha.clip(0, 1);
    alpha.correctBoundaryConditions();
    Info<< "Writing field " << alphaName << endl;
    ISstream::defaultPrecision(18);
    alpha.write();

    // --- Analyze the result
    Info<< nl << alphaName << ": volume fraction = "
        << alpha.weightedAverage(mesh.Vsc()).value()
        << " theoretical = " << volumeFraction
        << nl << endl;

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
