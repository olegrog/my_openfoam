/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
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

Description
    Set initial conditions for alpha field, which represent the powder bed on
    a substrate.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvCFD.H"      // for tmp<...>

tmp<volScalarField> generateBall(
    const volVectorField& coord,
    const dimensionedVector& center,
    const dimensionedScalar& radius)
{
    return min(
        2 * Foam::exp(-magSqr((coord - center) / radius)),
        scalar(1)
    );
}

int main(int argc, char *argv[])
{
    argList::validArgs.append("field");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    word alphaName = args.argRead<word>(1);

    // -- Read an alpha field
    Info<< "Reading field " << alphaName << endl;
    volScalarField alpha
    (
        IOobject
        (
            alphaName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // -- Read a dictionary
    const word dictName("powderBedProperties");
    Info<< "Reading " << dictName << endl;
    IOdictionary powderBedProperties(
        IOobject(
            dictName,
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar ballRadius(powderBedProperties.lookup("ballRadius"));
    dimensionedScalar substratePosition(powderBedProperties.lookup("substratePosition"));

    label seed(powderBedProperties.lookupOrDefault<label>("seed", 0));
    scalar amplitudeRadius(powderBedProperties.lookupOrDefault<scalar>("amplitudeRadius", 0));
    scalar amplitudePosition(powderBedProperties.lookupOrDefault<scalar>("amplitudePosition", 0));
    scalar latticeStep(powderBedProperties.lookupOrDefault<scalar>("latticeStep", 2));

    Random random(seed);
    for (int j = -2; j <= 2; j++) {
        for (int i = -2; i < 15-2; i++) {
            const dimensionedScalar& R = ballRadius * (1 + amplitudeRadius * random.scalarAB(-1, 1));
            alpha += generateBall(mesh.C(), dimensionedVector("center", dimless, vector(
                amplitudePosition * random.scalarAB(-1, 1) + latticeStep*i,
                amplitudePosition * random.scalarAB(-1, 1) + latticeStep*j,
                (substratePosition/R).value() + 1
            )) * R, R);
        }
    }
    //alpha = min(max(alpha, scalar(0)), scalar(1));

    Info<< "Writing field" << alphaName << endl;
    alpha.write();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
