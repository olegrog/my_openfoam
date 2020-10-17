/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020 Oleg Rogozin
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
    stereolithography

Description
    Calculate degree of polymerization after multi-layer stereolithography
    of polymer--ceramic suspension.
    The main laser scanning direction correspond to the x axis.
    The build direction coincides with the z axis.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "laserScanner.H"

using constant::mathematical::pi;

tmp<volScalarField> BeerLambert(const volScalarField& Z, const dimensionedScalar Dp)
{
    return neg0(Z)*exp(Z/Dp);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"
    #include "createFields.H"

    const volVectorField coord = mesh.C();
    const volScalarField Z = coord.component(2);

    const boundBox& bounds = mesh.bounds();
    const dimensionedScalar ymin("ymin", dimLength, bounds.min().y());
    const dimensionedScalar ymax("ymax", dimLength, bounds.max().y());
    const dimensionedScalar zmin("zmin", dimLength, bounds.min().z());
    const dimensionedScalar zmax("zmax", dimLength, bounds.max().z());
    const dimensionedScalar L = ymax - ymin;
    const dimensionedScalar small("small", dimLength, SMALL);

    const scalar Fscatt = aConst*log(laser.E()/Ec).value() + bConst;
    const dimensionedScalar PConst =
        -pow025(2*pi)*Foam::log(1 - pGel)/sqrt(Ec*laser.R()/laser.V(L));
    const dimensionedScalar radius = laser.R()*Fscatt;

    Info<< " -- Scattering factor = " << Fscatt << endl;
    Info<< " -- Polymerization constant = " << PConst.value() << endl;

    if (notEqual(zmin.value(), 0))
    {
        FatalError
            << "Coordinate zmin is not equal to zero."
            << exit(FatalError);
    }

    Info<< "Calculating totalExposure" << endl;
    label i;
    for (i = 1; laser.height(i) < zmax + small; i++)
    {
        volScalarField exposure = laser.E()*BeerLambert(Z - laser.height(i) - small, Dp);
        totalSqrtExposure += sqrt(exposure);
        totalExposure += exposure;
    }

    Info << " -- Total number of layers = " << --i << endl;

    if (notEqual(zmax.value(), laser.height(i).value()))
    {
        FatalError
            << "Coordinate zmax does not correspond to integer number of layers:"
            << "  zmax = " << zmax.value() << nl
            << "  laser.height(" << i << ") = " << laser.height(i).value()
            << exit(FatalError);
    }

    Info<< "Calculating degree of polymerization" << endl;
    polymerization = 1 - exp(-pow025(2*pi)*PConst*totalSqrtExposure*sqrt(radius/laser.V(L)));

    Info<< "Writing fields" << endl;
    totalExposure.write();
    totalSqrtExposure.dimensions().clear(); // ParaView cannot read fractional dimensions
    totalSqrtExposure.write();
    polymerization.write();

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
