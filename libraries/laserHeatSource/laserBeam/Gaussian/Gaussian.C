/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2021 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of laserHeatSource.

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

\*---------------------------------------------------------------------------*/

#include "Gaussian.H"

#include "addToRunTimeSelectionTable.H"
#include "constants.H"

#include "generateGeometricField.H"
#include "laserProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

Foam::scalar intensity
(
    Foam::scalar P,
    Foam::scalar R,
    const Foam::vector& x0,
    const Foam::vector& x,
    const Foam::vector& n
)
{
    using namespace Foam;
    using constant::mathematical::piByTwo;

    return P*Foam::exp(-2*magSqr(((x - x0) ^ n)/R))/piByTwo/sqr(R);
}


} // End namespace

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace beam
    {
        defineTypeName(Gaussian);
        addToRunTimeSelectionTable(laserBeam, Gaussian, dictionary);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::beam::Gaussian::I(const vector& x) const
{
    const scalar P = laser_.power().value();
    const scalar R = laser_.radius().value();
    const vector x0 = laser_.position().value();

    return intensity(P, R, x0, x, direction_);
}


Foam::tmp<Foam::volScalarField> Foam::beam::Gaussian::I() const
{
    const scalar P = laser_.power().value();
    const scalar R = laser_.radius().value();
    const vector x0 = laser_.position().value();
    const volVectorField& x = mesh_.C();

    return generateGeometricField<volScalarField>
    (
        "I",
        mesh_,
        dimPower/dimArea,
        [=](vector x)
        {
            return intensity(P, R, x0, x, direction_);
        },
        x
    );
}


// ************************************************************************* //
