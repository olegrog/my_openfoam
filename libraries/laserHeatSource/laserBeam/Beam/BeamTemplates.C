/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2022 Oleg Rogozin
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

#include "Beam.H"

#include "generateGeometricField.H"
#include "laserProperties.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::Beam<Type>::I(const vector& x) const
{
    const scalar P = laser_.power().value();
    const scalar R = laser_.radius().value();
    const vector x0 = laser_.position().value();

    return P*Type::intensity(x, x0, direction_, R);
}


template<class Type>
Foam::tmp<Foam::volScalarField> Foam::Beam<Type>::I() const
{
    const scalar P = laser_.power().value();
    const scalar R = laser_.radius().value();
    const vector x0 = laser_.position().value();
    const volVectorField& x = mesh_.C();

    return generateGeometricField<volScalarField>
    (
        "BeamIntensity",
        mesh_,
        dimPower/dimArea,
        [=](vector x)
        {
            return P*Type::intensity(x, x0, direction_, R);
        },
        x
    );
}


// ************************************************************************* //
