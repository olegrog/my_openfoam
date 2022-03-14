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

#include "error.H"
#include "sigFpe.H"

#include "generateGeometricField.H"
#include "laserProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Beam<Type>::Beam
(
    const dictionary& dict,
    const fvMesh& mesh,
    const laserProperties& laser
)
:
    laserBeam(dict, mesh, laser)
{
    const bool sigActive = sigFpe::active();
    // SIGFPE can be trapped, but here we have to disable this feature
    if (sigActive) sigFpe::unset();
    const List<Pair<scalar>> actualVsExpected
    {
        { Type::intensity(VGREAT), 0},
        { Type::integral(0), 0},
        { Type::integral(VGREAT), 1},
    };
    if (sigActive) sigFpe::set();

    for (const auto& pair : actualVsExpected)
    {
        if (notEqual(pair.first(), pair.second()))
        {
            FatalErrorInFunction
                << "Wrong implemented laser beam intensity model since "
                << pair.first() << " != " << pair.second()
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::Beam<Type>::I(scalar r) const
{
    using constant::mathematical::pi;

    const scalar P = laser_.power().value();
    const scalar R = laser_.radius().value();

    return P*Type::intensity(r/R)/pi/sqr(R);
}


template<class Type>
Foam::scalar Foam::Beam<Type>::integralI(scalar r) const
{
    const scalar P = laser_.power().value();
    const scalar R = laser_.radius().value();

    return P*Type::integral(r/R);
}


template<class Type>
Foam::tmp<Foam::volScalarField> Foam::Beam<Type>::I() const
{
    using constant::mathematical::pi;

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
            const scalar r = mag((x - x0) ^ direction_)/R;
            return P*Type::intensity(r)/pi/sqr(R);
        },
        x
    );
}


// ************************************************************************* //
