/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020-2021 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of gasMetalThermalProperties.

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

#include "generateGeometricField.H"
#include "geometricUniformField.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

Foam::scalar gasMetalAverage
(
    const Foam::Polynomial<2>& gasProperty,
    const Foam::Polynomial<2>& liquidProperty,
    const Foam::Polynomial<2>& solidProperty,
    Foam::scalar T,
    Foam::scalar liquidFraction,
    Foam::scalar gasFraction
)
{
    Foam::scalar solidFraction = 1 - liquidFraction - gasFraction;

    return
        solidProperty.value(T)*solidFraction
      + liquidProperty.value(T)*liquidFraction
      + gasProperty.value(T)*gasFraction;
}


} // End namespace

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T1, class T2, class T3>
Foam::tmp<T1> Foam::gasMetalThermo::Cp
(
    const T1& T,
    const T2& liquidFraction,
    const T3& gasFraction
) const
{
    const T2 sharpLiquidFraction
    (
        (1 - gasFraction)*pos(liquidFraction - (1 - gasFraction)/2)
    );

    return generateGeometricField<T1>
    (
        "Cp",
        mesh_,
        dimGasConstant,
        [this](scalar T, scalar phi, scalar alphaG)
        {
            return gasMetalAverage(gas_.Cp, liquid_.Cp, solid_.Cp, T, phi, alphaG);
        },
        T, sharpLiquidFraction, gasFraction
    );
}


template<class T1, class T2, class T3>
Foam::tmp<T1> Foam::gasMetalThermo::kappa
(
    const T1& T,
    const T2& liquidFraction,
    const T3& gasFraction
) const
{
    return generateGeometricField<T1>
    (
        "kappa",
        mesh_,
        dimPower/dimLength/dimTemperature,
        [this](scalar T, scalar phi, scalar alphaG)
        {
            return gasMetalAverage(gas_.kappa, liquid_.kappa, solid_.kappa, T, phi, alphaG);
        },
        T, liquidFraction, gasFraction
    );
}


template<class T1, class T2, class T3>
Foam::tmp<T3> Foam::gasMetalThermo::h
(
    const T1& T,
    const T2& liquidFraction,
    const T3& gasFraction,
    const word& name
) const
{
    return generateGeometricField<T3>
    (
        name,
        mesh_,
        dimEnergy/dimMass,
        [this](scalar T, scalar phi, scalar alphaG)
        {
            scalar alphaM = 1 - alphaG;
            scalar piecewise =
                T <= Tmelting_
              ? solid_.Cp.integral(0, T)
              : solid_.Cp.integral(0, Tmelting_) + liquid_.Cp.integral(Tmelting_, T);

            return alphaG*gas_.Cp.integral(0, T) + alphaM*piecewise
                + Hfusion_*phi;
        },
        T, liquidFraction, gasFraction
    );
}


template<class T1>
Foam::tmp<T1> Foam::gasMetalThermo::hAtMelting
(
    const T1& gasFraction
) const
{
    return h
    (
        geometricUniformField<scalar>(Tmelting_),
        ((1 - gasFraction)/2)(),
        gasFraction,
        "hAtMelting"
    );
}


template<class T1>
Foam::tmp<T1> Foam::gasMetalThermo::HsPrimeAlphaG
(
    const T1& T
) const
{
    return generateGeometricField<T1>
    (
        "HsPrimeAlphaG",
        mesh_,
        dimEnergy/dimMass,
        [this](scalar T)
        {
            scalar piecewise =
                T <= Tmelting_
              ? solid_.Cp.integral(Tmelting_, T)
              : liquid_.Cp.integral(Tmelting_, T);

            return gas_.Cp.integral(0, T) - solid_.Cp.integral(0, Tmelting_) - piecewise;
        },
        T
    );
}


template<class T1, class T2, class T3>
Foam::tmp<T1> Foam::gasMetalThermo::T
(
    const T1& h,
    const T1& hAtMelting,
    const T2& liquidFraction,
    const T3& gasFraction
) const
{
    return generateGeometricField<T1>
    (
        "T",
        mesh_,
        dimTemperature,
        [this](scalar h, scalar hAtMelting, scalar phi, scalar alphaG)
        {
            scalar alphaM = 1 - alphaG;
            scalar A = alphaG*gas_.Cp.derivative(0);
            scalar B = alphaG*gas_.Cp.value(0);
            scalar C = h - Hfusion_*phi - alphaM*solid_.Cp.integral(0, Tmelting_);

            A +=
                h < hAtMelting
              ? alphaM*solid_.Cp.derivative(0)
              : alphaM*liquid_.Cp.derivative(0);
            B +=
                h < hAtMelting
              ? alphaM*solid_.Cp.value(0)
              : alphaM*liquid_.Cp.value(0);
            C +=
                h < hAtMelting
              ? alphaM*solid_.Cp.integral(0, Tmelting_)
              : alphaM*liquid_.Cp.integral(0, Tmelting_);

            return
                mag(A) > SMALL
              ? (Foam::sqrt(sqr(B) + 2*A*C) - B)/A
              : C/B;
        },
        h, hAtMelting, liquidFraction, gasFraction
    );
}


// ************************************************************************* //
