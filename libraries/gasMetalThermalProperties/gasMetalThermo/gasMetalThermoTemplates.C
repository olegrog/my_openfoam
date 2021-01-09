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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
Foam::tmp<Foam::volScalarField> Foam::gasMetalThermo::Cp
(
    const T1& T,
    const T2& liquidFraction,
    const T3& gasFraction
) const
{
    const T2 sharpLiquidFraction
    (
        pos(liquidFraction - (1 - gasFraction)/2 - VSMALL)
    );

    return generateGeometricField<volScalarField>
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
Foam::tmp<Foam::volScalarField> Foam::gasMetalThermo::k
(
    const T1& T,
    const T2& liquidFraction,
    const T3& gasFraction
) const
{
    return generateGeometricField<volScalarField>
    (
        "k",
        mesh_,
        dimPower/dimLength/dimTemperature,
        [this](scalar T, scalar phi, scalar alphaG)
        {
            return gasMetalAverage(gas_.k, liquid_.k, solid_.k, T, phi, alphaG);
        },
        T, liquidFraction, gasFraction
    );
}


template<class T1, class T2, class T3>
Foam::tmp<Foam::volScalarField> Foam::gasMetalThermo::h
(
    const T1& T,
    const T2& liquidFraction,
    const T2& vapourFraction,
    const T3& gasFraction,
    const word& name
) const
{
    return generateGeometricField<volScalarField>
    (
        name,
        mesh_,
        dimEnergy/dimMass,
        [this](scalar T, scalar phi, scalar psi, scalar alphaG)
        {
            scalar alphaM = 1 - alphaG;
            scalar piecewise =
                T <= Tmelting_
              ? solid_.Cp.integral(0, T)
              : solid_.Cp.integral(0, Tmelting_) + liquid_.Cp.integral(Tmelting_, T);

            return alphaG*gas_.Cp.integral(0, T) + alphaM*piecewise
                + Hfusion_*phi + Hvapour_*psi;
        },
        T, liquidFraction, vapourFraction, gasFraction
    );
}


template<class T1>
Foam::tmp<Foam::volScalarField> Foam::gasMetalThermo::hAtMelting
(
    const T1& gasFraction
) const
{
    return h
    (
        geometricUniformField<scalar>(Tmelting_),
        geometricUniformField<scalar>(0.5),
        geometricUniformField<scalar>(0),
        gasFraction, "hAtMelting"
    );
}


template<class T1>
Foam::tmp<Foam::volScalarField> Foam::gasMetalThermo::hAtBoiling
(
    const T1& gasFraction
) const
{
    return h
    (
        geometricUniformField<scalar>(Tboiling_),
        geometricUniformField<scalar>(1),
        geometricUniformField<scalar>(0.5),
        gasFraction, "hAtBoiling"
    );
}


template<class T1>
Foam::tmp<Foam::volScalarField> Foam::gasMetalThermo::sensibleEnthalpyPrimeGasFraction
(
    const T1& T
) const
{
    return generateGeometricField<volScalarField>
    (
        "sensibleEnthalpyPrimeGasFraction",
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
Foam::tmp<Foam::volScalarField> Foam::gasMetalThermo::T
(
    const T1& h,
    const T1& hAtMelting,
    const T2& liquidFraction,
    const T2& vapourFraction,
    const T3& gasFraction
) const
{
    return generateGeometricField<volScalarField>
    (
        "T",
        mesh_,
        dimTemperature,
        [this](scalar h, scalar hAtMelting, scalar phi, scalar psi, scalar alphaG)
        {
            scalar alphaM = 1 - alphaG;
            scalar A = alphaG*gas_.Cp.derivative(0);
            scalar B = alphaG*gas_.Cp.value(0);
            scalar C = h - Hfusion_*phi - Hvapour_*psi - alphaM*solid_.Cp.integral(0, Tmelting_);

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
        h, hAtMelting, liquidFraction, vapourFraction, gasFraction
    );
}


template<class T1, class T2>
Foam::tmp<Foam::volScalarField> Foam::gasMetalThermo::metalPhaseFraction
(
    const T1& h,
    const T1& hAtPhaseTransition,
    const T2& metalFraction,
    const scalar latentHeat,
    const word& phaseName
) const
{
    return generateGeometricField<volScalarField>
    (
        phaseName + "Fraction",
        mesh_,
        dimless,
        [latentHeat](scalar h, scalar hAtPhaseTransition, scalar alphaM) -> scalar
        {
            scalar h1 = hAtPhaseTransition - alphaM*latentHeat/2;
            scalar h2 = hAtPhaseTransition + alphaM*latentHeat/2;

            if (h < h1) {
                return 0;
            } else if (h > h2) {
                return alphaM;
            } else {
                return (h - h1)/latentHeat;
            }
        },
        h, hAtPhaseTransition, metalFraction
    );
}


template<class T1, class T2>
Foam::tmp<Foam::volScalarField> Foam::gasMetalThermo::metalPhaseFractionPrimeEnthalpy
(
    const T1& h,
    const T1& hAtPhaseTransition,
    const T2& metalFraction,
    const scalar latentHeat,
    const word& phaseName
) const
{
    return generateGeometricField<volScalarField>
    (
        phaseName + "FractionPrimeEnthalpy",
        mesh_,
        dimMass/dimEnergy,
        [latentHeat](scalar h, scalar hAtPhaseTransition, scalar alphaM) -> scalar
        {
            scalar h1 = hAtPhaseTransition - alphaM*latentHeat/2;
            scalar h2 = hAtPhaseTransition + alphaM*latentHeat/2;

            if (h < h1 || h > h2) {
                return 0;
            } else {
                return 1/latentHeat;
            }
        },
        h, hAtPhaseTransition, metalFraction
    );
}


template<class T1, class T2>
Foam::tmp<Foam::volScalarField> Foam::gasMetalThermo::metalPhaseFractionPrimeGasFraction
(
    const T1& h,
    const T1& hAtPhaseTransition,
    const T2& metalFraction,
    const scalar latentHeat,
    const scalar hPrePhaseTransitionPrime,
    const word& phaseName
) const
{
    return generateGeometricField<volScalarField>
    (
        phaseName + "FractionPrimeGasFraction",
        mesh_,
        dimless,
        [=](scalar h, scalar hAtPhaseTransition, scalar alphaM) -> scalar
        {
            scalar h1 = hAtPhaseTransition - alphaM*latentHeat/2;
            scalar h2 = hAtPhaseTransition + alphaM*latentHeat/2;

            if (h < h1) {
                return 0;
            } else if (h > h2) {
                return -1;
            } else {
                return -hPrePhaseTransitionPrime/latentHeat;
            }
        },
        h, hAtPhaseTransition, metalFraction
    );
}


// ************************************************************************* //
