/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019-2021 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of solidificationFoam.

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


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class T1>
auto Foam::multicomponentAlloy::sumSqrDelta(const T1& T) const -> decltype(T+T)
{
    auto iter = components_.begin();

    auto result = sqr(iter().delta(T));

    for (++iter; iter != components_.end(); ++iter)
    {
        result = result + sqr(iter().delta(T));
    }

    return result;
}

template<class T1>
auto Foam::multicomponentAlloy::sumSqrDeltaPerDL(const T1& T) const -> decltype(T+T)
{
    auto iter = components_.begin();

    auto result = sqr(iter().delta(T))/iter().phase("liquid").diffusion();

    for (++iter; iter != components_.end(); ++iter)
    {
        result = result + sqr(iter().delta(T))/iter().phase("liquid").diffusion();
    }

    return result;
}


template<class T1>
auto Foam::multicomponentAlloy::factor(const word& phaseName, const T1& T) const -> decltype(T+T)
{
    auto iter = components_.begin();

    auto result = iter().delta(T)/iter().phase(phaseName).slope(T);

    for (++iter; iter != components_.end(); ++iter)
    {
        result = result + iter().delta(T)/iter().phase(phaseName).slope(T);
    }

    return entropyChange_/result;
}


template<class T1>
auto Foam::multicomponentAlloy::factorPrime
(
    const word& phaseName,
    const T1& T
) const -> decltype(T+T)
{
    auto iter = components_.begin();

    auto result = iter().deltaPrime(T)/iter().phase(phaseName).slope(T);

    for (++iter; iter != components_.end(); ++iter)
    {
        result = result + iter().deltaPrime(T)/iter().phase(phaseName).slope(T);
    }

    return sqr(factor(phaseName, T))/entropyChange_*result;
}


// ************************************************************************* //
