/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
icoReactingMultiphaseInterFoam | Copyright (C) 2017-2019 OpenCFD Ltd
            slmMeltPoolFoam | Copyright (C) 2021 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of slmMeltPoolFoam.

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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::rayTracingParticle::trackingData::trackingData
(
    Cloud<rayTracingParticle>& spc,
    const interpolationCell<scalar>& AInterp,
    const interpolationCellPoint<vector>& gradAlphaMInterp,
    const boolField& reflectedCell,
    const reflectionModel& reflection,
    volScalarField& Q
)
:
    particle::trackingData(spc),
    AInterp_(AInterp),
    gradAlphaMInterp_(gradAlphaMInterp),
    reflectingCells_(reflectedCell),
    reflection_(reflection),
    Q_(Q)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::interpolationCell<Foam::scalar>&
Foam::rayTracingParticle::trackingData::AInterp() const
{
    return AInterp_;
}


inline const Foam::interpolationCellPoint<Foam::vector>&
Foam::rayTracingParticle::trackingData::gradAlphaMInterp() const
{
    return gradAlphaMInterp_;
}

inline const Foam::boolField&
Foam::rayTracingParticle::trackingData::reflectingCells() const
{
    return reflectingCells_;
}


inline const Foam::reflectionModel&
Foam::rayTracingParticle::trackingData::reflection() const
{
    return reflection_;
}


inline Foam::scalar& Foam::rayTracingParticle::trackingData::Q(label celli)
{
    return Q_[celli];
}


inline const Foam::point& Foam::rayTracingParticle::p0() const
{
    return p0_;
}


inline const Foam::point& Foam::rayTracingParticle::p1() const
{
    return p1_;
}


inline Foam::scalar Foam::rayTracingParticle::I0() const
{
    return I0_;
}


inline Foam::scalar Foam::rayTracingParticle::I() const
{
    return I_;
}


inline Foam::scalar Foam::rayTracingParticle::dA() const
{
    return dA_;
}


inline Foam::scalar& Foam::rayTracingParticle::dA()
{
    return dA_;
}


inline Foam::point& Foam::rayTracingParticle::p1()
{
    return p1_;
}


inline Foam::scalar& Foam::rayTracingParticle::I0()
{
    return I0_;
}


inline Foam::scalar& Foam::rayTracingParticle::I()
{
    return I_;
}


// ************************************************************************* //
