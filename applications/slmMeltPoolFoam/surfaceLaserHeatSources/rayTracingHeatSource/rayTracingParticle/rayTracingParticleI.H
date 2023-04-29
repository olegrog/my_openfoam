/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
               DTRMParticle | Copyright (C) 2017-2019 OpenCFD Ltd
         rayTracingParticle | Copyright (C) 2021-2023 Oleg Rogozin
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
    Cloud<rayTracingParticle>& cloud,
    bool useSubCellData,
    const reconstructionSchemes& surf,
    const scatteringModel& scattering,
    const volScalarField& alphaM,
    const volVectorField& gradAlphaM,
    volScalarField& Q
)
:
    particle::trackingData(cloud),
    useSubCellData_(useSubCellData),
    surf_(surf),
    scatteringModel_(scattering),
    alphaM_(alphaM),
    gradAlphaM_(gradAlphaM),
    Q_(Q),
    markedCells_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool
Foam::rayTracingParticle::trackingData::useSubCellData() const
{
    return useSubCellData_;
}


inline const Foam::reconstructionSchemes&
Foam::rayTracingParticle::trackingData::surf() const
{
    return surf_;
}


inline const Foam::scatteringModel&
Foam::rayTracingParticle::trackingData::scattering() const
{
    return scatteringModel_;
}


inline Foam::scalar Foam::rayTracingParticle::trackingData::alphaM(label cellI) const
{
    return alphaM_[cellI];
}


inline Foam::vector Foam::rayTracingParticle::trackingData::gradAlphaM(label cellI) const
{
    return gradAlphaM_[cellI];
}


inline Foam::scalar& Foam::rayTracingParticle::trackingData::Q(label cellI)
{
    return Q_[cellI];
}


inline void Foam::rayTracingParticle::trackingData::markCell(label cellI, cellProperty property)
{
    markedCells_(cellI).set(property);
}


// ************************************************************************* //
