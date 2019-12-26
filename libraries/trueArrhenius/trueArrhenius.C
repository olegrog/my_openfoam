/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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


#include "trueArrhenius.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ViscousModel>
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::trueArrhenius<ViscousModel>::calcNu
(
    const volScalarField& field
) const
{
    return exp(-alpha_*(1/field - 1/Talpha_));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ViscousModel>
Foam::viscosityModels::trueArrhenius<ViscousModel>::trueArrhenius
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    ViscousModel(name, viscosityProperties, U, phi),
    trueArrheniusCoeffs_
    (
        viscosityProperties.subDict(typeName + "Coeffs")
    ),
    alpha_("alpha", dimTemperature, trueArrheniusCoeffs_),
    Talpha_("Talpha", dimTemperature, trueArrheniusCoeffs_),
    fieldName_(trueArrheniusCoeffs_.lookupOrDefault<word>("field", "T")),
    mesh_(U.mesh())
{
    const auto* fldPtr = mesh_.findObject<volScalarField>(fieldName_);

    if (fldPtr)
    {
        this->nu_ *= calcNu(*fldPtr);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ViscousModel>
bool Foam::viscosityModels::trueArrhenius<ViscousModel>::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    trueArrheniusCoeffs_ =
        viscosityProperties.subDict(typeName + "Coeffs");

    trueArrheniusCoeffs_.readEntry("alpha", alpha_);
    trueArrheniusCoeffs_.readEntry("Talpha", Talpha_);

    return true;
}


// ************************************************************************* //
