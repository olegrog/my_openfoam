/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                  Arrhenius | Copyright (C) 2017 OpenFOAM Foundation
              trueArrhenius | Copyright (C) 2019, 2022 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of trueArrhenius.

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

#include "updateGeometricField.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ViscosityModel>
Foam::scalar Foam::viscosityModels::trueArrhenius<ViscosityModel>::temperatureFactor
(
    scalar T
) const
{
    return min(exp(alpha_*(inv(T) - inv(Talpha_))), 1);
}


template<class ViscosityModel>
void Foam::viscosityModels::trueArrhenius<ViscosityModel>::calcTemperatureFactor
(
    const volScalarField& T
)
{
    updateGeometricField
    (
        temperatureFactor_,
        [this](scalar& factor, scalar T)
        {
            factor = temperatureFactor(T);
        },
        T
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ViscosityModel>
Foam::viscosityModels::trueArrhenius<ViscosityModel>::trueArrhenius
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    ViscosityModel(name, viscosityProperties, U, phi),
    trueArrheniusCoeffs_
    (
        viscosityProperties.subDict(typeName + "Coeffs")
    ),
    alpha_(trueArrheniusCoeffs_.get<scalar>("alpha")),
    Talpha_(trueArrheniusCoeffs_.get<scalar>("Talpha")),
    TName_(trueArrheniusCoeffs_.lookupOrDefault<word>("field", "T")),
    mesh_(U.mesh()),
    temperatureFactor_
    (
        IOobject
        (
            name + "TemperatureFactor",
            U.time().timeName(),
            U.db()
        ),
        mesh_,
        dimless
    ),
    isCorrected_(false)
{
    const scalar nu0 = gMax(ViscosityModel::nu()());

    Info<< " -- kinematic viscosity at Talpha = " << nu0 << endl
        << " -- kinematic viscosity at Tinfty = " << nu0*temperatureFactor(GREAT*Talpha_) << endl
        << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ViscosityModel>
void Foam::viscosityModels::trueArrhenius<ViscosityModel>::correct()
{
    const auto* TPtr = mesh_.findObject<volScalarField>(TName_);

    // The temperature field is not available in the first call
    if (TPtr)
    {
        if (isCorrected_)
        {
            // The temperature factor should be removed for viscous models without correction
            this->nu_ /= temperatureFactor_;
        }

        ViscosityModel::correct();

        calcTemperatureFactor(*TPtr);
        this->nu_ *= temperatureFactor_;
        isCorrected_ = true;
    }
}


template<class ViscosityModel>
bool Foam::viscosityModels::trueArrhenius<ViscosityModel>::read
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
