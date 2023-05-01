/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020-2023 Oleg Rogozin
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

#include "gasMetalThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gasMetalThermo::gasMetalThermo(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "thermalProperties",
            mesh.time().constant(),
            mesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    metalDict_(subDict("metal")),
    solidDict_(metalDict_.subDict("solid")),
    liquidDict_(metalDict_.subDict("liquid")),
    gasDict_(subDict("gas")),
    solid_(solidDict_),
    liquid_(liquidDict_),
    gas_(gasDict_),
    Tmelting_(metalDict_.get<scalar>("Tmelting")),
    Tboiling_(metalDict_.get<scalar>("Tboiling")),
    Hfusion_(metalDict_.get<scalar>("Hfusion")),
    Hvapour_(metalDict_.get<scalar>("Hvapour")),
    metalW_(metalDict_.get<scalar>("molWeight")),
    gasW_(gasDict_.get<scalar>("molWeight")),
    sigmoidPtr_(sigmoidFunction::New(*this, 0, 1))
{
    scalar hSolidus = solid_.Cp.integral(0, Tmelting_);
    scalar hLiquidus = hSolidus + Hfusion_;
    scalar hGasAtMelting = gas_.Cp.integral(0, Tmelting_);
    scalar hPreBoiling = hLiquidus + liquid_.Cp.integral(Tmelting_, Tboiling_);

    Info<< "\nThermal properties:" << endl
        << " -- Solidus enthalpy = " << hSolidus << endl
        << " -- Liquidus enthalpy = " << hLiquidus << endl
        << " -- Gas enthalpy at Tmelting = " << hGasAtMelting << endl
        << " -- Pre-boiling enthalpy = " << hPreBoiling << endl
        << " -- Solidus heat capacity = " << solid_.Cp.value(Tmelting_) << endl
        << " -- Liquidus heat capacity = " << liquid_.Cp.value(Tmelting_) << endl
        << " -- Solidus thermal conductivity = " << solid_.kappa.value(Tmelting_) << endl
        << " -- Liquidus thermal conductivity = " << liquid_.kappa.value(Tmelting_) << endl
        << " -- Gas thermal conductivity at Tmelting = " << gas_.kappa.value(Tmelting_) << endl
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::gasMetalThermo::hAtMeltingPrime() const
{
    return dimensionedScalar
    (
        "hAtMeltingPrime",
        dimEnergy/dimMass,
        gas_.Cp.integral(0, Tmelting_) - solid_.Cp.integral(0, Tmelting_) - Hfusion_/2
    );
}


// ************************************************************************* //
