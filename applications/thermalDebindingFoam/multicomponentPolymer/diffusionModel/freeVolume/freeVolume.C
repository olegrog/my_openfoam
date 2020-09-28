/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of thermalDebindingFoam.

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

#include "freeVolume.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace diffusion
    {
        defineTypeName(freeVolume);
        addToRunTimeSelectionTable(diffusionModel, freeVolume, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusion::freeVolume::freeVolume
(
    const dictionary& dict
)
:
    diffusionModel(dict),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    monomer_(coeffsDict_.subDict("monomer")),
    polymer_(coeffsDict_.subDict("polymer")),
    chi_(coeffsDict_.get<scalar>("chi")),
    xi_(coeffsDict_.get<scalar>("xi"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diffusion::freeVolume::D
(
    const volScalarField& T,
    const volScalarField& monomerVolumeFraction,
    const volScalarField& monomerMassFraction
) const
{
    volScalarField polymerVolumeFraction = 1 - monomerVolumeFraction;
    volScalarField polymerMassFraction = 1 - monomerMassFraction;

    volScalarField V = monomer_.K*(T - monomer_.Tg)*monomerMassFraction
        + polymer_.K*(T - polymer_.Tg)*polymerMassFraction;
    V.max(dimensionedScalar(inv(dimDensity), SMALL));

    return diffusionModel::D(T, monomerVolumeFraction, monomerMassFraction)
        *sqr(polymerVolumeFraction)*(1 - 2*chi_*monomerVolumeFraction)
        *Foam::exp(-(monomerMassFraction*monomer_.V + xi_*polymerMassFraction*polymer_.V)/V);
}


// ************************************************************************* //
