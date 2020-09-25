/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020-2020 Oleg Rogozin
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

#include "multicomponentPolymer.H"

#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multicomponentPolymer::multicomponentPolymer(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "polymerProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    initialVolumeFraction_("initialVolumeFraction", dimless, *this),
    rhoPolymer_("rhoPolymer", dimDensity, *this),
    diffusionConstant_("diffusionConstant", dimViscosity, *this),
    activationEnergy_("activationEnergy", dimEnergy/dimMoles, *this),
    components_(lookup("components"), polymerComponent::iNew(mesh))
{}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::multicomponentPolymer::diffusion
(
    const volScalarField& rho,
    const volScalarField& T
) const
{
    using constant::physicoChemical::R;
    return diffusionConstant_*exp(-activationEnergy_/R/T);
}


// ************************************************************************* //
