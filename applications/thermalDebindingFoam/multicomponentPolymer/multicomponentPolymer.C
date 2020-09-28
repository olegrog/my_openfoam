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
    polymerRho_("polymerRho", dimDensity, *this),
    monomerW_("monomerW", dimMass/dimMoles, *this),
    totalVolumeFraction_("totalVolumeFraction", dimless, *this),
    initialPorosity_("initialPorosity", dimless, *this),
    diffusionModelPtr_(diffusionModel::New(subDict("diffusionModel"))),
    components_(lookup("components"), polymerComponent::iNew(mesh))
{
    // --- Checks

    volScalarField total = volumeFraction()/(totalVolumeFraction_ - initialPorosity_);

    if (gMax(total) > 1)
    {
        FatalError
            << "The total sum of components is more than 1."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::multicomponentPolymer::volumeFraction() const
{
    auto iter = components_.begin();

    tmp<volScalarField> result = iter();

    for (++iter; iter != components_.end(); ++iter)
    {
        result = result() + iter();
    }

    return result()*(totalVolumeFraction_ - initialPorosity_);
}


Foam::tmp<Foam::volScalarField> Foam::multicomponentPolymer::pressure
(
    const volScalarField& rho,
    const volScalarField& T
) const
{
    using constant::physicoChemical::R;

    dimensionedScalar small(inv(dimless), SMALL);
    return rho*R*T/monomerW_/(poresFraction() + small);
}


Foam::tmp<Foam::volScalarField> Foam::multicomponentPolymer::diffusion
(
    const volScalarField& rho,
    const volScalarField& T
) const
{
    // Fractions in the polymer--monomer system
    const volScalarField monomerVolumeFraction = 1 - volumeFraction()/totalVolumeFraction_;
    const volScalarField monomerMassFraction = rho/(rho + volumeFraction()*polymerRho_);

    return diffusionModelPtr_->D(T, monomerVolumeFraction, monomerMassFraction)
        *Foam::pow(totalVolumeFraction_, 1.5);
}


// ************************************************************************* //
