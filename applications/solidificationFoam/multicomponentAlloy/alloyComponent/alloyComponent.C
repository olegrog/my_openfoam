/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019-2020 Oleg Rogozin
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

#include "alloyComponent.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::alloyComponent::alloyComponent
(
    const word& name,
    const dictionary& alloyComponentDict,
    const fvMesh& mesh,
    const dimensionedScalar& Tmelting
)
:
    volScalarField
    (
        IOobject
        (
            "concentration" + name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.lookupObject<volScalarField>("phase")
    ),
    name_(name),
    alloyComponentDict_(alloyComponentDict),
    molarMass_("molarMass", dimMass/dimMoles, alloyComponentDict_),
    equilibriumS_("equilibriumS", dimless, alloyComponentDict_),
    equilibriumL_("equilibriumL", dimless, alloyComponentDict_),
    slopeS_("slopeS", dimTemperature, alloyComponentDict_),
    slopeL_("slopeL", dimTemperature, alloyComponentDict_),
    diffusionS_("diffusionS", dimArea/dimTime, alloyComponentDict_),
    diffusionL_("diffusionL", dimArea/dimTime, alloyComponentDict_),
    Tmelting_(Tmelting)
{}


// ************************************************************************* //
