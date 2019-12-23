/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "alloyComponent.H"

namespace Foam {

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alloyComponent::alloyComponent(
    const word& name,
    const dictionary& alloyComponentDict,
    const fvMesh& mesh,
    const dimensionedScalar& meltingTemp
) :
    volScalarField(
        IOobject(
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
    densityMelting_("densityMelting", alloyComponentDict_),
    molarMass_("molarMass", alloyComponentDict_),
    equilibriumS_("equilibriumS", alloyComponentDict_),
    equilibriumL_("equilibriumL", alloyComponentDict_),
    slopeS_("slopeS", alloyComponentDict_),
    slopeL_("slopeL", alloyComponentDict_),
    diffusionS_("diffusionS", alloyComponentDict_),
    diffusionL_("diffusionL", alloyComponentDict_),
    meltingTemp_(meltingTemp)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

autoPtr<alloyComponent> alloyComponent::clone() const
{
    notImplemented("alloyComponent::clone() const");
    return autoPtr<alloyComponent>();
}

} // End namespace Foam
