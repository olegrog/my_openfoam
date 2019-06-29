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
    densityMelting_(alloyComponentDict_.lookup("densityMelting")),
    molarMass_(alloyComponentDict_.lookup("molarMass")),
    equilibriumS_(alloyComponentDict_.lookup("equilibriumS")),
    equilibriumL_(alloyComponentDict_.lookup("equilibriumL")),
    slopeS_(alloyComponentDict_.lookup("slopeS")),
    slopeL_(alloyComponentDict_.lookup("slopeL")),
    diffusionS_(alloyComponentDict_.lookup("diffusionS")),
    diffusionL_(alloyComponentDict_.lookup("diffusionL")),
    meltingTemp_(meltingTemp)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

autoPtr<alloyComponent> alloyComponent::clone() const
{
    notImplemented("alloyComponent::clone() const");
    return autoPtr<alloyComponent>(NULL);
}

} // End namespace Foam
