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
    This file is part of beamFunction.

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

#include "addToRunTimeSelectionTable.H"

#include "Beam.H"

#include "Gaussian.H"
#include "topHat.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeLaserBeam(beamType, type)                                          \
                                                                               \
namespace Foam                                                                 \
{                                                                              \
    typedef beamType<beam::type> beamType##type;                               \
                                                                               \
    addNamedToRunTimeSelectionTable                                            \
    (                                                                          \
        laserBeam,                                                             \
        beamType##type,                                                        \
        dictionary,                                                            \
        type                                                                   \
    );                                                                         \
                                                                               \
    defineTemplateTypeNameWithName                                             \
    (                                                                          \
        beamType##type,                                                        \
        #beamType"<"#type">"                                                   \
    );                                                                         \
}

makeLaserBeam(Beam, Gaussian)
makeLaserBeam(Beam, topHat)

// ************************************************************************* //
