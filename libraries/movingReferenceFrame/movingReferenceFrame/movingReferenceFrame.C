/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2022 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of movingReferenceFrame.

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

#include "movingReferenceFrame.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(movingReferenceFrame);
    defineRunTimeSelectionTable(movingReferenceFrame, mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::movingReferenceFrame> Foam::movingReferenceFrame::New
(
    const fvMesh& mesh
)
{
    const IOdictionary dict
    (
        IOobject
        (
            "movingFrameDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    const word frameType(dict.get<word>("type"));

    Info<< "Selecting moving reference frame type " << frameType << endl;

    const auto cstrIter = meshConstructorTablePtr_->cfind(frameType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "movingReferenceFrame",
            frameType,
            *meshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<movingReferenceFrame>(cstrIter()(mesh));
}


Foam::movingReferenceFrame::movingReferenceFrame(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "movingFrameDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    UrelDict_
    (
        IOobject
        (
            "Urel",
            mesh.time().timeName(),
            "uniform",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    Urel_("Urel", dimVelocity, vector::zero, UrelDict_),
    phiRel_("phiRel", Urel_ & mesh.Sf())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingReferenceFrame::correct(bool oriented)
{
    evaluate();

    if (mag(Urel_.value()) > 0)
    {
        UrelDict_.set(Urel_.name(), Urel_.value());
        Info<< "Urel = " << Urel_.value() << endl;
    }

    phiRel_ = Urel_ & mesh_.Sf();
    phiRel_.setOriented(oriented);
}


// ************************************************************************* //
