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

#include "adaptiveFrame.H"

#include "addToRunTimeSelectionTable.H"
#include "fvcVolumeIntegrate.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adaptiveFrame, 0);
    addToRunTimeSelectionTable(movingReferenceFrame, adaptiveFrame, mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adaptiveFrame::adaptiveFrame(const fvMesh& mesh)
:
    movingReferenceFrame(mesh),
    dict_(subDict(typeName + "Frame")),
    direction_(dict_.get<vector>("direction")),
    fieldName_(dict_.get<word>("fieldName")),
    field_(mesh_.lookupObject<volScalarField>(fieldName_)),
    meanValue_(dict_.get<scalar>("meanValue")),
    meanValuePrev_
    (
        UrelDict_.getOrDefault
        (
            "meanValue",
            fvc::domainIntegrate(field_).value()/gSum(mesh.V())
        )
    ),
    valueFactor_(dict_.get<scalar>("valueFactor")),
    derivativeFactor_(dict_.get<scalar>("derivativeFactor"))
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::adaptiveFrame::evaluate()
{
    const scalar dt = mesh_.time().deltaTValue();
    const scalar dt0 = mesh_.time().deltaT0Value();
    const scalar meanValueCurr = fvc::domainIntegrate(field_).value()/gSum(mesh_.V());
    const scalar delta = meanValueCurr - meanValue_;
    const scalar dotMean = (meanValueCurr - meanValuePrev_)/dt0;

    meanValuePrev_ = meanValueCurr;
    UrelDict_.set("meanValue", meanValuePrev_);

    if (debug)
    {
        Info<< "delta = " << delta << " dotMean = " << dotMean << endl;
        Info<< "meanValue = " << meanValueCurr << endl;
    }

    const scalar dUrelDt = valueFactor_*delta + derivativeFactor_*dotMean;
    Urel_.value() += direction_*dUrelDt*dt;
}


// ************************************************************************* //
