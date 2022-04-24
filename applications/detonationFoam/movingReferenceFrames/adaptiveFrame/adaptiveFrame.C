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
    This file is part of detonationFoam.

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
    meanValue_(dict_.get<scalar>("meanValue")),
    valueFactor_(dict_.get<scalar>("valueFactor")),
    derivativeFactor_(dict_.get<scalar>("derivativeFactor")),
    prevFieldPtr_(nullptr),
    prevUrel_(vector::zero)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::adaptiveFrame::evaluate()
{
    const scalar dt = mesh_.time().deltaTValue();
    const scalar meshVolume = gSum(mesh_.V());
    const volScalarField& field = mesh_.lookupObject<volScalarField>(fieldName_);
    const scalar meanValueCurr = fvc::domainIntegrate(field).value()/meshVolume;
    const scalar delta = meanValueCurr - meanValue_;

    // Estimate the time derivative
    scalar dotMean = 0;
    if (prevFieldPtr_)
    {
        const scalar dt0 = mesh_.time().deltaT0Value();
        dotMean = fvc::domainIntegrate(field - *prevFieldPtr_).value()/meshVolume/dt0;
        *prevFieldPtr_ == field;
    }
    else
    {
        prevFieldPtr_.reset(new volScalarField(field.name() + "Prev", field));
    }

    if (debug)
    {
        Info<< "delta = " << delta << " dotMean = " << dotMean << endl;
        Info<< "meanValue = " << meanValueCurr << endl;
    }

    const scalar dUrelDt = valueFactor_*delta + derivativeFactor_*dotMean;
    prevUrel_ = Urel_.value();
    Urel_.value() = prevUrel_ + direction_*dUrelDt*dt;
}


// ************************************************************************* //
