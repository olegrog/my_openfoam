/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2021 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of laserHeatSource.

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

Class
    Foam::beam::Gaussian

Description
    Abstract base class for laser beam properties.

SourceFiles
    Gaussian.C

\*---------------------------------------------------------------------------*/

#ifndef beam_Gaussian_H
#define beam_Gaussian_H

#include "laserBeam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace beam
{

/*---------------------------------------------------------------------------*\
                          Class Gaussian Declaration
\*---------------------------------------------------------------------------*/

class Gaussian
:
    public laserBeam
{

public:

    //- Runtime type information
    TypeName("Gaussian");

    //- Constructor
    Gaussian
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const laserProperties& laser
    )
    :
        laserBeam(dict, mesh, laser) {}

    //- Destructor
    ~Gaussian() override = default;


    // Member Functions

        //- Return the intensity as a function of coordinate [W/m2]
        scalar I(const vector& r) const override;

        //- Return the intensity distribution [W/m2]
        tmp<volScalarField> I() const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace beam
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
