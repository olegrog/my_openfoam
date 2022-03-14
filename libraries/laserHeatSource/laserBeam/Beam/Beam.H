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
    Foam::Beam

Description
    Class template for arbitrary laser beam intensity.

SourceFiles
    BeamTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Beam_H
#define Beam_H

#include "laserBeam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class Beam Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Beam
:
    public laserBeam
{

public:

    //- Runtime type information
    TypeName("Beam");

    //- Constructor
    Beam
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const laserProperties& laser
    );

    //- Destructor
    ~Beam() final = default;


    // Member Functions

        //- Intensity as a function of radial coordinate [W/m2]
        scalar I(scalar r) const final;

        //- Integral of the intensity over the disk of radius r [W]
        scalar integralI(scalar r) const final;

        //- Intensity field [W/m2]
        tmp<volScalarField> I() const final;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BeamTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
