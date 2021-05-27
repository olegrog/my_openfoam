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
    This file is part of slmMeltPoolFoam.

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
    Foam::scattering::FresnelNoRefraction

Description
    General Fresnel's relations without refraction.

SourceFiles
    FresnelNoRefraction.C

\*---------------------------------------------------------------------------*/

#ifndef scattering_FresnelNoRefraction_H
#define scattering_FresnelNoRefraction_H

#include "Fresnel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace scattering
{

/*---------------------------------------------------------------------------*\
                    Class FresnelNoRefraction Declaration
\*---------------------------------------------------------------------------*/

class FresnelNoRefraction
:
    public Fresnel
{
public:

    //- Runtime type information
    TypeName("FresnelNoRefraction");

    //- Constructor
    FresnelNoRefraction(const dictionary& dict) : Fresnel(dict) {}

    //- Destructor
    ~FresnelNoRefraction() final = default;


    // Member Functions

        //- Return direction of refraction
        vector refraction(const vector& i, const vector& n) const final;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace scattering
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
