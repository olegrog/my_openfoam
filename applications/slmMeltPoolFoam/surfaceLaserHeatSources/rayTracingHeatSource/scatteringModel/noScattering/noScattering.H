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
    Foam::scattering::noScattering

Description
    Dummy scattering model for 'none'

SourceFiles
    noScattering.C

\*---------------------------------------------------------------------------*/

#ifndef scattering_noScattering_H
#define scattering_noScattering_H

#include "scatteringModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace scattering
{

/*---------------------------------------------------------------------------*\
                         Class noScattering Declaration
\*---------------------------------------------------------------------------*/

class noScattering
:
    public scatteringModel
{

public:

    //- Runtime type information
    TypeName("none");

    //- Constructor
    noScattering(const dictionary& dict) : scatteringModel(dict) {}

    //- Destructor
    ~noScattering() final = default;


    // Member Functions

        //- Return direction of reflection
        vector reflection(const vector& i, const vector& n) const final
        {
            return vector::zero;
        }

        //- Return direction of refraction
        vector refraction(const vector& i, const vector& n) const final
        {
            return i;
        }

        //- Return reflectivity
        scalar R(scalar cosTheta) const final
        {
            return 0;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace scattering
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
