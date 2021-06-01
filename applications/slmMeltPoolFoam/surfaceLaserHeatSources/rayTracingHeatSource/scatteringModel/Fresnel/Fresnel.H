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
    Foam::scattering::Fresnel

Description
    General Fresnel's relations for vacuum and absorbing media.

    J.A. Stratton. Electromagnetic theory. 1941. Chapter 9.
    M. Modest. Radiative heat transfer. 2013. Chapter 2.5.

SourceFiles
    Fresnel.C

\*---------------------------------------------------------------------------*/

#ifndef scattering_Fresnel_H
#define scattering_Fresnel_H

#include "scatteringModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace scattering
{

/*---------------------------------------------------------------------------*\
                           Class Fresnel Declaration
\*---------------------------------------------------------------------------*/

class Fresnel
:
    public scatteringModel
{
    // Private data

        //- Refractive index of metal
        const scalar n_;

        //- Absorptive index of metal
        const scalar k_;

public:

    //- Runtime type information
    TypeName("Fresnel");

    //- Constructor
    Fresnel(const dictionary& dict);

    //- Destructor
    ~Fresnel() override = default;


    // Member Functions

        //- Return direction of reflection
        vector reflection(const vector& i, const vector& n) const override;

        //- Return direction of refraction
        vector refraction(const vector& i, const vector& n) const override;

        //- Return reflectivity
        scalar R(scalar cosTheta) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace scattering
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
