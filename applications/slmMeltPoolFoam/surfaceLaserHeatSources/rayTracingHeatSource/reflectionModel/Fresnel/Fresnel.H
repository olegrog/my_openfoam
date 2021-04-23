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
    Foam::reflection::Fresnel

Description
    General Fresnel reflection model between gas and absorbing media.

    J.A. Stratton. Electromagnetic theory. 1941. Chapter 9.

SourceFiles
    Fresnel.C

\*---------------------------------------------------------------------------*/

#ifndef reflection_Fresnel_H
#define reflection_Fresnel_H

#include "reflectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reflection
{

/*---------------------------------------------------------------------------*\
                           Class Fresnel Declaration
\*---------------------------------------------------------------------------*/

class Fresnel
:
    public reflectionModel
{
    // Private data

        //- Refractive index of metal
        const scalar n_;

        //- Extinction coefficient of metal
        const scalar k_;

public:

    //- Runtime type information
    TypeName("Fresnel");

    //- Constructor
    Fresnel(const dictionary& dict);

    //- Destructor
    ~Fresnel() final = default;


    // Member Functions

        //- Return reflected direction
        vector R(const vector& incident, const vector& n) const final;

        //- Return reflectivity
        scalar rho(const scalar incidentAngle) const final;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reflection
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
