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
    Foam::absorptivity::constant

Description
    Constant absorptivity.

\*---------------------------------------------------------------------------*/

#ifndef absorptivity_constant_H
#define absorptivity_constant_H

#include "absorptivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace absorptivity
{

/*---------------------------------------------------------------------------*\
                        Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public absorptivityModel
{

public:

    //- Runtime type information
    TypeName("constant");

    //- Constructor
    constant(const dictionary& dict) : absorptivityModel(dict) {}

    //- Destructor
    virtual ~constant() = default;

    // Member Functions

        //- Return the absorptivity field
        tmp<volScalarField> value
        (
            const volScalarField& alphaM,
            const volVectorField& gradAlphaM,
            const interfacialLaserHeatSource& laser
        ) const final
        {
            return A_*mag(gradAlphaM);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptivity
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
