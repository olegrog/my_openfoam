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
    Foam::absorption::constant

Description
    Constant absorptivity.

\*---------------------------------------------------------------------------*/

#ifndef absorption_constant_H
#define absorption_constant_H

#include "absorptionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace absorption
{

/*---------------------------------------------------------------------------*\
                        Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public absorptionModel
{

public:

    //- Runtime type information
    TypeName("constant");

    //- Constructor
    constant(const dictionary& dict) : absorptionModel(dict) {}

    //- Destructor
    ~constant() final = default;


    // Member Functions

        //- Return absorptivity []
        tmp<volScalarField> A
        (
            const volVectorField& gradAlphaM,
            const laserProperties& laser
        ) const final
        {
            return A_*mag(gradAlphaM);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorption
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
