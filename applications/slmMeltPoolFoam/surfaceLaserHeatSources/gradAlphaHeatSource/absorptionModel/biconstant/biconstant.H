/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2023 Oleg Rogozin
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
    Foam::absorption::biconstant

Description
    Constant absorptivity but different in liquid and solid.

\*---------------------------------------------------------------------------*/

#ifndef absorption_biconstant_H
#define absorption_biconstant_H

#include "absorptionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace absorption
{

/*---------------------------------------------------------------------------*\
                        Class biconstant Declaration
\*---------------------------------------------------------------------------*/

class biconstant
:
    public absorptionModel
{
    // Private data

        //- Absorptivity coefficient in solid []
        const dimensionedScalar Asolid_;

public:

    //- Runtime type information
    TypeName("biconstant");

    //- Constructor
    biconstant(const dictionary& dict);

    //- Destructor
    ~biconstant() final = default;


    // Member Functions

        //- Return absorptivity []
        tmp<volScalarField> A
        (
            const volVectorField& gradAlphaM,
            const laserProperties& laser
        ) const final;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorption
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
