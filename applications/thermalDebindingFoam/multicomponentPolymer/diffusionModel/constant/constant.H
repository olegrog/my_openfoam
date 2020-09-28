/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of thermalDebindingFoam.

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
    Foam::diffusion::constant

Description
    Diffusion coefficient independent of polymer volume fraction.

SourceFiles
    constant.C

\*---------------------------------------------------------------------------*/

#ifndef diffusion_constant_H
#define diffusion_constant_H

#include "diffusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusion
{

/*---------------------------------------------------------------------------*\
                         Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public diffusionModel
{

public:

    //- Runtime type information
    TypeName("constant");

    //- Constructor
    constant(const dictionary& dict) : diffusionModel(dict) {}

    //- Destructor
    virtual ~constant() = default;

    // Member Functions

        //- Return the value of the diffusion coefficient
        tmp<volScalarField> D
        (
            const volScalarField& T,
            const volScalarField& monomerVolumeFraction,
            const volScalarField& monomerMassFraction
        ) const final;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusion
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
