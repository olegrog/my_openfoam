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
    Foam::diffusion::freeVolume

Description
    Diffusion coefficient based on the free-volume (Vrentas--Duda) theory.
    See Duda, Vrentas, Ju, and Liu. AIChE Journal 28, no. 2 (1982): 279-285.

SourceFiles
    freeVolume.C

\*---------------------------------------------------------------------------*/

#ifndef diffusion_freeVolume_H
#define diffusion_freeVolume_H

#include "diffusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusion
{

/*---------------------------------------------------------------------------*\
                         Class freeVolume Declaration
\*---------------------------------------------------------------------------*/

class freeVolume
:
    public diffusionModel
{

    struct TVK
    {
        //- Glass transition temperature
        const dimensionedScalar Tg;

        //- Speficic volume at zero temperature
        const dimensionedScalar V;

        //- Free-volume parameter
        const dimensionedScalar K;

        //- Constructor
        TVK(const dictionary& dict)
        :
            Tg("Tg", dimTemperature, dict),
            V("V", inv(dimDensity), dict),
            K("K", inv(dimDensity*dimTemperature), dict)
        {}
    };

    // Private data

        //- Coefficients dictionary
        const dictionary coeffsDict_;

        //- Properties of all phases
        const TVK monomer_, polymer_;

        //- Interaction parameter of the Flory--Huggins theory
        const scalar chi_;

        //- Ratio of critical molar volumes of hopping units
        const scalar xi_;

public:

    //- Runtime type information
    TypeName("freeVolume");

    //- Constructor
    freeVolume(const dictionary& dict);

    //- Destructor
    virtual ~freeVolume() = default;

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
