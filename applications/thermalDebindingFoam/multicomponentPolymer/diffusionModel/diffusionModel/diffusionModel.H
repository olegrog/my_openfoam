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
    Foam::diffusionModel

Description
    Abstract base class for diffusion model of monomer--polymer mixture.

\*---------------------------------------------------------------------------*/

#ifndef diffusionModel_H
#define diffusionModel_H

#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class diffusionModel Declaration
\*---------------------------------------------------------------------------*/

class diffusionModel
{
    // Private data

        //- Diffusion constant [m^2/s]
        const dimensionedScalar diffusionConstant_;

        //- Activation energy for diffusion [J/mol]
        const dimensionedScalar activationEnergy_;

public:

    //- Runtime type information
    TypeName("diffusionModel");

    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        diffusionModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    //- Constructor
    diffusionModel(const dictionary& dict);

    //- Selector
    static autoPtr<diffusionModel> New(const dictionary& dict);

    //- Destructor
    virtual ~diffusionModel() = default;

    // Member Functions

        //- Return the value of the diffusion coefficient
        virtual tmp<volScalarField> D
        (
            const volScalarField& T,
            const volScalarField& monomerVolumeFraction,
            const volScalarField& monomerMassFraction
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
