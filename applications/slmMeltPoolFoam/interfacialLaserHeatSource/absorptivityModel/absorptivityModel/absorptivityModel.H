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
    Foam::absorptivityModel

Description
    Abstract base class for absorptivity model for interfacial heat source.

\*---------------------------------------------------------------------------*/

#ifndef absorptivityModel_H
#define absorptivityModel_H

#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class interfacialLaserHeatSource;

/*---------------------------------------------------------------------------*\
                        Class absorptivityModel Declaration
\*---------------------------------------------------------------------------*/

class absorptivityModel
{

protected:

    // Protected data

        //- Absorptivity coefficient []
        const dimensionedScalar A_;

public:

    //- Runtime type information
    TypeName("absorptivityModel");

    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        absorptivityModel,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    //- Constructor
    absorptivityModel(const dictionary& dict);

    //- Selector
    static autoPtr<absorptivityModel> New(const dictionary& dict);

    //- Destructor
    virtual ~absorptivityModel() = default;

    // Member Functions

        //- Return the absorptivity field
        virtual tmp<volScalarField> value
        (
            const volScalarField& alphaM,
            const volVectorField& gradAlphaM,
            const interfacialLaserHeatSource& laser
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
