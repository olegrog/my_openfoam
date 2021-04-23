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
    Foam::absorptionModel

Description
    Abstract base class for absorption model in surface heat sources.

\*---------------------------------------------------------------------------*/

#ifndef absorptionModel_H
#define absorptionModel_H

#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class laserProperties;

/*---------------------------------------------------------------------------*\
                        Class absorptionModel Declaration
\*---------------------------------------------------------------------------*/

class absorptionModel
{

protected:

    // Protected data

        //- Absorptivity coefficient []
        const dimensionedScalar A_;

public:

    //- Runtime type information
    TypeName("absorptionModel");

    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        absorptionModel,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    //- Constructor
    absorptionModel(const dictionary& dict);

    //- Selector
    static autoPtr<absorptionModel> New(const dictionary& dict);

    //- Destructor
    virtual ~absorptionModel() = default;


    // Member Functions

        //- Return absorptivity []
        virtual tmp<volScalarField> A
        (
            const volVectorField& gradAlphaM,
            const laserProperties& laser
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
