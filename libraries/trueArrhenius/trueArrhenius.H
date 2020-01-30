/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::viscosityModels::trueArrhenius

Description
    trueArrhenius type of dependency on a given scalar field name. Most likely
    temperature. The expression is as follow:
    \verbatim
        mu = exp(alpha_*(inv(T) - inv(Talpha_)))
    \endverbatim

SourceFiles
    trueArrhenius.C

\*---------------------------------------------------------------------------*/

#ifndef trueArrhenius_H
#define trueArrhenius_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class trueArrhenius Declaration
\*---------------------------------------------------------------------------*/

template<class ViscousModel>
class trueArrhenius
:
    public ViscousModel
{
    // Private data

        dictionary trueArrheniusCoeffs_;
        volScalarField temperatureFactor_;

        // Model coefficients
        dimensionedScalar alpha_;
        dimensionedScalar Talpha_;

        //- Name of the field used as temperature
        const word TName_;

        //- Reference to the temperature field
        const volScalarField& T_;


    // Private Member Functions

        //- Calculate the temperature factor for laminar viscosity
        void calcTemperatureFactor();


public:

    //- Runtime type information
    TypeName("trueArrhenius");


    // Constructors

        //- Construct from components
        trueArrhenius
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~trueArrhenius()
    {}


    // Member Functions

        //- Correct the laminar viscosity
        virtual void correct()
        {
            // the temperature factor should be removed for viscous models without correction
            this->nu_ /= temperatureFactor_;
            ViscousModel::correct();
            calcTemperatureFactor();
            this->nu_ *= temperatureFactor_;
        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "trueArrhenius.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
