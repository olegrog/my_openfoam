/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of trueArrhenius.

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

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Flag to indicate whether the viscosity field is multiplied by a temperature factor
        bool isCorrected_;


    // Private Member Functions

        //- Calculate the temperature factor for laminar viscosity
        void calcTemperatureFactor(const volScalarField& T);


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
    ~trueArrhenius() override;


    // Member Functions

        //- Correct the laminar viscosity
        void correct() override;

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties) override;
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
