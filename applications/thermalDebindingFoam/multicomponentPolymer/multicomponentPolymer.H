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
    Foam::multicomponentPolymer

Description
    System of a multicomponent polymer with its monomer.
    Volume of the system is constant.

SourceFiles
    multicomponentPolymer.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentPolymer_H
#define multicomponentPolymer_H

#include "volFields.H"
#include "PtrDictionary.H"

#include "polymerComponent.H"
#include "diffusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class multicomponentPolymer Declaration
\*---------------------------------------------------------------------------*/

class multicomponentPolymer
:
    public IOdictionary
{
    // Private data

        //- Density of the polymer [kg/m^3]
        const dimensionedScalar polymerRho_;

        //- Molar mass of the monomer [kg/mol]
        const dimensionedScalar monomerW_;

        //- Total volume fraction of polymer and pores []
        const dimensionedScalar totalVolumeFraction_;

        //- Initial porosity []
        const dimensionedScalar initialPorosity_;

        //- The diffusion model
        autoPtr<diffusionModel> diffusionModelPtr_;

        //- Pointer to the dictionary of polymer components
        PtrDictionary<polymerComponent> components_;

public:

    //- Constructor
    multicomponentPolymer(const fvMesh& mesh);

    //- Destructor
    ~multicomponentPolymer() = default;

    // Member Functions

        //- Return the polymer components
        PtrDictionary<polymerComponent>& components()
        {
            return components_;
        }

        // Access

            //- Return the polymer density [kg/m^3]
            dimensionedScalar rho() const
            {
                return polymerRho_;
            }

            //- Return the initial volume fraction of polymer
            dimensionedScalar initialVolumeFraction() const
            {
                return totalVolumeFraction_ - initialPorosity_;
            }

            //- Return the volume fraction of polymer
            tmp<volScalarField> volumeFraction() const;

            //- Return the volume fraction of pores
            tmp<volScalarField> poresFraction() const
            {
                return totalVolumeFraction_ - volumeFraction();
            }

            //- Return the total volume fraction
            dimensionedScalar totalVolumeFraction() const
            {
                return totalVolumeFraction_;
            }


        // Evaluation

            //- Return the monomer pressure [Pa]
            tmp<volScalarField> pressure
            (
                const volScalarField& rho,
                const volScalarField& T
            ) const;

            //- Return the diffusion coefficient [m^2/s]
            tmp<volScalarField> diffusion
            (
                const volScalarField& rho,
                const volScalarField& T
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
