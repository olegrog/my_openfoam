/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2020-2020 Oleg Rogozin
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
    Multicomponent polymer.

SourceFiles
    multicomponentPolymer.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentPolymer_H
#define multicomponentPolymer_H

#include "volFields.H"
#include "PtrDictionary.H"

#include "polymerComponent.H"

namespace Foam
{

class multicomponentPolymer
:
    public IOdictionary
{
    // Private data

        //- Initial volume fraction of the polymer []
        const dimensionedScalar initialVolumeFraction_;

        //- Density of the polymer [kg/m^3]
        const dimensionedScalar rhoPolymer_;

        //- Diffusion constant [m^2/s]
        const dimensionedScalar diffusionConstant_;

        //- Activation energy for diffusion [J/mol]
        const dimensionedScalar activationEnergy_;

        //- Pointer to the dictionary of polymer components
        PtrDictionary<polymerComponent> components_;

    // Private Member functions

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

            //- Return the initial volume fraction []
            dimensionedScalar initialVolumeFraction() const
            {
                return initialVolumeFraction_;
            }

            //- Return the polymer density [kg/m^3]
            dimensionedScalar rho() const
            {
                return rhoPolymer_;
            }

            //- Return the diffusion coefficient [m^2/s]
            tmp<volScalarField> diffusion
            (
                const volScalarField& rho,
                const volScalarField& T
            ) const;
};


} // End namespace Foam

#endif
