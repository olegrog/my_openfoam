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
    Foam::polymerComponent

Description
    Mass fraction of a single polymer component.
    Used as part of the multicomponentPolymer for thermal debinding
    simulations.

SourceFiles
    polymerComponent.C

\*---------------------------------------------------------------------------*/

#ifndef polymerComponent_H
#define polymerComponent_H

#include "volFields.H"
#include "dictionaryEntry.H"
#include "constants.H"

namespace Foam
{

class polymerComponent
:
    public volScalarField
{
    // Private data

        word name_;
        dictionary polymerComponentDict_;
        dimensionedScalar rateConstant_;
        dimensionedScalar activationEnergy_;

public:

    // Constructors

        //- Construct from components
        polymerComponent
        (
            const word& name,
            const dictionary& polymerComponentDict,
            const fvMesh& mesh
        );

        //- Return a pointer to a new polymer component created on freestore from Istream
        class iNew
        {
            const fvMesh& mesh_;

        public:

            iNew(const fvMesh& mesh) : mesh_(mesh) {}

            autoPtr<polymerComponent> operator()(Istream& is) const
            {
                dictionaryEntry ent(dictionary::null, is);
                return autoPtr<polymerComponent>
                (
                    new polymerComponent(ent.keyword(), ent, mesh_)
                );
            }
        };

        //- Return clone (need for PtrList)
        autoPtr<polymerComponent> clone() const
        {
            notImplemented("polymerComponent::clone() const");
            return autoPtr<polymerComponent>();
        }

    // Member Functions

        // Access

            const word& name() const
            {
                return name_;
            }

            const word& keyword() const
            {
                return name_;
            }

            //- Return the degradation rate [1/s]
            tmp<volScalarField> degradationRate(const volScalarField& T) const
            {
                using constant::physicoChemical::R;
                return rateConstant_*exp(-activationEnergy_/R/T);
            }
};


} // End namespace Foam

#endif
