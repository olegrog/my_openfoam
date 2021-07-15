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
    This file is part of solidificationFoam.

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
    Foam::componentPhase

Description
    Abstract base class for an alloy component phase.

\*---------------------------------------------------------------------------*/

#ifndef componentPhase_H
#define componentPhase_H

#include "volFields.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class alloyComponent;

/*---------------------------------------------------------------------------*\
                           Class componentPhase Declaration
\*---------------------------------------------------------------------------*/

class componentPhase
{
    // Private Data

        //- Name of the phase
        const word name_;

        //- Diffusion coefficient [m2/s]
        const dimensionedScalar diffusion_;

        //- Reference to the alloy component
        const alloyComponent& component_;

protected:

    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

public:

    //- Runtime type information
    TypeName("componentPhase");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        componentPhase,
        dictionary,
        (
            const dictionaryEntry& entry,
            const fvMesh& mesh,
            const alloyComponent& component
        ),
        (entry, mesh, component)
    );


    // Constructors

        //- Construct from components
        componentPhase
        (
            const dictionaryEntry& entry,
            const fvMesh& mesh,
            const alloyComponent& component
        );

        //- Helper for construction of componentPhase PtrList/PtrDictionary
        class iNew
        {
            const word& modelType_;
            const fvMesh& mesh_;
            const alloyComponent& component_;

        public:

            iNew
            (
                const word& modelType,
                const fvMesh& mesh,
                const alloyComponent& component
            )
            :
                modelType_(modelType),
                mesh_(mesh),
                component_(component)
            {}

            autoPtr<componentPhase> operator()(Istream& is) const
            {
                dictionaryEntry entry(dictionary::null, is);
                return autoPtr<componentPhase>
                (
                    componentPhase::New(entry, modelType_, mesh_, component_)
                );
            }
        };

        //- Return clone (need for PtrList)
        autoPtr<componentPhase> clone() const
        {
            notImplemented("componentPhase::clone() const");
            return autoPtr<componentPhase>();
        }


    //- Selector
    static autoPtr<componentPhase> New
    (
        const dictionaryEntry& entry,
        const word& modelType,
        const fvMesh& mesh,
        const alloyComponent& component
    );

    //- Destructor
    virtual ~componentPhase() = default;


    // Member Functions

        // Access

            //- Return name of the phase (need for PtrDictionary)
            const word& keyword() const
            {
                return name_;
            }

            //- Return the diffusion coefficients (solid and liquid) [m2/s]
            const dimensionedScalar& diffusion() const
            {
                return diffusion_;
            }

            //- Return the reference to the component
            const alloyComponent& component() const
            {
                return component_;
            }

        // Evaluation

            //- Return alloy component phase name
            word name() const;

            //- Return the equilibrium concentration field []
            virtual tmp<volScalarField> equilibrium(const volScalarField& T) const = 0;

            //- Return the equilibrium concentration []
            virtual dimensionedScalar equilibrium(const dimensionedScalar& T) const = 0;

            //- Return the slope along the solidification path field [K]
            virtual tmp<volScalarField> slope(const volScalarField& T) const = 0;

            //- Return the slope along the solidification path [K]
            virtual dimensionedScalar slope(const dimensionedScalar& T) const = 0;

            //- Return the derivative of the equilibrium concentration [1/K]
            template<class T1>
            auto equilibriumPrime(const T1& T) const -> decltype(T+T)
            {
                return 1/slope(T);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
