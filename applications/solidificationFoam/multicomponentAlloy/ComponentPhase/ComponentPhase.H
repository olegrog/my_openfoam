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
    Foam::ComponentPhase

Description
    Class template for an arbitrary phase boundary of the alloy component.

\*---------------------------------------------------------------------------*/

#ifndef ComponentPhase_H
#define ComponentPhase_H

#include "componentPhase.H"
#include "generateGeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ComponentPhase Declaration
\*---------------------------------------------------------------------------*/

template<class PhaseBoundary>
class ComponentPhase
:
    public componentPhase
{
    // Private Data

        //- Phase boundary
        const PhaseBoundary phaseBoundary_;

public:

    //- Runtime type information
    TypeName("ComponentPhase");

    //- Constructor
    ComponentPhase
    (
        const dictionaryEntry& entry,
        const fvMesh& mesh,
        const alloyComponent& component
    )
    :
        componentPhase(entry, mesh, component),
        phaseBoundary_(entry, component)
    {}

    //- Destructor
    ~ComponentPhase() final = default;


    // Member Functions

        //- Return the equilibrium concentration field []
        tmp<volScalarField> equilibrium(const volScalarField& T) const final
        {
            return generateGeometricField<volScalarField>
            (
                "equilibrium_" + name(),
                mesh_,
                dimless,
                [this](scalar T)
                {
                    return phaseBoundary_.equilibrium(T);
                },
                T
            );
        }

        //- Return the equilibrium concentration []
        dimensionedScalar equilibrium(const dimensionedScalar& T) const final
        {
            return dimensionedScalar
            (
                "equilibrium_" + name(),
                dimless,
                phaseBoundary_.equilibrium(T.value())
            );
        }

        //- Return the slope along the solidification path field [K]
        tmp<volScalarField> slope(const volScalarField& T) const final
        {
            return generateGeometricField<volScalarField>
            (
                "slope_" + name(),
                mesh_,
                dimTemperature,
                [this](scalar T)
                {
                    return phaseBoundary_.slope(T);
                },
                T
            );
        }

        //- Return the slope along the solidification path [K]
        dimensionedScalar slope(const dimensionedScalar& T) const final
        {
            return dimensionedScalar
            (
                "slope_" + name(),
                dimTemperature,
                phaseBoundary_.slope(T.value())
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
