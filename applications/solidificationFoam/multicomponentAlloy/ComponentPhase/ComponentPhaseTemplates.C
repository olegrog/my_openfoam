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

\*---------------------------------------------------------------------------*/

#include "ComponentPhase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class PhaseBoundary>
Foam::ComponentPhase<PhaseBoundary>::ComponentPhase
(
    const dictionaryEntry& entry,
    const fvMesh& mesh,
    const alloyComponent& component
)
:
    componentPhase(entry, mesh, component),
    phaseBoundary_(entry, *this)
{
    const label nPoints = 100;
    const scalar S = component.alloy().get<dimensionedScalar>("solidus").value();
    const scalar L = component.alloy().liquidus().value();

    // Scan the typical temperature range
    for (label i = -nPoints; i <= nPoints; i++)
    {
        const scalar T = S + (L - S)*i/nPoints;
        const scalar Ceq = phaseBoundary_.equilibrium(T);

        if (componentPhase::debug || debug)
        {
            Info<< " -- Ceq( " << T << " ) = " << Ceq << endl;
        }
    }
}


// ************************************************************************* //
