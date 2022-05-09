/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2022 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of setZNDsolution.

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
    Foam::ZNDsystem

Description
    The ODE system for finding the ZND solution

SourceFiles
    ZNDsystem.C

\*---------------------------------------------------------------------------*/

#ifndef ZNDsystem_H
#define ZNDsystem_H

#include "ODESystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ZNDsystem Declaration
\*---------------------------------------------------------------------------*/

class ZNDsystem
:
    public ODESystem
{

private:

    // Private data

        //- Specific heat release [J/kg]
        const scalar Q_;

        //- Activation energy [J/kg]
        const scalar E_;

        //- Rate constant [1/s]
        const scalar k_;

        //- Specific heat ratio []
        const scalar g_;

        //- Chapman--Jouguet velocity [m/s]
        const scalar D_;

public:

    //- Constructor
    ZNDsystem(const dictionary& dict);

    //- Destructor
    virtual ~ZNDsystem() = default;


    // Member functions

        //- Return the velocity
        scalar U(scalar lambda) const;

        //- Return the pressure
        scalar p(scalar lambda) const;

        //- Return the density
        scalar rho(scalar lambda) const;

        //- Return the number of equations in the system
        virtual label nEqns() const;

        //- Calculate the derivatives in dydx
        virtual void derivatives
        (
            const scalar x,
            const scalarField& y,
            scalarField& dydx
        ) const final;

        //- Calculate the Jacobian of the system
        //  Need by the stiff-system solvers
        virtual void jacobian
        (
            const scalar x,
            const scalarField& y,
            scalarField& dfdx,
            scalarSquareMatrix& dfdy
        ) const final
        {
            NotImplemented;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
