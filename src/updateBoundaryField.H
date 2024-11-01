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
    This file is a missed part of OpenFOAM.

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

Description
    Function for updating GeometricField<Type>.boundaryField() based on Type
    functor with normal .

\*---------------------------------------------------------------------------*/

#ifndef updateBoundaryField_H
#define updateBoundaryField_H

namespace Foam
{

    //- Update the boundary field only.
    //  Typical usage:
    //      updateBoundaryField
    //      (
    //          f, [](scalar& f, const vector& n, scalar p1, scalar p2)
    //          {
    //              f = Op(n, p1, p2);
    //          },
    //          p1, p2
    //      );
    //
    //  \note Dimensions of fields are not checked.
    template
    <
        class Field,
        class Func,
        class... Args
    >
    void updateBoundaryField
    (
        Field& f,
        Func calc,
        const Args&... args
    )
    {
        auto& fBf = f.boundaryFieldRef();
        forAll(fBf, patchi)
        {
            const vectorField n(fBf[patchi].patch().nf());
            forAll(fBf[patchi], faceI)
            {
                calc(fBf[patchi][faceI], n[faceI], (args.boundaryField()[patchi][faceI])...);
            }
        }
    }

} // End namespace Foam

#endif
