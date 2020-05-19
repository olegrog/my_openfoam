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
    Function for updating GeometricField<Type> based on Type functor.

\*---------------------------------------------------------------------------*/

#ifndef updateGeometricField_H
#define updateGeometricField_H

namespace Foam
{

    //- Update internal and boundary fields
    //  Typical usage:
    //      updateGeometricField
    //      (
    //          f, [](scalar& f, scalar p1, scalar p2)
    //          {
    //              f = Op(p1, p2);
    //          },
    //          p1, p2
    //      );
    template
    <
        class Field,
        class Func,
        class... Args
    >
    void updateGeometricField
    (
        Field& f,
        Func calc,
        const Args&... args
    )
    {
        forAll(f, cellI)
        {
            calc(f[cellI], (args[cellI])...);
        }
        auto& fBf = f.boundaryFieldRef();
        forAll(fBf, patchi)
        {
            forAll(fBf[patchi], faceI)
            {
                calc(fBf[patchi][faceI], (args.boundaryField()[patchi][faceI])...);
            }
        }
    }

} // End namespace Foam

#endif
