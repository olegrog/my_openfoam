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
    Function for generating GeometricField<Type> based on Type functor.

\*---------------------------------------------------------------------------*/

#ifndef generateGeometricField_H
#define generateGeometricField_H

#include "tmp.H"
#include "dimensionSet.H"

namespace Foam
{

    //- Generate internal and boundary fields.
    //  Typical usage:
    //      generateGeometricField
    //      (
    //          "field",
    //          mesh,
    //          dimless,
    //          [](scalar p1, scalar p2)
    //          {
    //              return Op(p1, p2);
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
    tmp<Field> generateGeometricField
    (
        const word& name,
        const typename Field::Mesh& mesh,
        const dimensionSet& ds,
        Func calc,
        const Args&... args
    )
    {
        auto tres = Field::New(name, mesh, ds);

        auto& f = tres.ref();
        forAll(f, cellI)
        {
            f[cellI] = calc((args[cellI])...);
        }

        auto& fBf = f.boundaryFieldRef();
        forAll(fBf, patchi)
        {
            forAll(fBf[patchi], faceI)
            {
                fBf[patchi][faceI] = calc((args.boundaryField()[patchi][faceI])...);
            }
        }

        return tres;
    }

} // End namespace Foam

#endif
