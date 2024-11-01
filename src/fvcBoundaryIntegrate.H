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

InNamespace
    Foam::fvc

Description
    Integrate volField over the whole boundary creating a dimensioned\<Type\>

\*---------------------------------------------------------------------------*/

#ifndef fvcBoundaryIntegrate_H
#define fvcBoundaryIntegrate_H

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    template<class Type>
    dimensioned<Type> boundaryIntegrate
    (
        const GeometricField<Type, fvPatchField, volMesh>& vf
    )
    {
        Type integral = Zero;
        const auto& vfBf = vf.boundaryField();
        const auto& meshMagSfBf = vf.mesh().magSf().boundaryField();
        forAll(vfBf, patchi)
        {
            integral += gSum(vfBf[patchi]*meshMagSfBf[patchi]);
        }
        return dimensioned<Type>
        (
            "boundaryIntegrate(" + vf.name() + ')',
            dimArea*vf.dimensions(),
            integral
        );
    }

    template<class Type>
    dimensioned<Type> boundaryIntegrate
    (
        const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
    )
    {
        dimensioned<Type> integral = boundaryIntegrate(tvf());
        tvf.clear();
        return integral;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
