/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "thermoElasticPlastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermoElasticPlastic, 0);
    addToRunTimeSelectionTable(rheologyLaw, thermoElasticPlastic, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermoElasticPlastic::thermoElasticPlastic
(
    const word& name,
    const volSymmTensorField& sigma,
    const dictionary& dict
)
:
    elasticPlastic(name, sigma, dict),
    T_(sigma.mesh().lookupObject<volScalarField>(dict.lookup("temperatureField"))),
    tempVsSigmaY_(dict.subDict("sigmaYInterpolationTable"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::thermoElasticPlastic::sigmaY() const
{
    tmp<volScalarField> tresult(elasticPlastic::sigmaY());

    forAll(T_.internalField(), cellI) {
        tresult().internalField()[cellI] *= tempVsSigmaY_(T_.internalField()[cellI]);
    }

    tresult().correctBoundaryConditions();

    // use this line for debug: tresult().write();
    return tresult;
}

Foam::scalar
Foam::thermoElasticPlastic::sigmaY(const scalar epsilonPEq, const label cellID) const
{
    return elasticPlastic::sigmaY(epsilonPEq, cellID) *
        tempVsSigmaY_(T_.internalField()[cellID]);
}

// ************************************************************************* //
