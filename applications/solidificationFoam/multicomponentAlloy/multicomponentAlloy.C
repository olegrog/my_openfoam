/*---------------------------------------------------------------------------*\
  =========                 |
  \\     / F ield         | foam-extend: Open Source CFD
   \\   /  O peration     | Version:     4.0
    \\ /   A nd           | Web:         http://www.foam-extend.org
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

#include "multicomponentAlloy.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multicomponentAlloy::multicomponentAlloy(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "alloyProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    entropyChange_(get<dimensionedScalar>("entropyChange")),
    interfaceEnergy_(get<dimensionedScalar>("interfaceEnergy")),
    solidus_(get<dimensionedScalar>("solidus")),
    liquidus_(get<dimensionedScalar>("liquidus")),
    components_(lookup("components"), alloyComponent::iNew(mesh, liquidus_)),
    factorS_(calcFactor<0>()),
    factorL_(calcFactor<1>())
{
    for (auto iter = components_.begin(); iter != components_.end(); ++iter)
    {
        Info<< "Component " << iter().name() << ": D_L = " << iter().diffusion<1>().value()
            << ", deltaA = " << (iter().deltaA()).value()
            << ", sqr(deltaA)/D_L = " << (sqr(iter().deltaA())/iter().diffusion<1>()).value()
            << endl;
    }

    Info << nl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::multicomponentAlloy::relaxationTime() const
{
    auto iter = components_.begin();

    dimensionedScalar result = sqr(iter().deltaA())/iter().diffusion<1>();

    for (++iter; iter != components_.end(); ++iter)
    {
        result += sqr(iter().deltaA())/iter().diffusion<1>();
    }

    return result*factorL_/interfaceEnergy_;
}


Foam::dimensionedScalar Foam::multicomponentAlloy::capillaryLength() const
{
    auto iter = components_.begin();

    dimensionedScalar result = iter().diffusion<1>();

    for (++iter; iter != components_.end(); ++iter)
    {
        result += iter().diffusion<1>();
    }

    return 1./result/relaxationTime();
}


Foam::dimensionedScalar Foam::multicomponentAlloy::diffusionL() const
{
    auto iter = components_.begin();

    dimensionedScalar result = sqr(iter().deltaA());

    for (++iter; iter != components_.end(); ++iter)
    {
        result += sqr(iter().deltaA());
    }

    return result*factorL_/relaxationTime()/interfaceEnergy_;
}


Foam::tmp<Foam::volScalarField> Foam::multicomponentAlloy::chemicalDrivingForce
(
    const volScalarField& phase,
    const volScalarField& T
) const
{
    auto iter = components_.begin();

    tmp<volScalarField> result = iter().deltaA()*(iter() - iter().equilibrium(phase, T));

    for (++iter; iter != components_.end(); ++iter)
    {
        result = result() + iter().deltaA() * (iter() - iter().equilibrium(phase, T));
    }

    return result()*factorL_/partition(phase);
}


// ************************************************************************* //
