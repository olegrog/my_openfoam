/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019-2020 Oleg Rogozin
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

#include "multicomponentAlloy.H"

#include "IOmanip.H"

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
    baseMolarMass_("baseMolarMass", dimMass/dimMoles, *this),
    latentHeat_("latentHeat", dimEnergy/dimMass, *this),
    interfaceEnergy_("interfaceEnergy", dimEnergy/dimArea, *this),
    liquidus_("liquidus", dimTemperature, *this),
    rhoLiquidus_("rhoLiquidus", dimDensity, *this),
    entropyChange_("entropyChange", latentHeat_*rhoLiquidus_/liquidus_),
    components_(lookup("components"), alloyComponent::iNew(mesh, liquidus_)),
    solidus_(liquidus_ - deltaTemp())
{
    Info<< "Alloy properties:" << endl
        << " -- liquidus (K) = " << liquidus_.value() << endl
        << " -- specified solidus (K) = " << get<dimensionedScalar>("solidus").value() << endl
        << " -- calculated solidus (K) = " << solidus_.value() << endl
        << " -- solidification interval (K) = " << (liquidus_ - solidus_).value() << endl
        << " -- molar mass (kg/mol) = " << molarMass().value() << endl
        << " -- entropy change (J/mol/K) = "
        << (entropyChange_*molarMass()/rhoLiquidus_).value() << endl
        << " -- entropy change (J/m^3/K) = " << entropyChange_.value() << endl
        << " -- Gibbs--Thomson coefficient (Km) = "
        << (interfaceEnergy_/entropyChange_).value() << endl
        << " -- liquid thermodynamic factor (J/m^3) = " << factor<0>().value() << endl
        << " -- solid thermodynamic factor (J/m^3) = " << factor<1>().value() << endl
        << " -- mean diffusion in liquid (m^2/s) = " << diffusionL().value() << endl
        << nl
        << "Component properties:" << endl;

    for (auto iter = components_.begin(); iter != components_.end(); ++iter)
    {
        Info<< " -- " << setw(2) << iter().name()
            << ": D_L = " << iter().diffusion<0>().value()
            << ", deltaA = " << (iter().deltaA()).value()
            << ", sqr(deltaA)/D_L = " << (sqr(iter().deltaA())/iter().diffusion<0>()).value()
            << endl;
    }

    Info<< nl;
}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::multicomponentAlloy::sumSqrA() const
{
    auto iter = components_.begin();

    dimensionedScalar result = sqr(iter().deltaA());

    for (++iter; iter != components_.end(); ++iter)
    {
        result += sqr(iter().deltaA());
    }

    return result;
}


Foam::dimensionedScalar Foam::multicomponentAlloy::sumD() const
{
    auto iter = components_.begin();

    dimensionedScalar result = iter().diffusion<0>();

    for (++iter; iter != components_.end(); ++iter)
    {
        result += iter().diffusion<0>();
    }

    return result;
}


Foam::dimensionedScalar Foam::multicomponentAlloy::sumSqrAperD() const
{
    auto iter = components_.begin();

    dimensionedScalar result = sqr(iter().deltaA())/iter().diffusion<0>();

    for (++iter; iter != components_.end(); ++iter)
    {
        result += sqr(iter().deltaA())/iter().diffusion<0>();
    }

    return result;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::multicomponentAlloy::molarMass() const
{
    auto iter = components_.begin();

    dimensionedScalar result = iter().equilibrium(0, liquidus_)/iter().molarMass();
    dimensionedScalar baseCeq = dimensionedScalar(1) - iter().equilibrium(0, liquidus_);

    for (++iter; iter != components_.end(); ++iter)
    {
        result += iter().equilibrium(0, liquidus_)/iter().molarMass();
        baseCeq -= iter().equilibrium(0, liquidus_);
    }

    return 1/(result + baseCeq/baseMolarMass_);
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
        result = result() + iter().deltaA()*(iter() - iter().equilibrium(phase, T));
    }

    return result()*factor<0>()/partition(phase);
}


// ************************************************************************* //
