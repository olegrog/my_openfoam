/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019-2021 Oleg Rogozin
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
    rhoSolid_("rhoSolid", dimDensity, *this),
    entropyChange_("entropyChange", latentHeat_*rhoSolid_/liquidus_),
    phaseNames_(get<wordList>("phases")),
    components_(lookup("components"), alloyComponent::iNew(mesh, *this)),
    solidus_(liquidus_ - deltaTemp())
{
    Info<< "\nAlloy properties:" << endl
        << " -- liquidus (K) = " << liquidus_.value() << endl
        << " -- specified solidus (K) = " << get<dimensionedScalar>("solidus").value() << endl
        << " -- calculated solidus (K) = " << solidus_.value() << endl
        << " -- solidification interval (K) = " << (liquidus_ - solidus_).value() << endl
        << " -- molar mass (kg/mol) = " << molarMass().value() << endl
        << " -- entropy change (J/mol/K) = "
        << (entropyChange_*molarMass()/rhoSolid_).value() << endl
        << " -- entropy change (J/m^3/K) = " << entropyChange_.value() << endl
        << " -- Gibbs--Thomson coefficient (Km) = "
        << (interfaceEnergy_/entropyChange_).value() << endl;

    for (const word& phaseName : phaseNames_)
    {
        Info<< " -- thermodynamic factor in " << phaseName << " at liquidus (J/m^3) = "
            << factor(phaseName, liquidus_).value() << endl
            << " -- thermodynamic factor in " << phaseName << " at solidus (J/m^3) = "
            << factor(phaseName, solidus_).value() << endl;
    }

    PtrList<Tuple2<word, const dimensionedScalar&>> pairs;
    pairs.emplace_append("solidus", solidus_);
    pairs.emplace_append("liquidus", liquidus_);

    for (auto pair : pairs)
    {
        Info<< " -- mean diffusion in liquid at " << pair.first() << " (m^2/s) = "
            << diffusionL(pair.second()).value() << endl;
    }

    for (auto pair : pairs)
    {
        Info<< "\nComponent properties at " << pair.first() << ":" << endl;
        for (const alloyComponent& component : components_)
        {
            const dimensionedScalar& T = pair.second();
            const scalar delta = component.delta(T).value();
            const scalar D_L = component.phase("liquid").diffusion().value();
            const scalar slopeS = component.phase("solid").slope(T).value();
            const scalar slopeL = component.phase("liquid").slope(T).value();

            Info<< " -- " << setw(2) << component.name()
                << ": D_L = " << D_L
                << ", delta = " << delta
                << ", sqr(delta)/D_L = " << sqr(delta)/D_L
                << ", delta/mS = " << delta/slopeS
                << ", delta/mL = " << delta/slopeL
                << endl;
        }
    }

    Info<< nl;
}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::multicomponentAlloy::sumDL() const
{
    auto iter = components_.begin();

    dimensionedScalar result = iter().phase("liquid").diffusion();

    for (++iter; iter != components_.end(); ++iter)
    {
        result += iter().phase("liquid").diffusion();
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

    auto result = iter().delta(T)*(iter() - iter().equilibrium(phase, T));

    for (++iter; iter != components_.end(); ++iter)
    {
        result = result() + iter().delta(T)*(iter() - iter().equilibrium(phase, T));
    }

    return result()*factor("liquid", T)/partition(phase, T);
}


// ************************************************************************* //
