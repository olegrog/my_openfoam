/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
                isoAdvector | Copyright (C) 2016 DHI
              Modified work | Copyright (C) 2018 Johan Roenby
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
    Solver for thermo-fluid-dynamic model of SLM based on interIsoFoam.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

#include "../solidificationFoam/multicomponentAlloy/multicomponentAlloy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Func>
void calcGeomField(volScalarField& f, Func calc) {
    forAll(f, cellI) {
        calc(f.primitiveFieldRef(), cellI);
    }
    forAll(f.boundaryFieldRef(), patchi) {
        forAll(f.boundaryField()[patchi], faceI) {
            calc(f.boundaryFieldRef(false)[patchi], faceI);
        }
    }
}

void calcLiquidFraction(
    volScalarField& lf,
    const volScalarField& field,
    const dimensionedScalar& solidus,
    const dimensionedScalar& liquidus
)
{
    const scalar width = 0.5 * (liquidus - solidus).value();
    const scalar field1 = 0.5 * (liquidus + solidus).value();
    calcGeomField(lf, [&](scalarField& f, label i) {
        f[i] = 0.5 * Foam::tanh((field[i] - field1) / width) + 0.5;
    });
}

void calcLiquidFractionDer(
    volScalarField& dlf,
    const volScalarField& field,
    const dimensionedScalar& solidus,
    const dimensionedScalar& liquidus
)
{
    const scalar width = 0.5 * (liquidus - solidus).value();
    const scalar field1 = 0.5 * (liquidus + solidus).value();
    calcGeomField(dlf, [&](scalarField& f, label i) {
        f[i] = 0.5 / width / sqr(Foam::cosh((field[i] - field1) / width));
    });
}

volScalarField surfaceGaussian(
    const volVectorField& x,
    const dimensionedVector& x0,
    const dimensionedScalar& radius
)
{
    volVectorField r = x - x0;
    r.replace(2, 0);
    return 2 * exp(-2*magSqr((x - x0) / radius)) / constant::mathematical::pi / sqr(radius);
}

template<class T>
auto threePhaseParameter(const T& temp, const T& phi, const T& alpha,
    dimensionedScalar A_sol, dimensionedScalar A_liq,
    dimensionedScalar dA_sol, dimensionedScalar dA_liq,
    dimensionedScalar A_gas
) -> decltype(
    // can be removed in C++14
    ((A_sol + dA_sol*temp)*(1-phi) + (A_liq + dA_liq*temp)*phi)*(1-alpha) + A_gas*alpha
)
{
    return ((A_sol + dA_sol*temp)*(1-phi) + (A_liq + dA_liq*temp)*phi)*(1-alpha) + A_gas*alpha;
}

void calcEnthalpy(volScalarField& he, const volScalarField& temp, const volScalarField& phi,
    dimensionedScalar Cp_sol, dimensionedScalar Cp_liq,
    dimensionedScalar dCp_sol, dimensionedScalar dCp_liq,
    dimensionedScalar T_solidus, dimensionedScalar T_liquidus,
    dimensionedScalar enthalpyFusion
)
{
    const scalar T1 = (T_solidus + T_liquidus).value() / 2;
    auto he_S = [&](scalar temp) {
        return Cp_sol.value()*temp + dCp_sol.value()*sqr(temp)/2;
    };
    auto he_L = [&](scalar temp) {
        return Cp_liq.value()*(temp-T1) + dCp_liq.value()*sqr(temp-T1)/2;
    };
    calcGeomField(he, [&](scalarField& f, label i) {
        if (temp[i] <= T1) {
            f[i] = he_S(temp[i]);
        } else {
            f[i] = he_S(T1) + he_L(temp[i]);
        }
        f[i] += phi[i] * enthalpyFusion.value();
    });
    he.correctBoundaryConditions();
}

void calcTemperature(volScalarField& temp, const volScalarField& he, const volScalarField& phi,
    dimensionedScalar Cp_sol, dimensionedScalar Cp_liq,
    dimensionedScalar dCp_sol, dimensionedScalar dCp_liq,
    dimensionedScalar T_solidus, dimensionedScalar T_liquidus,
    dimensionedScalar enthalpyFusion
)
{
    const dimensionedScalar T1 = (T_solidus + T_liquidus) / 2;
    auto he_S = [=](dimensionedScalar temp) {
        return Cp_sol*temp + dCp_sol*sqr(temp)/2;
    };
    const dimensionedScalar he1 = he_S(T1) + enthalpyFusion/2;
    forAll(temp, cellI) {
        scalar A, B, C = he[cellI] - phi[cellI]*enthalpyFusion.value();
        if (he[cellI] <= he1.value()) {
            A = dCp_sol.value();
            B = Cp_sol.value();
        } else {
            A = dCp_liq.value();
            B = (Cp_liq - dCp_liq*T1).value();
            C += ((Cp_liq - Cp_sol)*T1 - (dCp_sol + dCp_liq)*sqr(T1)/2).value();
        }
        if (A > SMALL) {
            temp[cellI] = (Foam::sqrt(sqr(B) + 2*A*C) - B) / A;
        } else {
            temp[cellI] = C / B;
        }
    }
    temp.correctBoundaryConditions();
}

tmp<volVectorField> calcNormal(const volVectorField& vField) {
    dimensionedVector smallVector("small", vField.dimensions(), vector(0, 0, SMALL));
    return (vField + smallVector) / mag(vField + smallVector);
}

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    calcTemperature(T, he, liquidFraction, Cp_sol, Cp_liq, dCp_sol, dCp_liq, T_solidus, T_liquidus, enthalpyFusion);

    // -- Initial conditions for concentrations
    if (segregation) {
        forAllIter(PtrDictionary<alloyComponent>, pAlloy->components(), iter) {
            alloyComponent& C = iter();
            C == C.equilibrium(0, pAlloy->liquidus());
        }
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }


                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            vaporPressure = ambientPressure * exp(molarMass * enthalpyBoiling
                / constant::physicoChemical::R * (1/T_boiling - 1/T));

            #include "heEqn.H"

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        if (segregation) {
            #include "CEqn.H"
        }

        kinematicViscosity = mixture.nu();
        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
