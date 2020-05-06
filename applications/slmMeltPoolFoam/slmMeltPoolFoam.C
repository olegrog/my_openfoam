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
            slmMeltPoolFoam | Copyright (C) 2019-2020 Oleg Rogozin
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

#include "LiquidFraction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// We use this function instead of GeometricField::operator= to
//  1) change boundary values along with internal ones,
//  2) optimize piecewise functions
template<class Func, class... Args>
void calcGeomField(volScalarField& f, Func calc, const Args&... args)
{
    forAll(f, cellI)
    {
        calc(f.primitiveFieldRef(), cellI, (args.primitiveField())...);
    }
    forAll(f.boundaryField(), patchi)
    {
        forAll(f.boundaryField()[patchi], faceI)
        {
            calc(f.boundaryFieldRef(false)[patchi], faceI,
                (args.boundaryField()[patchi])...);
        }
    }
}

volScalarField surfaceGaussian
(
    const volVectorField& x,
    dimensionedVector x0,
    dimensionedScalar radius
)
{
    volVectorField r = x - x0;
    r.replace(2, 0);
    return 2*exp(-2*magSqr((x - x0)/radius))/constant::mathematical::pi/sqr(radius);
}

template<class T1, class T2, class T3>
auto threePhaseParameter
(
    const T1& temp, const T2& phi, const T3& alpha,
    dimensionedScalar A_sol, dimensionedScalar A_liq,
    dimensionedScalar dA_sol, dimensionedScalar dA_liq,
    dimensionedScalar A_gas
) -> decltype(temp + phi + alpha) // can be removed in C++14
{
    return (A_sol + dA_sol*temp)*(1-phi) + (A_liq + dA_liq*temp)*phi + A_gas*alpha;
}

void calcEnthalpy
(
    volScalarField& he, const volScalarField& temp,
    const LiquidFraction& phi, const volScalarField& alpha,
    dimensionedScalar Cp_sol, dimensionedScalar dCp_sol,
    dimensionedScalar Cp_liq, dimensionedScalar dCp_liq,
    dimensionedScalar Cp_gas, dimensionedScalar T_melting,
    dimensionedScalar enthalpyFusion, bool dAlpha = false
)
{
    // To check dimensions
    (Cp_gas + Cp_sol + Cp_liq + temp*(dCp_sol + dCp_liq))*T_melting/(he+enthalpyFusion) + phi + alpha;
    const scalar T_M = T_melting.value();
    const scalar he_fus = enthalpyFusion.value();
    auto he_S = [=](scalar temp) {
        return Cp_sol.value()*temp + dCp_sol.value()*sqr(temp)/2;
    };
    auto he_L = [=](scalar temp) {
        return Cp_liq.value()*temp + dCp_liq.value()*sqr(temp)/2;
    };
    const volScalarField he_G = Cp_gas*temp;

    calcGeomField(he, [=](scalarField& f, label i, const scalarField& temp,
        const scalarField& phi, const scalarField& alpha, const scalarField& he_G)
    {
        scalar piecewise;
        if (temp[i] <= T_M) {
            piecewise = he_S(temp[i]);
        } else {
            piecewise = he_S(T_M) + he_L(temp[i]) - he_L(T_M);
        }
        if (!dAlpha) {
            f[i] = alpha[i]*he_G[i] + he_fus*phi[i] + (1-alpha[i])*piecewise;
        } else {
            f[i] = he_G[i] - piecewise;
        }
    }, temp, phi, alpha, he_G);
    he.correctBoundaryConditions();
}

void calcTemperature
(
    volScalarField& temp, const volScalarField& he,
    const LiquidFraction& phi, const volScalarField& alpha,
    dimensionedScalar Cp_sol, dimensionedScalar dCp_sol,
    dimensionedScalar Cp_liq, dimensionedScalar dCp_liq,
    dimensionedScalar Cp_gas, dimensionedScalar T_melting,
    dimensionedScalar enthalpyFusion
)
{
    // to check dimensions
    (Cp_gas + Cp_sol + Cp_liq + temp*(dCp_sol + dCp_liq))*T_melting/(he+enthalpyFusion) + phi + alpha;
    const scalar T_M = T_melting.value();
    const scalar he_fus = enthalpyFusion.value();
    auto he_S = [=](scalar temp) {
        return Cp_sol.value()*temp + dCp_sol.value()*sqr(temp)/2;
    };
    auto he_L = [=](scalar temp) {
        return Cp_liq.value()*temp + dCp_liq.value()*sqr(temp)/2;
    };
    auto he_G = [=](scalar temp) {
        return Cp_gas.value()*temp;
    };

    calcGeomField(temp, [=](scalarField& f, label i, const scalarField& he,
        const scalarField& phi, const scalarField& alpha)
    {
        scalar A, B = alpha[i]*Cp_gas.value();
        scalar C = he[i] - he_fus*phi[i] - (1-alpha[i])*he_S(T_M);
        scalar he_M = alpha[i]*he_G(T_M) + (1-alpha[i])*(he_S(T_M) + he_fus/2);
        if (he[i] < he_M) {
            A = (1-alpha[i])*dCp_sol.value();
            B += (1-alpha[i])*Cp_sol.value();
            C += (1-alpha[i])*he_S(T_M);
        } else {
            A = (1-alpha[i])*dCp_liq.value();
            B += (1-alpha[i])*Cp_liq.value();
            C += (1-alpha[i])*he_L(T_M);
        }
        if (mag(A) > SMALL) {
            f[i] = (Foam::sqrt(sqr(B) + 2*A*C) - B)/A;
        } else {
            f[i] = C/B;
        }
    }, he, phi, alpha);
    temp.correctBoundaryConditions();
}

tmp<volVectorField> calcNormal(const volVectorField& vField)
{
    const dimensionedVector smallVector("small", vField.dimensions(), vector(0, 0, SMALL));
    return (vField + smallVector) / mag(vField + smallVector);
}

tmp<volScalarField> magModifiedGradAlpha
(
    const volScalarField& temp, const LiquidFraction& phi, const volScalarField& alpha,
    dimensionedScalar A_sol, dimensionedScalar A_liq,
    dimensionedScalar dA_sol, dimensionedScalar dA_liq,
    dimensionedScalar A_gas
)
{
    const dimensionedScalar zero(0), one(1);
    return 2*mag(fvc::grad(alpha)) *
        threePhaseParameter(temp, phi, alpha, A_sol, A_liq, dA_sol, dA_liq, A_gas)
    / (
        threePhaseParameter(temp, phi, zero, A_sol, A_liq, dA_sol, dA_liq, A_gas)
        +
        threePhaseParameter(temp, phi, one, A_sol, A_liq, dA_sol, dA_liq, A_gas)
    );
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

    // For debug
    auto pp = [](const volScalarField& f)
    {
        Info<< f.name() << ": " << min(f).value() << " " << max(f).value() << endl;
    };

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

        const dimensionedVector laserCoordinate
        (
            "laserCoordinate", coordStart + laserVelocity * mesh.time()
        );
        const dimensionedScalar totalEnthalpy = fvc::domainIntegrate(rho*he);
        // vector beamDirection(0, 0, -1);
        // mag(fvc::grad(alpha2) & beamDirection) can be used instead of mag(fvc::grad(alpha2))

        // --- Enthalpy--pressure--velocity PIMPLE corrector loop
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

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            // Enthalpy-independent quantities
            heAtFusion = he_melting(alpha2);
            laserHeatSource = (runTime < timeStop) * absorptivity * laserPower
                * surfaceGaussian(mesh.C(), laserCoordinate, laserRadius);

            while (pimple.correct())
            {
                #include "heEqn.H"
            }

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

        // Update passive quantities
        kinematicViscosity = mixture.nu();
        liquidFraction.finalUpdate();

        Info<< "Real energy input = " << (fvc::domainIntegrate(rho*he) - totalEnthalpy).value()
            << ", laser input = "
            << fvc::domainIntegrate(laserHeatSource * mag(fvc::grad(alpha2)) * runTime.deltaT()).value()
            << ", theoretical value = " << (absorptivity * laserPower * runTime.deltaT()).value()
            << endl;

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
