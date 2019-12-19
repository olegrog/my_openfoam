/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
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
    Solver for macroscopic model of SLM based on interFoam.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
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

template<class T>
T tanhSmooth(const T& x, const dimensionedScalar& x0, const dimensionedScalar& width)
{
    return 0.5 * tanh((x - x0) / width) + 0.5;
}

volScalarField gaussian(const volVectorField& x, const dimensionedVector& x0, const dimensionedScalar& radius)
{
    volVectorField r = x - x0;
    r.replace(2, 0);
    return 2 * exp(-2*magSqr((x - x0) / radius)) / constant::mathematical::pi / sqr(radius);
}

volScalarField fourParameterModel (const volScalarField& phi, const volScalarField& T,
    dimensionedScalar A_sol, dimensionedScalar A_liq,
    dimensionedScalar dA_sol, dimensionedScalar dA_liq, dimensionedScalar T0)
{
    return (A_sol+dA_sol*(T-T0))+phi*((A_liq-A_sol)+(dA_liq-dA_sol)*(T-T0));
}

template<class T>
T enthalpyCalc(const T& temp, const T& phi,
    dimensionedScalar Cp_sol, dimensionedScalar Cp_liq,
    dimensionedScalar dCp_sol, dimensionedScalar dCp_liq,
    dimensionedScalar T_solidus, dimensionedScalar T_liquidus,
    dimensionedScalar enthalpyFusion)
{
    return ((0.5*temp*(2*Cp_sol+dCp_sol*temp))*(scalar(1)-phi)
    + phi*(0.5*(temp-T_liquidus)*(2*Cp_liq+dCp_liq*(temp)) + enthalpyFusion + 0.5*T_solidus*(2*Cp_sol+dCp_sol*T_solidus)));
}

template<class T>
T temperatureCalc (const T& he, const T& phi,
    dimensionedScalar Cp_sol, dimensionedScalar Cp_liq,
    dimensionedScalar dCp_sol, dimensionedScalar dCp_liq,
    dimensionedScalar T_solidus, dimensionedScalar T_liquidus,
    dimensionedScalar enthalpyFusion)
{
    volScalarField c (phi*(enthalpyFusion + 0.5*T_solidus*(2*Cp_sol + dCp_sol*T_solidus) - T_liquidus*Cp_liq) - he);
    volScalarField b (phi*(Cp_liq - (T_liquidus*dCp_liq)/2) - Cp_sol*(phi - 1));
    volScalarField a (((dCp_liq*phi)/2 - (dCp_sol*(phi - 1))/2));
    volScalarField D (b*b - 4*a*c);
    volScalarField ans ((-b+Foam::sqrt(D))/(2*a));
    return ans;
}

tmp<volVectorField> calcNormal(const volVectorField& vField) {
    dimensionedVector smallVector("small", vField.dimensions(), vector(0, 0, SMALL));
    return (vField + smallVector) / mag(vField + smallVector);
}

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    dimensionedScalar phi0 = tanhSmooth(T0, (T_solidus + T_liquidus)/2, (T_liquidus - T_solidus)/2);
    Info<< "Initial enthalpy: "
        << enthalpyCalc(T0, phi0, Cp_sol, Cp_liq, dCp_sol, dCp_liq, T_solidus, T_liquidus, enthalpyFusion).value()
        << endl;

    T = temperatureCalc(he, liquidFraction, Cp_sol, Cp_liq, dCp_sol, dCp_liq, T_solidus, T_liquidus, enthalpyFusion);

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

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Temperature equation
        volScalarField Cp = fourParameterModel(liquidFraction, T, Cp_sol, Cp_liq, dCp_sol, dCp_liq, T0);
        volScalarField k = fourParameterModel(liquidFraction, T, k_sol, k_liq, dk_sol, dk_liq, T0) * alpha1;
        diffusivity = k / Cp;

        dimensionedVector laserCoordinate("laserCoordinate", coordStart + laserVelocity * mesh.time());
        volScalarField laserHeatSource = mag(fvc::grad(alpha1))
            * absorptivity * laserPower * gaussian(mesh.C(), laserCoordinate, laserRadius);
        volScalarField fusionTerm = fvc::laplacian(diffusivity * enthalpyFusion, liquidFraction);

        fvScalarMatrix heEqn
        (
            fvm::ddt(rho, he)
            + fvm::div(rhoPhi, he)
            - fvm::laplacian(diffusivity, he)
            - fvc::div(rho * U & turbulence->devReff())
            - (runTime < timeStop) * laserHeatSource
            + fusionTerm
        );
        heEqn.solve();

        liquidFraction = tanhSmooth(he, (he_solidus + he_liquidus)/2, (he_liquidus - he_solidus)/2);
        wasMelted = max(wasMelted, liquidFraction);
        T = temperatureCalc(he, liquidFraction, Cp_sol, Cp_liq, dCp_sol, dCp_liq, T_solidus, T_liquidus, enthalpyFusion);

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

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
                    }

                    mixture.correct();

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            // -- Calculate concentrations
            if (segregation) {
                volScalarField phase = 1 - liquidFraction;
                volScalarField h = pAlloy->partition(phase);
                forAllIter(PtrDictionary<alloyComponent>, pAlloy->components(), iter) {
                    alloyComponent& C = iter();
                    surfaceScalarField phi = fvc::snGrad(phase) * mesh.magSf()
                        * fvc::interpolate(C.diffusion(phase) / sqr(h)) * pAlloy->partitionPrime();
                    fvScalarMatrix CEqn(
                        fvm::ddt(C) == fvm::laplacian(C.diffusion(phase) / h, C) - fvm::div(phi, C)
                        - fvc::laplacian(C.diffusion(phase), C.equilibrium(phase, T) / h)
                    );
                    CEqn.solve();
                }
            }
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

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

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
