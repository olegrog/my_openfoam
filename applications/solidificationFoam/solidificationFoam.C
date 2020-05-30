/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.2.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
-------------------------------------------------------------------------------
Application
    solidificationFoam

Description
    Solver for the phase-field equations for the solidification of binary alloys.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "zeroGradientFvPatchField.H"

#include "multicomponentAlloy/multicomponentAlloy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using constant::mathematical::pi;

tmp<volScalarField> fPrime(const volScalarField& phase)
{
    return 2*phase*(1 - phase)*(1 - 2*phase);
}

tmp<volScalarField> fPrime2(const volScalarField& phase)
{
    return 3*sqr(1 - 2*phase) - 1;
}

tmp<volScalarField> gPrime(const volScalarField& phase)
{
    return 30*sqr(phase)*sqr(1 - phase);
}

tmp<volScalarField> gPrime2(const volScalarField& phase)
{
    return 30*fPrime(phase);
}

tmp<volScalarField> generateSeed
(
    const volVectorField& coord,
    const dimensionedVector& center,
    const dimensionedScalar& radius
)
{
    return min
    (
        2*Foam::exp(-magSqr((coord - center)/radius)),
        scalar(1)
    );
}

void addGrain(volVectorField& grain, const volScalarField& phase, label nGrain, label nGrains)
{
    forAll(grain, cellI)
    {
        scalar argument = 2*pi*sign(phase[cellI])*nGrain/nGrains;
        scalar magnitude = fabs(phase[cellI]);
        grain[cellI].x() += magnitude*Foam::cos(argument);
        grain[cellI].y() += magnitude*Foam::sin(argument);
    }
}

void calcNGrain(volScalarField& nGrain, const volVectorField& grain, label nGrains)
{
    forAll(grain, cellI)
    {
        nGrain[cellI] = Foam::atan2(grain[cellI].y(), grain[cellI].x())/2/pi*nGrains;
    }
}

tmp<volVectorField> calcNormal(const volVectorField& vField)
{
    dimensionedVector smallVector("small", vField.dimensions(), vector(0, SMALL, 0));
    return (vField + smallVector)/mag(vField + smallVector);
}

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // --- Additional variables

    volScalarField theta = 0*phase;
    volScalarField theta0 = 0*phase;
    const tensor rot(0, -1, 0, 1, 0, 0, 0, 0, 1);

    // Derived quantities
    const dimensionedScalar tau = a1*a2*pow3(interfaceWidth)*alloy.relaxationTime();
    dimensionedScalar tipVelocity = coolingRate/tempGradient;
    const label nGrains = crystallographicAngles.size();

    // Coordinates-related constants and variables
    const volVectorField coord = mesh.C();
    const boundBox& bounds = mesh.bounds();
    const dimensionedScalar xmax("xmax", dimLength, bounds.max().x());
    const dimensionedScalar xmin("xmin", dimLength, bounds.min().x());
    const dimensionedScalar ymax("ymax", dimLength, bounds.max().y());
    const dimensionedScalar ymin("ymin", dimLength, bounds.min().y());
    const dimensionedScalar height("height", ymax - ymin);
    const dimensionedScalar width("width", xmax - xmin);
    const dimensionedVector center("center", dimLength, (bounds.max() + bounds.min())/2);
    const dimensionedScalar frontPosition = ymin + frontPositionRel*height;
    const dimensionedScalar initialWidth = interfaceWidth/interfaceNarrowing;
    dimensionedScalar tipPosition = frontPosition;
    dimensionedScalar tipPositionPrev = tipPosition;
    label tipCell(0), tipCellPrev(0);
    dimensionedScalar tipTimeInterval("tipTimeInterval", dimTime, 0);

    // --- Initial conditions

    // phase + grain
    phase = Foam::atan(pow3((frontPosition - coord.component(vector::Y))/initialWidth))/pi + .5;
    const dimensionedScalar radius = width/nSeeds/seedNarrowing;
    for (int i = 0; i < nSeeds; i++)
    {
        dimensionedVector position = center;
        position.replace(vector::X, xmin + (i+.5)/nSeeds*width);
        position.replace(vector::Y, frontPosition);
        phase = max(phase, generateSeed(coord, position, radius));
    }
    phase.clip(0, 1);

    const volScalarField grainNumber
    (
        sign((coord.component(vector::X) - center.component(vector::X))/width)
    );
    addGrain(grain, phase*grainNumber, 1, nGrains);

    if (addRandomSeeds) {
        dimensionedVector position = center;
        position.replace(vector::Y, ymin/3 + 2*ymax/3);
        volScalarField seed = generateSeed(coord, position, radius/2);
        phase += seed;
        addGrain(grain, seed, 0, nGrains);
    }

    // number of grain
    calcNGrain(nGrain, grain, nGrains);

    // temperature
    T = alloy.liquidus() - undercooling
        + tempGradient*(coord.component(vector::Y) - ymin/3 - 2*ymax/3);

    // concentrations
    forAllIter(PtrDictionary<alloyComponent>, alloy.components(), iter)
    {
        alloyComponent& C = iter();
        const volScalarField temp = alloy.solidus()*phase + alloy.liquidus()*(1 - phase);
        C == C.equilibrium(phase, temp);
    }

    // --- Print reference parameters

    const scalar minMeshStep = 1/gMax(mesh.surfaceInterpolation::deltaCoeffs());
    const scalar sigmaW =
    (
        alloy.diffusionL()*alloy.capillaryLength()/tipVelocity/sqr(interfaceWidth)
    ).value();
    const scalar minSigmaW(runTime.controlDict().get<scalar>("minSigmaW"));
    const dimensionedScalar Tmin("Tmin", dimTemperature, gMin(T));
    const dimensionedScalar Tmax("Tmax", dimTemperature, gMax(T));
    const dimensionedScalar PDAS = (xmax - xmin)/nSeeds;

    Info<< "Dimensionless parameters:" << endl
        << " -- minimal undercooling = "
        << alloy.undercooling(Tmin).value() << endl
        << " -- maximum undercooling = "
        << alloy.undercooling(Tmax).value() << endl
        << " -- interface width / capillary length = "
        << (interfaceWidth/alloy.capillaryLength()).value() << endl
        << " -- mesh step / interface width = "
        << (minMeshStep/interfaceWidth).value() << endl
        << " -- theoretical tip velocity = "
        << (tipVelocity*alloy.capillaryLength()/alloy.diffusionL()).value() << endl
        << " -- interface stability parameter = "
        << sigmaW << endl
        << " -- Reynolds number = "
        << max(mag(U)/laminarTransport.nu()*PDAS).value() << endl
        << nl;

    Info<< "Dimensioned parameters:" << endl
        << " -- generated PDAS = "
        << PDAS.value() << endl
        << " -- theoretical tip velocity (m/s) = "
        << tipVelocity.value() << endl
        << " -- characteristic velocity (m/s) = "
        << (alloy.diffusionL()/alloy.capillaryLength()).value() << endl
        << " -- capillary length (m) = "
        << alloy.capillaryLength().value() << endl
        << " -- relaxation time (s) = "
        << tau.value() << endl;

    if (sigmaW < minSigmaW)
    {
        FatalError
            << "Interface stability parameter is too low:" << endl
            << sigmaW << " < " << minSigmaW
            << abort(FatalError);
    }

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "phaseCourantNo.H"
        #include "setDeltaT.H"

        scalar minDeltaT =
            runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", 0);
        if (runTime.deltaTValue() < minDeltaT)
        {
            Info<< "Time step becomes too small!" << endl;
            runTime.writeAndEnd();
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "tipVelocity.H"

        // --- Calculate temperature

        T = alloy.liquidus() - undercooling
            + tempGradient*(coord.component(vector::Y) - ymin/3 - 2*ymax/3)
            - coolingRate*runTime;

        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    if (correctPhi)
                    {
                        #include "correctPhi.H"
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            label nCorrPhase(readLabel(pimple.dict().lookup("nPhaseCorrectors")));
            for (label corrPhase = 1; corrPhase <= nCorrPhase; ++corrPhase)
            {
                #include "phaseEqn.H"
            }

            Info<< "Solid fraction: min = " << gMin(phase)
                << " max = 1 + " << gMax(phase) - 1 << endl;

            if (gMin(phase) < -ROOTSMALL || gMax(phase) > 1 + ROOTSMALL)
            {
                FatalError
                    << "Phase is out of bounds."
                    << abort(FatalError);
            }

            phase.clip(0, 1);

            if (!pimple.frozenFlow())
            {
                #include "UEqn.H"

                while (pimple.correct())
                {
                    #include "pEqn.H"
                }
            }

            #include "CEqn.H"
        }

        #include "grainEqn.H"

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
