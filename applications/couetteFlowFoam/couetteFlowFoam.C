/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.2.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
-------------------------------------------------------------------------------
Application
    couetteFlowFoam

Description
    The asymptotic equations for the planar Couette problem.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    simpleControl simple(mesh);

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        solve
        (
            fvm::laplacian(pow(T0,s), U0)
        );

        solve
        (
            fvm::laplacian(1.25*gamma2*pow(T0,s), T0) + gamma1*pow(T0,s)*magSqr(fvc::grad(U0))
        );

        volVectorField p_xy = -gamma1*pow(T0,s)*fvc::grad(U0);
        p1 = p1 * 0;
        p1.internalField().replace(symmTensor::XY, p_xy.internalField().component(vector::Y));
        p1.internalField().replace(symmTensor::XX, 1);
        p1.internalField().replace(symmTensor::YY, 1);
        p1.internalField().replace(symmTensor::ZZ, 1);
        
        q1 = -1.25*gamma2*pow(T0,s)*fvc::grad(T0);

        volScalarField DDT = tr(fvc::grad(fvc::grad(T0)));
        volScalarField DDU = tr(fvc::grad(fvc::grad(U0)));
        volScalarField DT2 = magSqr(fvc::grad(T0));
        volScalarField DU2 = magSqr(fvc::grad(U0));

        volScalarField p_xx = (2*(gamma8+gamma9)*T0*DU2 - gamma3*T0*DDT - gamma7*DT2) / 3;
        volScalarField p_yy = 2*((gamma8-2*gamma9)*T0*DU2 + gamma3*T0*DDT + gamma7*DT2) / 3;
        volScalarField p_zz = (2*(gamma9-2*gamma8)*T0*DU2 - gamma3*T0*DDT - gamma7*DT2) / 3;
        volScalarField q_x = 0.5*gamma3 * sqr(T0) * DDU + 4*gamma10*T0*( fvc::grad(T0) & fvc::grad(U0) );
        
        p2 = p2 * 0;
        p2.internalField().replace(symmTensor::XX, p_xx.internalField());
        p2.internalField().replace(symmTensor::YY, p_yy.internalField());
        p2.internalField().replace(symmTensor::ZZ, p_zz.internalField());

        q2 = q2 * 0;
        q2.internalField().replace(vector::X, q_x.internalField());

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
