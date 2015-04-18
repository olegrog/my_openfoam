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
#include "wallFvPatch.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeWallDerivativeField(
    const std::string name,
    const volScalarField& field,
    const Foam::fvMesh& mesh,
    const Foam::Time& runTime
)
{
    dimensionedScalar zeroField(
        "zero",
        field.dimensions(),
        0
    );
    volScalarField wallField(
        IOobject("wallD" + name, runTime.timeName(), mesh),
        mesh,
        zeroField
    );
    const surfaceScalarField derivative = fvc::snGrad(field);
    forAll(wallField.boundaryField(), patchi) {
        wallField.boundaryField()[patchi] = derivative.boundaryField()[patchi];
    }
    wallField.write();
}

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
        // find U0
        solve
        (
            fvm::laplacian(pow(T0,s), U0)
        );
        // find T0
        solve
        (
            gamma1*pow(T0,s)*magSqr(fvc::grad(U0)) + 1.25*gamma2*fvm::laplacian(pow(T0,s), T0)
        );
        // find U1
        solve
        (
            fvm::laplacian(pow(T0,s), U1) + fvc::laplacian(s*pow(T0,s-1)*T1, U0)
        );
        // find T1
        solve
        (
            gamma1*( fvc::grad(U0) && (2*pow(T0,s)*fvc::grad(U1) + s*pow(T0,s-1)*T1*fvc::grad(U0)) )
            + 1.25*gamma2*( fvm::laplacian(pow(T0,s), T1) + fvc::laplacian(s*pow(T0,s-1)*T1, T0) )
        );

        press0 = sum(mesh.V()) / fvc::domainIntegrate(1/T0);
        press1 = sqr(press0) * fvc::domainIntegrate(T1/sqr(T0)) / sum(mesh.V());

        volVectorField P_xy1 = -gamma1*pow(T0,s)*fvc::grad(U0);
        P1 = P1 * 0;
        P1.internalField().replace(symmTensor::XY, P_xy1.internalField().component(vector::Y));
        P1.internalField().replace(symmTensor::XX, 1);
        P1.internalField().replace(symmTensor::YY, 1);
        P1.internalField().replace(symmTensor::ZZ, 1);
        
        q1 = -1.25*gamma2*pow(T0,s)*fvc::grad(T0);

        volScalarField DDT = tr(fvc::grad(fvc::grad(T0)));
        volScalarField DDU = tr(fvc::grad(fvc::grad(U0)));
        volScalarField DT2 = magSqr(fvc::grad(T0));
        volScalarField DU2 = magSqr(fvc::grad(U0));

        volVectorField P_xy2 = -gamma1*( pow(T0,s)*fvc::grad(U1) + s*pow(T0,s-1)*T1*fvc::grad(U0) );
        volScalarField P_xx = (2*(gamma8+gamma9)*T0*DU2 - gamma3*T0*DDT - gamma7*DT2) / 3 / press0;
        volScalarField P_yy = 2*((gamma8-2*gamma9)*T0*DU2 + gamma3*T0*DDT + gamma7*DT2) / 3 / press0;
        volScalarField P_zz = (2*(gamma9-2*gamma8)*T0*DU2 - gamma3*T0*DDT - gamma7*DT2) / 3 / press0;
        volScalarField q_x = (0.5*gamma3 * sqr(T0) * DDU + 4*gamma10*T0*( fvc::grad(T0) & fvc::grad(U0) )) / press0;
        
        P2 = P2 * 0;
        P2.internalField().replace(symmTensor::XY, P_xy2.internalField().component(vector::Y));
        P2.internalField().replace(symmTensor::XX, P_xx.internalField());
        P2.internalField().replace(symmTensor::YY, P_yy.internalField());
        P2.internalField().replace(symmTensor::ZZ, P_zz.internalField());

        q2 = -1.25*gamma2*( pow(T0,s)*fvc::grad(T1) + s*pow(T0,s-1)*T1*fvc::grad(T0) );
        q2.internalField().replace(vector::X, q_x.internalField());

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        runTime.write();
    }
    writeWallDerivativeField("U0", U0, mesh, runTime);
    writeWallDerivativeField("T0", T0, mesh, runTime);

    dimensionedScalar k = kn*Foam::sqrt(constant::mathematical::pi)/2;
    Info<< "Overall integration:" << endl
        << "\tM = " << (fvc::domainIntegrate(U0) + k*fvc::domainIntegrate(U1)).value() << endl
        << "\tT = " << (fvc::domainIntegrate(T0) + k*fvc::domainIntegrate(T1)).value() << endl
        << "\tp = " << (fvc::domainIntegrate(press0) + k*fvc::domainIntegrate(press1)).value() << endl
        << "\tP_ij = " << (k*fvc::domainIntegrate(P1) + k*k*fvc::domainIntegrate(P2)).value() << endl
        << "\tq_i = " << (k*fvc::domainIntegrate(q1) + k*k*fvc::domainIntegrate(q2)).value() << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
