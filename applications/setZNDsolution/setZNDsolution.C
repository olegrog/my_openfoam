/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2022 Oleg Rogozin
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

Application
    setZNDsolution

Description
    Set initial conditions as the ZND solution for detonation problems.
    Works with dynamicRefineFvMesh.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicRefineFvMesh.H"
#include "psiThermo.H"
#include "ODESolver.H"
#include "interpolationTable.H"

#include "updateGeometricField.H"

#include "ZNDsystem/ZNDsystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"

    // 1. Find the ZND solution
    scalarField y(1, 0);  // lambda
    scalar x = 0;
    scalar dx = odeSolver->relTol()[0];

    using pair = Tuple2<scalar, scalar>;
    DynamicList<pair> lambda1, U1, p1, rho1, T1;

    auto append = [&](scalar x, scalar y)
    {
        lambda1.append(pair(x, y));
        U1.append(pair(x, ode.U(y)));
        p1.append(pair(x, ode.p(y)));
        T1.append(pair(x, WbyR*ode.p(y)/ode.rho(y)));
        rho1.append(pair(x, ode.rho(y)));
        Info<< fixed << x << tab << y << tab << ode.U(y) << tab << ode.p(y)
            << tab << ode.rho(y) << tab << T1.rbegin()->second() << endl;
    };

    Info().precision(6);
    Info<< "\nx\t\tlambda\t\tU\t\tp\t\trho\t\tT" << endl;
    append(x, y[0]);
    do
    {
        odeSolver->solve(x, y, dx);
        append(x, y[0]);
    }
    while (1 - y[0] > odeSolver->relTol()[0]);
    append(x + dx, 1.);

    // 2. Create interpolants
    using interp = interpolationTable<scalar>;
    const interp lambda1d(lambda1, bounds::repeatableBounding::CLAMP, "");
    const interp U1d(U1, bounds::repeatableBounding::CLAMP, "");
    const interp p1d(p1, bounds::repeatableBounding::CLAMP, "");
    const interp T1d(T1, bounds::repeatableBounding::CLAMP, "");
    const interp rho1d(rho1, bounds::repeatableBounding::CLAMP, "");

    // 3. Find fields from the interpolants
    // NB: mesh.update() works only for timeIndex > 0 && timeIndex % refineInterval == 0
    runTime.setTime(0, timeIndex);

    auto update = [](volScalarField& f, const volScalarField& x, const interp& interp1d)
    {
        updateGeometricField
        (
            f, [&](scalar& f, scalar x) { if (x > 0) f = interp1d(x); }, x
        );
        f.correctBoundaryConditions();
    };

    label prevMeshSize;
    Info<< nl << "Initial Mesh size = " << mesh.cells().size() << endl;

    do
    {
        const volScalarField x = -(mesh.C() - origin) & nHat;

        update(lambda, x, lambda1d);
        update(p, x, p1d);
        update(T, x, T1d);
        update(rho, x, rho1d);
        updateGeometricField
        (
            U, [&](vector& U, scalar x) { if (x > 0) U = U1d(x)*nHat; }, x
        );
        U.correctBoundaryConditions();

        if (mesh.dynamic())
        {
            volScalarField& normalisedGradRho = *normalisedGradRhoPtr;
            normalisedGradRho = mag(fvc::grad(rho));
            normalisedGradRho /= gMax(normalisedGradRho);
        }

        prevMeshSize = mesh.cells().size();
        mesh.update();
        Info<< "Mesh size = " << mesh.cells().size() << endl;
    }
    while (mesh.changing() && mesh.cells().size() > prevMeshSize);

    // 4. Save the result
    Info<< nl << "Writing fields" << endl;
    ISstream::defaultPrecision(18);
    runTime.setTime(0, 0);
    runTime.writeNow();

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
