{
    const int oldPrecision = Info.stream().precision();
    const int newPrecision = 2;
    Info<< setprecision(newPrecision);

    if (isA<dynamicRefineFvMesh>(mesh))
    {
        Info<< mesh.type() << ": time consumption = "
            << 100*meshUpdateTime/runTime.elapsedCpuTime() << "%" << endl;
    }

    Info<< "isoAdvection: time consumption = "
        << 100*advector.advectionTime()/runTime.elapsedCpuTime() << "%" << endl;

    Info<< laserHeatSource->type() << "HeatSource: time consumption = "
        << 100*laserHeatSource->elapsedTime()/runTime.elapsedCpuTime() << "%" << endl;

    Info<< setprecision(oldPrecision);
}
