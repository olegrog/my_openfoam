Info<< "isoAdvection: time consumption = "
    << label(100*advector.advectionTime()/runTime.elapsedCpuTime()) << "%" << endl;

Info<< laserHeatSource->type() << "HeatSource: time consumption = "
    << label(100*laserHeatSource->elapsedTime()/runTime.elapsedCpuTime()) << "%" << endl;
