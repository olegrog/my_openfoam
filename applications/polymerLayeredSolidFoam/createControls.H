const dictionary& stressControl = mesh.solutionDict().subDict("stressAnalysis");
const label nCorr = stressControl.getOrDefault<label>("nCorrectors", 1000);
const scalar convergenceTolerance = stressControl.get<scalar>("D");
const bool accumulateInterlayer = stressControl.getOrDefault<bool>("accumulateInterlayer", true);

const label totalIter = accumulateInterlayer ? nLayer : 1;

if (runTime.endTime() > totalIter*runTime.deltaT())
{
    runTime.setEndTime(totalIter*runTime.deltaT());
}
