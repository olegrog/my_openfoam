const dictionary& stressControl = mesh.solutionDict().subDict("stressAnalysis");
const label nCorr = stressControl.getOrDefault<label>("nCorrectors", 1000);
const scalar convergenceTolerance = stressControl.get<scalar>("D");
const bool accumulateResidual = stressControl.getOrDefault<bool>("accumulateResidual", true);
const bool accumulateInterlayer = stressControl.getOrDefault<bool>("accumulateInterlayer", true);
const label nResidualCorr =
    accumulateResidual
  ? stressControl.getOrDefault<label>("nResidualCorrectors", 1)
  : 1;

const label totalIter = nResidualCorr*nLayer;

if (runTime.endTime() > totalIter*runTime.deltaT())
{
    runTime.setEndTime(totalIter*runTime.deltaT());
}
