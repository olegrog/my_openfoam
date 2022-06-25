const dictionary& stressControl = mesh.solutionDict().subDict("stressAnalysis");
const label nCorr = stressControl.getOrDefault<label>("nCorrectors", 1000);
const label nResidualCorr = stressControl.getOrDefault<label>("nResidualCorrectors", 1);
const scalar convergenceTolerance = stressControl.get<scalar>("D");
const label totalIter = nResidualCorr*nLayer;
