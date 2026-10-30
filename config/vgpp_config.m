function params = vgpp_config(repoRoot)
%VGPP_CONFIG Validated fixed-final-radius benchmark settings (MATLAB R2024b).
% Settings derived from code_refined(8).m; the reference sweep uses 20 seeds
% at each of six nominal area fractions for both VGPP and RSA-PBC.

params.Lx = 10.0; params.Ly = 10.0;
params.Nx = 250; params.Ny = 250;
params.particleMode = "monodisperse";
params.rMono = 0.18;
params.rMean = 0.18; params.rStd = 0.035;
params.rMin = 0.10; params.rMax = 0.28;
params.minGap = 0;
params.batchSize = 80;
params.nCandidates = 4000;
params.distanceUpdateInterval = 5;
params.useWeightedSampling = true;
% Keep legacy sampling for reproducibility. Switch only for a new validated study.
params.weightedSamplingAlgorithm = 'legacy';
params.useLargestVoidInsertion = true;
params.largestVoidFraction = 0.70;
params.largestVoidMinSeparation = 0.75 * params.rMono;
params.maxCycles = 10000;
params.maxFailedCycles = 500;
params.useRelaxation = true;
params.relaxationCycles = 20;
params.relaxationMovesPerCycle = 5000;
params.initialMoveAmplitude = 0.15 * params.rMono;
params.finalMoveAmplitude = 0.02 * params.rMono;
params.relaxationInsertionAttempts = 10;
params.maxRelaxationFailedCycles = 5;
params.fixedRadiusSeedVF = 0.50;
params.fixedRadiusGrowthScale = 1.0025;
params.fixedRadiusMaxGrowthSteps = 300;
params.fixedRadiusMaxBacktracks = 14;
params.fixedRadiusMinGrowthScale = 1.00002;
params.overlapAuditTolerance = 1.0e-8;
params.growthCompressionOverlapTolerance = params.overlapAuditTolerance;
params.growthCompressionRelaxIterations = 800;
params.growthCompressionRelaxFactor = 0.45;
params.useGrowthCompression = true;
params.useHexagonalFallback = false;
params.stopAtSpecifiedTarget = true;
params.targetVFTolerance = 1.0e-6;
params.rsaMaxTrials = 5000000;
params.rsaMaxFailedTrials = 200000;
params.matrixRichBins = 8;
params.matrixRichLowVFThresholdFactor = 0.75;
params.targetVFSweepList = [0.40 0.50 0.55 0.60 0.65 0.70];
params.targetVFSweepRuns = 20;
params.runRSAInTargetSweep = true;
params.stopIfSpecifiedVFNotReached = false;
params.resultsFolder = fullfile(repoRoot,'outputs','benchmark');
params.figuresFolder = fullfile(repoRoot,'outputs','figures');
params.referenceFile = fullfile(repoRoot,'inputs','reference','Table_R1_fixedRadius_raw.csv');
params.validationTolerance = 1.0e-7;
params.enableCheckpoints = true;
params.resumeCheckpoints = true;
params.checkpointFolder = fullfile(repoRoot,'outputs','checkpoints');
end
