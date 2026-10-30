%% VGPP: Multi-Target Periodic RVE Validation Study
% -------------------------------------------------------------------------
% This code generates a 2D periodic particle distribution using a
% largest-void-guided relaxed distance-transform parallel adsorption
% strategy.
%
% IMPORTANT:
%   Revision-oriented version for manuscript validation.
%   It performs target-VF and multiple-seed comparisons between VGPP and
%   RSA-PBC using the same physical particle settings.
%
% Main outputs:
%   1. Multi-target, multi-seed raw data and summary tables
%   2. Spatial-uniformity, nearest-neighbor and clearance statistics
%   3. Runtime and target-reachability statistics
%   4. Final manuscript figures saved in the folder 'pics':
%      Fig. 3 target reachability
%      Fig. 4 spatial heterogeneity
%      Fig. 5 periodic overlap audit
%      Fig. 6 ordering and runtime
%      Fig. 7 representative fixed-radius configurations
%
% Main features:
%   1. No interphase effect
%   2. Periodic boundary condition in x and y directions
%   3. Particles crossing one edge reappear from the opposite edge
%   4. Largest-void-guided candidate selection
%   5. Optional stochastic relaxation and reinsertion
%   6. Optional RSA-PBC comparison without parameter optimization
%   7. Optional growth-compression stage to push the packing toward a
%      user-defined high VF such as 0.82
%
% Method name:
%   VGPP
%   Void-Guided Particle Packing
%
% Required MATLAB toolbox:
%   Image Processing Toolbox for bwdist() and imregionalmax()
% -------------------------------------------------------------------------

clc;
clear;
close all;

% Output folders for the revised manuscript study
picsFolder = 'pics\R2_C35';
resultsFolder = 'revision_results\R2_C35';
if ~exist(picsFolder,'dir'), mkdir(picsFolder); end
if ~exist(resultsFolder,'dir'), mkdir(resultsFolder); end

% Unified publication style for all generated figures.
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',14);
set(groot,'defaultAxesFontWeight','bold');
set(groot,'defaultAxesLineWidth',2);
set(groot,'defaultLineLineWidth',2);
set(groot,'defaultTextFontSize',14);
set(groot,'defaultTextFontWeight','bold');

%% ------------------------------------------------------------------------
%                           USER INPUTS
% -------------------------------------------------------------------------

% Domain size
Lx = 10.0;
Ly = 10.0;

% Grid resolution
Nx = 250;
Ny = 250;

% Target area fraction for reference only
targetAreaFraction = 0.80;

% Particle-size mode
particleMode = "monodisperse";   % "monodisperse" or "polydisperse"

% Monodisperse particle radius
rMono = 0.18;

% Polydisperse particle radius parameters
rMean = 0.18;
rStd  = 0.035;
rMin  = 0.10;
rMax  = 0.28;

% Minimum surface-to-surface gap between particles
minGap = 0.00;

% Parallel insertion parameters
batchSize = 80;

% Number of candidate locations examined per cycle
nCandidates = 4000;

% Distance-transform update frequency
% 1 = update every cycle
% 5 = update once every 5 cycles
distanceUpdateInterval = 5;

% Candidate selection
useWeightedSampling = true;

% -------------------------------------------------------------------------
% Largest-void-guided insertion parameters
% -------------------------------------------------------------------------

useLargestVoidInsertion = true;
largestVoidFraction = 0.70;
largestVoidMinSeparation = 0.75 * rMono;

% -------------------------------------------------------------------------
% Initial adsorption stopping criteria
% -------------------------------------------------------------------------

maxCycles = 10000;
maxFailedCycles = 500;

% -------------------------------------------------------------------------
% Relaxation and reinsertion parameters
% -------------------------------------------------------------------------

useRelaxation = true;
relaxationCycles = 20;
relaxationMovesPerCycle = 5000;
initialMoveAmplitude = 0.15 * rMono;
finalMoveAmplitude   = 0.02 * rMono;
relaxationInsertionAttempts = 10;
maxRelaxationFailedCycles = 5;

% -------------------------------------------------------------------------
% Growth-compression VF maximization parameters
% -------------------------------------------------------------------------
% This stage is designed to overcome RSA-type jamming. It gradually grows
% the already generated particles and resolves the induced overlaps by
% collective periodic repulsive relaxation. This can push a 2D circle
% packing beyond the usual RSA-like limit toward dense packing values such
% as 0.82.
%
% IMPORTANT: In monodisperse mode, this step increases the final particle
% radius above rMono. If you need the radius to remain exactly rMono, set
% useGrowthCompression = false and use a hexagonal/randomized dense
% initialization instead.

useGrowthCompression = true;
growthCompressionTargetVF = 0.80;
growthCompressionMaxSteps = 350;
growthCompressionStepScale = 1.0025;
growthCompressionRelaxIterations = 800;
growthCompressionRelaxFactor = 0.45;
growthCompressionOverlapTolerance = 1.0e-5;
growthCompressionMinStepScale = 1.00005;
growthCompressionBacktrackFactor = 0.50;
growthCompressionMaxBacktracks = 12;

% -------------------------------------------------------------------------
% Early-stop and final-report parameters
% -------------------------------------------------------------------------
% If the requested VF is not reached after the allowed growth-compression
% iterations or after repeated no-progress steps, the script prints and saves
% a final summary table and terminates before plotting/RSA comparison.

stopIfSpecifiedVFNotReached = true;
growthCompressionNoProgressLimit = 30;
growthCompressionMinVFImprovement = 1.0e-5;

% -------------------------------------------------------------------------
% Growth-compression convergence/stagnation trigger
% -------------------------------------------------------------------------
% If the VF increase over the last N accepted/attempted GC steps becomes
% smaller than the selected tolerance, the GC stage is considered converged
% or practically jammed. The code then exits GC and proceeds to the next
% stage, usually the hexagonal fallback if the target VF has not been reached.

growthCompressionConvergenceWindow = 5;       % recent GC steps used for convergence check
growthCompressionConvergenceTolerance = 1.0e-2; % minimum VF gain over the window
growthCompressionMinVFForConvergenceCheck = 0.70; % avoid premature switching at low VF

% -------------------------------------------------------------------------
% Hexagonal-seeded dense fallback parameters
% -------------------------------------------------------------------------
% If the stochastic growth-compression stage becomes jammed below the
% requested VF, this fallback constructs a periodic hexagonal-seeded dense
% configuration at the requested VF. This is the robust route to 0.82 VF.
% The previous stochastic configuration is replaced only if it cannot reach
% the specified VF.

useHexagonalFallback = false;
hexagonalFallbackTargetVF = growthCompressionTargetVF;
hexagonalFallbackJitterFraction = 0.00;   % 0.00 guarantees the cleanest high-VF packing
hexagonalFallbackRelaxIterations = 5000;
hexagonalFallbackRelaxFactor = 0.45;

% -------------------------------------------------------------------------
% Optional RSA comparison settings
% -------------------------------------------------------------------------
% The RSA comparison uses the same physical settings only:
% domain size, particle-size model, minimum gap, and PBC.
% It does not use LV-R-DT-PAA algorithmic parameters.

runRSAComparison = false;
rsaMaxTrials = 5000000;
rsaMaxFailedTrials = 200000;

% Plot largest residual void centers after final generation
plotLargestVoidCenters = false;
nLargestVoidsToShow = 30;

% -------------------------------------------------------------------------
% Packing-evolution visualization randomness parameters
% -------------------------------------------------------------------------
% These parameters affect only the reconstructed evolution figure. They do
% not modify the final accepted particle coordinates or the overlap audit.
% Larger evolutionRandomnessStrength gives a more natural stochastic-looking
% reveal sequence, while the blue-noise term avoids clustered partial maps.

evolutionRandomSeed = 2026;
evolutionRandomnessStrength = 0.75;      % 0 = very uniform, 1 = more random
evolutionBlueNoiseStrength = 0.65;       % promotes homogeneous stochastic reveal
evolutionCandidatePoolSize = 90;         % candidates tested at each reveal step
evolutionGridBins = 10;                  % spatial balance bins for reveal ordering

% -------------------------------------------------------------------------
% Literature-inspired additional result blocks
% -------------------------------------------------------------------------
% These blocks are optional. The Voronoi/Delaunay and local-VF window plots
% are lightweight and useful for the current final VGPP/RSA result. The VF
% sweep and repeatability study require many independent simulations, so they
% are disabled by default and can be activated for the final manuscript data.

runVoronoiDelaunayBlock = false;
runLocalVFWindowBlock = false;
localVFWindowBinsList = [3 4 5 6 8 10 12 15];

runLegacySingleCase = false;  % legacy single-case analysis is not used in the revised manuscript
runTargetVFSweepBlock = true;   % main revised-manuscript study
targetVFSweepList = [0.40 0.50 0.55 0.60 0.65 0.70];
targetVFSweepRuns = 20;

% Set true only to check that the script runs correctly before the final study.
% In quick-test mode, only 2 seeds are used per target.
quickTestMode = false;   % final manuscript study: 20 seeds per target
if quickTestMode
    targetVFSweepRuns = 2;
end
runRSAInTargetSweep = true;

runRepeatabilityBlock = false;  % redundant because sweep already uses 20 seeds/target
repeatabilityRuns = 10;
repeatabilityTargetVF = growthCompressionTargetVF;
runRSAInRepeatability = true;

% -------------------------------------------------------------------------
% Final selected supportive manuscript blocks
% -------------------------------------------------------------------------
% Optional legacy/supportive blocks retained for compatibility.
% The revised manuscript reports ordering metrics explicitly (NND CV and RDF)
% rather than suppressing regularity-related information.

runSelectedMatrixRichZoneBlock = false;
runSelectedClearanceSafetyBlock = false;
runSelectedStageContributionBlock = false;

% Heavy multi-run blocks. Activate only when preparing final manuscript data.
runSelectedTargetReachabilityBlock = false;
runSelectedRepeatabilitySupportBlock = false;

matrixRichBins = 8;
matrixRichLowVFThresholdFactor = 0.75;   % low-local-VF cells below 0.75*mean are highlighted
nearContactThresholdFactor = 0.02;       % near-contact threshold = 0.02*rMono

% Random seed for repeatability
rng(10);

% During the multi-target study, insertion is capped at the requested VF.
% Growth-compression may then close the very small residual VF gap.
stopAtSpecifiedTarget = true;
targetVFTolerance = 1.0e-6;

% Fixed-final-radius densification settings for the revised study.
% All accepted final VGPP particles have radius exactly rMono.
fixedRadiusSeedVF = 0.50;          % initial shrunken-particle packing level
fixedRadiusGrowthScale = 1.0025;  % gradual radius growth toward rMono
fixedRadiusMaxGrowthSteps = 300;
fixedRadiusMaxBacktracks = 14;
fixedRadiusMinGrowthScale = 1.00002;

% One tolerance is used by both overlap resolution and final audit.
overlapAuditTolerance = 1.0e-8;

%% ------------------------------------------------------------------------
%                           PARAMETER STRUCTURE
% -------------------------------------------------------------------------

params.Lx = Lx;
params.Ly = Ly;
params.Nx = Nx;
params.Ny = Ny;

params.targetAreaFraction = targetAreaFraction;

params.particleMode = particleMode;
params.rMono = rMono;
params.rMean = rMean;
params.rStd = rStd;
params.rMin = rMin;
params.rMax = rMax;

params.minGap = minGap;

params.batchSize = batchSize;
params.nCandidates = nCandidates;
params.useWeightedSampling = useWeightedSampling;

params.useLargestVoidInsertion = useLargestVoidInsertion;
params.largestVoidFraction = largestVoidFraction;
params.largestVoidMinSeparation = largestVoidMinSeparation;

params.distanceUpdateInterval = distanceUpdateInterval;

params.maxCycles = maxCycles;
params.maxFailedCycles = maxFailedCycles;

params.useRelaxation = useRelaxation;
params.relaxationCycles = relaxationCycles;
params.relaxationMovesPerCycle = relaxationMovesPerCycle;
params.initialMoveAmplitude = initialMoveAmplitude;
params.finalMoveAmplitude = finalMoveAmplitude;
params.relaxationInsertionAttempts = relaxationInsertionAttempts;
params.maxRelaxationFailedCycles = maxRelaxationFailedCycles;

params.useGrowthCompression = useGrowthCompression;
params.growthCompressionTargetVF = growthCompressionTargetVF;
params.growthCompressionMaxSteps = growthCompressionMaxSteps;
params.growthCompressionStepScale = growthCompressionStepScale;
params.growthCompressionRelaxIterations = growthCompressionRelaxIterations;
params.growthCompressionRelaxFactor = growthCompressionRelaxFactor;
params.growthCompressionOverlapTolerance = growthCompressionOverlapTolerance;
params.growthCompressionMinStepScale = growthCompressionMinStepScale;
params.growthCompressionBacktrackFactor = growthCompressionBacktrackFactor;
params.growthCompressionMaxBacktracks = growthCompressionMaxBacktracks;
params.stopIfSpecifiedVFNotReached = stopIfSpecifiedVFNotReached;
params.growthCompressionNoProgressLimit = growthCompressionNoProgressLimit;
params.growthCompressionMinVFImprovement = growthCompressionMinVFImprovement;
params.growthCompressionConvergenceWindow = growthCompressionConvergenceWindow;
params.growthCompressionConvergenceTolerance = growthCompressionConvergenceTolerance;
params.growthCompressionMinVFForConvergenceCheck = growthCompressionMinVFForConvergenceCheck;

params.useHexagonalFallback = useHexagonalFallback;
params.hexagonalFallbackTargetVF = hexagonalFallbackTargetVF;
params.hexagonalFallbackJitterFraction = hexagonalFallbackJitterFraction;
params.hexagonalFallbackRelaxIterations = hexagonalFallbackRelaxIterations;
params.hexagonalFallbackRelaxFactor = hexagonalFallbackRelaxFactor;

params.evolutionRandomSeed = evolutionRandomSeed;
params.evolutionRandomnessStrength = evolutionRandomnessStrength;
params.evolutionBlueNoiseStrength = evolutionBlueNoiseStrength;
params.evolutionCandidatePoolSize = evolutionCandidatePoolSize;
params.evolutionGridBins = evolutionGridBins;

params.runVoronoiDelaunayBlock = runVoronoiDelaunayBlock;
params.runLocalVFWindowBlock = runLocalVFWindowBlock;
params.localVFWindowBinsList = localVFWindowBinsList;
params.runTargetVFSweepBlock = runTargetVFSweepBlock;
params.targetVFSweepList = targetVFSweepList;
params.targetVFSweepRuns = targetVFSweepRuns;
params.runRSAInTargetSweep = runRSAInTargetSweep;
params.runRepeatabilityBlock = runRepeatabilityBlock;
params.repeatabilityRuns = repeatabilityRuns;
params.repeatabilityTargetVF = repeatabilityTargetVF;
params.runRSAInRepeatability = runRSAInRepeatability;

params.runSelectedMatrixRichZoneBlock = runSelectedMatrixRichZoneBlock;
params.runSelectedClearanceSafetyBlock = runSelectedClearanceSafetyBlock;
params.runSelectedStageContributionBlock = runSelectedStageContributionBlock;
params.runSelectedTargetReachabilityBlock = runSelectedTargetReachabilityBlock;
params.runSelectedRepeatabilitySupportBlock = runSelectedRepeatabilitySupportBlock;
params.matrixRichBins = matrixRichBins;
params.matrixRichLowVFThresholdFactor = matrixRichLowVFThresholdFactor;
params.nearContactThresholdFactor = nearContactThresholdFactor;

params.stopAtSpecifiedTarget = stopAtSpecifiedTarget;
params.targetVFTolerance = targetVFTolerance;
params.rsaMaxTrials = rsaMaxTrials;
params.rsaMaxFailedTrials = rsaMaxFailedTrials;

params.fixedRadiusSeedVF = fixedRadiusSeedVF;
params.fixedRadiusGrowthScale = fixedRadiusGrowthScale;
params.fixedRadiusMaxGrowthSteps = fixedRadiusMaxGrowthSteps;
params.fixedRadiusMaxBacktracks = fixedRadiusMaxBacktracks;
params.fixedRadiusMinGrowthScale = fixedRadiusMinGrowthScale;
params.overlapAuditTolerance = overlapAuditTolerance;

% Use the same strict tolerance in repulsive relaxation.
params.growthCompressionOverlapTolerance = overlapAuditTolerance;


%% R2 C35: independent robustness / confidence / one-factor sensitivity
% Full scientific audit uses original run_fixed_final_radius_VGPP() logic.
% Grid spacing uses L/(Nx-1). Changing L with r fixed alters N and r/L;
% changing r with L fixed likewise changes N. Interpret effects accordingly.
% MATLAB R2024b + Image Processing Toolbox (bwdist, imregionalmax).
% Exact Clopper-Pearson intervals use betaincinv (Statistics Toolbox).

quickTestMode = false;     % first run: 2 scenarios x one target x one seed
nSeedsFull = 10;          % targeted sensitivity, NOT proof of high reliability
nSeedsQuick = 1;
baseSeed = 35000;
targetsFull = [0.50,0.70];
if quickTestMode
    targets = 0.50;
    nSeeds = nSeedsQuick;
else
    targets = targetsFull;
    nSeeds = nSeedsFull;
end

% One-factor changes relative to the 10x10, r=0.18, 250x250 baseline.
caseName = ["baseline";"larger_domain";"smaller_domain"; ...
            "larger_radius";"smaller_radius";"finer_grid";"coarser_grid"];
caseL = [10; 15; 7.5; 10; 10; 10; 10];
caseR = [0.18; 0.18; 0.18; 0.24; 0.14; 0.18; 0.18];
caseGrid = [250;375;188;250;250;375;188];
% `larger_domain` and `smaller_domain` preserve baseline approximate h,
% while finer/coarser-grid scenarios deliberately change h at fixed L.
if quickTestMode
    activeCases = [1 6];
else
    activeCases = 1:numel(caseName);
end

outFolder = resultsFolder;
if ~exist(outFolder,'dir'), mkdir(outFolder); end
picsFolder35 = picsFolder;
if ~exist(picsFolder35,'dir'), mkdir(picsFolder35); end

% Portable checkpoint: each finished run is appended to disk immediately.
%% R2 C35: RESUME FROM PREVIOUS CHECKPOINT

rawName = fullfile(outFolder, 'R2_C35_raw.csv');

if ~isfile(rawName)
    error(['Checkpoint CSV not found: %s\n' ...
           'Copy the original R2_C35_raw.csv to this location first.'], ...
           rawName);
end

Tprevious = readtable(rawName, 'TextType', 'string');
oldRows = table2struct(Tprevious);

fprintf('\n============================================\n');
fprintf('R2 C35 - RESUMING PREVIOUS SIMULATIONS\n');
fprintf('Previously completed runs: %d\n', height(Tprevious));
fprintf('Expected total runs: %d\n', ...
    numel(activeCases)*numel(targets)*nSeeds);
fprintf('Remaining runs: %d\n', ...
    numel(activeCases)*numel(targets)*nSeeds - height(Tprevious));
fprintf('============================================\n\n');


rows = struct('Scenario',{},'RequestedAF',{},'Run',{},'Seed',{}, ...
    'Lx',{},'Ly',{},'FinalRadius',{},'GridN',{},'GridSpacing',{}, ...
    'RelativeGridSpacing',{},'RadiusDomainRatio',{},'Ntarget',{}, ...
    'Nfinal',{},'ActualTargetAF',{},'FinalAF',{},'TargetReached',{}, ...
    'RadiusCorrect',{},'OverlapAuditPassed',{},'MaxOverlap',{}, ...
    'MinClearance',{},'Runtime_s',{},'GrowthSteps',{},'SeedRadius',{}, ...
    'ErrorMessage',{});
kk = 0;
if ~isempty(oldRows)
    rows = oldRows(:);
    kk = numel(rows);
end
for ia = 1:numel(activeCases)
    ic = activeCases(ia);
    for it = 1:numel(targets)
        nominalAF = targets(it);
        for ir = 1:nSeeds
            p = params;
            p.Lx = caseL(ic); p.Ly = caseL(ic);
            p.rMono = caseR(ic); p.rMean = p.rMono;
            p.largestVoidMinSeparation = 0.75*p.rMono;
            p.initialMoveAmplitude = 0.15*p.rMono;
            p.finalMoveAmplitude = 0.02*p.rMono;
            p.Nx = caseGrid(ic); p.Ny = caseGrid(ic);
            p.targetAreaFraction = nominalAF;
            p.targetParticleCount = max(1,round(nominalAF*p.Lx*p.Ly/(pi*p.rMono^2)));
            % Keep same algorithmic budgets and stopping rules as base code.
            % No lattice/hexagonal fallback allowed.
            p.useHexagonalFallback = false;
            p.overlapAuditTolerance = 1e-8;
            p.growthCompressionOverlapTolerance = 1e-8;
            seedNow = baseSeed + 10000*ic + 100*it + ir;
            if kk>0 && any(strcmp(string({rows.Scenario}),caseName(ic)) & ...
                    [rows.Seed]==seedNow & abs([rows.RequestedAF]-nominalAF)<1e-12)
                fprintf('Skipping completed scenario %s, AF %.2f, seed %d\n',caseName(ic),nominalAF,seedNow);
                continue;
            end
            rng(seedNow,'twister');
            err = "";
            tic;
            try
                res = run_fixed_final_radius_VGPP(p);
                elapsed = toc;
                nfinal = height(res.particleTable);
                radiusOK = nfinal>0 && all(abs(res.particleTable.radius-p.rMono) <= 1e-10);
                [auditedMaxOverlap,auditedMinGap,auditedViolations] = ...
                    C35_coordinate_audit(res.particleTable,p.Lx,p.Ly,p.overlapAuditTolerance);
                auditOK = auditedViolations==0;
                targetOK = logical(res.reachedTarget) && nfinal==p.targetParticleCount && radiusOK && auditOK;
                finalAF = sum(pi*res.particleTable.radius.^2)/(p.Lx*p.Ly);
                maxOver = auditedMaxOverlap;
                minClear = auditedMinGap;
                gSteps = res.growthSteps;
                seedRadius = res.seedRadius;
            catch ME
                elapsed=toc;
                nfinal=NaN; radiusOK=false; auditOK=false; targetOK=false;
                finalAF=NaN; maxOver=NaN; minClear=NaN;
                gSteps=NaN; seedRadius=NaN;
                err=string(ME.identifier)+": "+string(ME.message);
                warning('R2_C35:RunFailed','Scenario %s target %.2f seed %d: %s', ...
                    char(caseName(ic)),nominalAF,seedNow,char(err));
            end
            kk=kk+1;
            rows(kk).Scenario = caseName(ic);
            rows(kk).RequestedAF = nominalAF;
            rows(kk).Run = ir;
            rows(kk).Seed = seedNow;
            rows(kk).Lx = p.Lx;
            rows(kk).Ly = p.Ly;
            rows(kk).FinalRadius = p.rMono;
            rows(kk).GridN = p.Nx;
            rows(kk).GridSpacing = p.Lx/(p.Nx-1);
            rows(kk).RelativeGridSpacing = (p.Lx/(p.Nx-1))/p.rMono;
            rows(kk).RadiusDomainRatio = p.rMono/p.Lx;
            rows(kk).Ntarget = p.targetParticleCount;
            rows(kk).Nfinal = nfinal;
            rows(kk).ActualTargetAF = p.targetParticleCount*pi*p.rMono^2/(p.Lx*p.Ly);
            rows(kk).FinalAF = finalAF;
            rows(kk).TargetReached = double(targetOK);
            rows(kk).RadiusCorrect = double(radiusOK);
            rows(kk).OverlapAuditPassed = double(auditOK);
            rows(kk).MaxOverlap = maxOver;
            rows(kk).MinClearance = minClear;
            rows(kk).Runtime_s = elapsed;
            rows(kk).GrowthSteps = gSteps;
            rows(kk).SeedRadius = seedRadius;
            rows(kk).ErrorMessage = err;
            writetable(struct2table(rows),rawName);
            fprintf('R2 C35 [%s / AF=%.2f] %d/%d success=%d N=%g/%d runtime=%.1fs\n', ...
                char(caseName(ic)),nominalAF,ir,nSeeds,targetOK,nfinal,p.targetParticleCount,elapsed);
        end
    end
end

T = struct2table(rows);
S = table();
for ic = activeCases(:)'
    for it = 1:numel(targets)
        filt = T.Scenario==caseName(ic) & abs(T.RequestedAF-targets(it))<1e-12;
        R=T(filt,:);
        n=height(R); k=sum(R.TargetReached);
        alpha=0.05;
        if k==0, lower=0; else, lower=betaincinv(alpha/2,k,n-k+1); end
        if k==n, upper=1; else, upper=betaincinv(1-alpha/2,k+1,n-k); end
        rr=table(caseName(ic),targets(it),n,k,k/n,lower,upper, ...
            mean(R.Runtime_s),median(R.Runtime_s),std(R.Runtime_s), ...
            mean(R.Nfinal,'omitnan'),R.Ntarget(1),R.GridSpacing(1), ...
            R.RadiusDomainRatio(1),sum(~ismissing(string(R.ErrorMessage)) & strlength(strtrim(string(R.ErrorMessage)))>0), ...
            'VariableNames',{'Scenario','RequestedAF','Runs','Successes', ...
            'SuccessRate','CI95_Lower','CI95_Upper','RuntimeMean_s', ...
            'RuntimeMedian_s','RuntimeSD_s','MeanNfinal','Ntarget', ...
            'GridSpacing','RadiusDomainRatio','RuntimeErrors'});
        S=[S;rr]; %#ok<AGROW>
    end
end
writetable(S,fullfile(outFolder,'R2_C35_summary.csv'));

% Explicit statistical note: the interval is per-condition and assumes
% independent stochastic runs. It does not assess out-of-sample settings.
F=figure('Color','w','Position',[100 100 1200 550]);
hold on;
for i=1:height(S)
    x=i;
    line([x x],[S.CI95_Lower(i) S.CI95_Upper(i)],'Color','k','LineWidth',1.4);
    plot(x,S.SuccessRate(i),'ko','MarkerFaceColor','k','MarkerSize',5);
end
ylim([0 1.03]); xlim([0 height(S)+1]);
set(gca,'XTick',1:height(S),'XTickLabel',strcat(S.Scenario," / ",string(S.RequestedAF)), ...
    'XTickLabelRotation',40,'FontName','Times New Roman');
ylabel('Observed success rate and exact 95% CI','FontName','Times New Roman');
grid on; box on;
exportgraphics(F,fullfile(picsFolder35,'R2_C35_success_confidence.png'),'Resolution',250);
close(F);

F=figure('Color','w','Position',[100 100 1100 500]);
boxplot(T.Runtime_s,cellstr(strcat(T.Scenario," / ",string(T.RequestedAF))));
set(gca,'XTickLabelRotation',40,'FontName','Times New Roman');
ylabel('Runtime (s)','FontName','Times New Roman');
grid on;
exportgraphics(F,fullfile(picsFolder35,'R2_C35_runtime_distributions.png'),'Resolution',250);
close(F);

disp(S);
fprintf('Done. Saved outputs to %s\n',outFolder);

function summaryTable = make_final_summary_table(result, params)

    % Compact final report. This table is displayed and saved even when the
    % code terminates early because the specified VF was not reached.

    if isfield(params, 'useGrowthCompression') && params.useGrowthCompression
        specifiedVF = params.growthCompressionTargetVF;
    else
        specifiedVF = params.targetAreaFraction;
    end

    finalVF = result.maxAreaFraction;
    absoluteError = specifiedVF - finalVF;
    percentOfTarget = 100 * finalVF / specifiedVF;

    if height(result.particleTable) > 0
        meanRadius = mean(result.particleTable.radius);
        minRadius = min(result.particleTable.radius);
        maxRadius = max(result.particleTable.radius);
    else
        meanRadius = NaN;
        minRadius = NaN;
        maxRadius = NaN;
    end

    summaryTable = table( ...
        specifiedVF, ...
        finalVF, ...
        absoluteError, ...
        percentOfTarget, ...
        result.growthCompressionReachedTarget, ...
        string(result.growthCompressionStopReason), ...
        result.hexagonalFallbackUsed, ...
        string(result.hexagonalFallbackStopReason), ...
        height(result.particleTable), ...
        meanRadius, ...
        minRadius, ...
        maxRadius, ...
        result.cycles, ...
        result.failedCycles, ...
        result.relaxationCyclesUsed, ...
        result.relaxationInsertedParticles, ...
        result.growthCompressionCyclesUsed, ...
        result.growthCompressionAcceptedSteps, ...
        'VariableNames', { ...
        'SpecifiedVF', ...
        'FinalVF', ...
        'VFShortfall', ...
        'PercentOfSpecifiedVF', ...
        'ReachedSpecifiedVF', ...
        'StopReason', ...
        'HexagonalFallbackUsed', ...
        'HexagonalFallbackStopReason', ...
        'NumberOfParticles', ...
        'MeanRadius', ...
        'MinRadius', ...
        'MaxRadius', ...
        'InitialInsertionCycles', ...
        'InitialFailedCycles', ...
        'RelaxationCyclesUsed', ...
        'ParticlesInsertedDuringRelaxation', ...
        'GrowthCompressionCyclesUsed', ...
        'GrowthCompressionAcceptedSteps'});
end


function result = run_LV_R_DTPAA_PBC(params)

    %% Domain and grid

    Lx = params.Lx;
    Ly = params.Ly;
    Nx = params.Nx;
    Ny = params.Ny;

    x = linspace(0, Lx, Nx);
    y = linspace(0, Ly, Ny);

    [X, Y] = meshgrid(x, y);

    dxGrid = Lx/(Nx - 1);
    dyGrid = Ly/(Ny - 1);
    pixelSize = mean([dxGrid, dyGrid]);

    domainArea = Lx * Ly;
    if isfield(params,'activeTargetVF') && isfinite(params.activeTargetVF)
        targetVFactive = params.activeTargetVF;
    elseif isfield(params,'growthCompressionTargetVF')
        targetVFactive = params.growthCompressionTargetVF;
    else
        targetVFactive = params.targetAreaFraction;
    end
    targetArea = targetVFactive * domainArea;

    %% Occupancy map

    occCore = false(Ny, Nx);

    %% Particle storage

    cxList = [];
    cyList = [];
    rList  = [];

    %% History storage

    cycleHistory = [];
    stageHistory = strings(0,1);
    particleHistory = [];
    areaFractionHistory = [];

    %% Initialization

    currentArea = 0.0;
    cycles = 0;
    failedCycles = 0;

    targetReached = false;
    targetAreaFractionReached = NaN;

    distOcc = inf(Ny, Nx);

    %% --------------------------------------------------------------------
    %  Stage 1: Initial DT-PAA-PBC adsorption
    % ---------------------------------------------------------------------

    while cycles < params.maxCycles && failedCycles < params.maxFailedCycles

        cycles = cycles + 1;

        updateDistanceMap = cycles == 1 || ...
            mod(cycles - 1, params.distanceUpdateInterval) == 0;

        [cxList, cyList, rList, occCore, distOcc, nInserted, currentArea] = ...
            attempt_DTPAA_batch_insertion(cxList, cyList, rList, occCore, distOcc, ...
            currentArea, X, Y, pixelSize, params, updateDistanceMap);

        if nInserted == 0
            % If the remaining VF gap is smaller than approximately one
            % particle increment, stop insertion and let growth-compression
            % close the small residual gap exactly.
            if isfield(params,'stopAtSpecifiedTarget') && params.stopAtSpecifiedTarget && ...
                    isfield(params,'activeTargetVF') && isfinite(params.activeTargetVF)
                currentVFtmp = currentArea / domainArea;
                if params.particleMode == "monodisperse"
                    oneParticleVF = pi*params.rMono^2/domainArea;
                else
                    oneParticleVF = pi*params.rMin^2/domainArea;
                end
                if params.activeTargetVF - currentVFtmp <= 1.05*oneParticleVF
                    targetReached = true;
                    targetAreaFractionReached = currentVFtmp;
                    fprintf('Insertion stopped near target VF: %.5f (target %.5f).\n', ...
                        currentVFtmp, params.activeTargetVF);
                    break;
                end
            end

            failedCycles = failedCycles + 1;
            continue;
        else
            failedCycles = 0;
        end

        currentAreaFraction = currentArea / domainArea;

        cycleHistory(end+1,1) = cycles;
        stageHistory(end+1,1) = "Initial adsorption";
        particleHistory(end+1,1) = numel(rList);
        areaFractionHistory(end+1,1) = currentAreaFraction;

        if ~targetReached && currentArea >= targetArea
            targetReached = true;
            targetAreaFractionReached = currentAreaFraction;
        elseif isfield(params,'stopAtSpecifiedTarget') && params.stopAtSpecifiedTarget && ...
                isfield(params,'activeTargetVF') && isfinite(params.activeTargetVF)
            if params.particleMode == "monodisperse"
                oneParticleVF = pi*params.rMono^2/domainArea;
            else
                oneParticleVF = pi*params.rMin^2/domainArea;
            end
            if params.activeTargetVF - currentAreaFraction <= 1.05*oneParticleVF
                targetReached = true;
                targetAreaFractionReached = currentAreaFraction;
                break;
            end
        end
    end

    %% --------------------------------------------------------------------
    %  Stage 2: Optional stochastic relaxation and reinsertion
    % ---------------------------------------------------------------------

    relaxationCyclesUsed = 0;
    relaxationInsertedParticles = 0;
    relaxationFailedCycles = 0;

    if isfield(params, 'useRelaxation') && params.useRelaxation && numel(rList) > 1

        fprintf('\nStarting relaxation-reinsertion stage...\n');

        for relaxCycle = 1:params.relaxationCycles

            relaxationCyclesUsed = relaxationCyclesUsed + 1;

            %% Gradually reduce movement amplitude
            alpha = (relaxCycle - 1) / max(params.relaxationCycles - 1, 1);
            moveAmplitude = (1 - alpha) * params.initialMoveAmplitude + ...
                             alpha * params.finalMoveAmplitude;

            %% Stochastic relaxation of existing particles
            [cxList, cyList] = stochastic_relaxation_PBC(cxList, cyList, rList, ...
                params, moveAmplitude, params.relaxationMovesPerCycle);

            %% Rebuild occupancy map after particle motion
            occCore = rebuild_periodic_occupancy(X, Y, cxList, cyList, rList, params);
            distOcc = inf(size(occCore));

            %% Reinsertion after relaxation
            insertedThisRelaxCycle = 0;

            for attempt = 1:params.relaxationInsertionAttempts

                [cxList, cyList, rList, occCore, distOcc, nInserted, currentArea] = ...
                    attempt_DTPAA_batch_insertion(cxList, cyList, rList, occCore, distOcc, ...
                    currentArea, X, Y, pixelSize, params, true);

                insertedThisRelaxCycle = insertedThisRelaxCycle + nInserted;
                relaxationInsertedParticles = relaxationInsertedParticles + nInserted;

                if nInserted == 0
                    break;
                end
            end

            currentAreaFraction = currentArea / domainArea;

            cycleHistory(end+1,1) = cycles + relaxCycle;
            stageHistory(end+1,1) = "Relaxation-reinsertion";
            particleHistory(end+1,1) = numel(rList);
            areaFractionHistory(end+1,1) = currentAreaFraction;

            fprintf('Relaxation cycle %d/%d: inserted = %d, maxVF = %.5f\n', ...
                relaxCycle, params.relaxationCycles, insertedThisRelaxCycle, currentAreaFraction);

            if insertedThisRelaxCycle == 0
                relaxationFailedCycles = relaxationFailedCycles + 1;
            else
                relaxationFailedCycles = 0;
            end

            if relaxationFailedCycles >= params.maxRelaxationFailedCycles
                fprintf('Relaxation stage stopped because no new particles were inserted.\n');
                break;
            end
        end
    end

    %% --------------------------------------------------------------------
    %  Stage 3: Growth-compression VF maximization
    % ---------------------------------------------------------------------

    growthCompressionCyclesUsed = 0;
    growthCompressionAcceptedSteps = 0;
    growthCompressionReachedTarget = false;
    growthCompressionStopReason = "Growth-compression disabled or not attempted";

    hexagonalFallbackUsed = false;
    hexagonalFallbackStopReason = "Hexagonal fallback disabled or not attempted";

    if isfield(params, 'useGrowthCompression') && params.useGrowthCompression && numel(rList) > 1

        fprintf('\nStarting growth-compression VF maximization stage...\n');

        [cxList, cyList, rList, occCore, currentArea, gcInfo, gcHistory] = ...
            growth_compression_PBC(cxList, cyList, rList, X, Y, params, domainArea);

        growthCompressionCyclesUsed = gcInfo.cyclesUsed;
        growthCompressionAcceptedSteps = gcInfo.acceptedSteps;
        growthCompressionReachedTarget = gcInfo.reachedTarget;
        growthCompressionStopReason = gcInfo.stopReason;

        if ~isempty(gcHistory.areaFraction)
            nGC = numel(gcHistory.areaFraction);
            cycleHistory = [cycleHistory; cycles + relaxationCyclesUsed + (1:nGC)'];
            stageHistory = [stageHistory; repmat("Growth-compression", nGC, 1)];
            particleHistory = [particleHistory; repmat(numel(rList), nGC, 1)];
            areaFractionHistory = [areaFractionHistory; gcHistory.areaFraction(:)];
        end
    end

    %% --------------------------------------------------------------------
    %  Stage 4: Hexagonal-seeded dense fallback for specified VF
    % ---------------------------------------------------------------------

    currentVF_afterGC = currentArea / domainArea;

    if isfield(params, 'useHexagonalFallback') && params.useHexagonalFallback && ...
            currentVF_afterGC < params.hexagonalFallbackTargetVF

        fprintf('\nGrowth-compression jammed below the specified VF.\n');
        fprintf('Starting hexagonal-seeded dense fallback stage...\n');

        [cxHex, cyHex, rHex, occHex, hexInfo] = ...
            generate_hexagonal_seeded_dense_PBC(params, X, Y, domainArea);

        cxList = cxHex;
        cyList = cyHex;
        rList  = rHex;
        occCore = occHex;
        currentArea = sum(pi * rList.^2);

        hexagonalFallbackUsed = true;
        hexagonalFallbackStopReason = hexInfo.stopReason;

        currentVF_hex = currentArea / domainArea;

        cycleHistory(end+1,1) = cycles + relaxationCyclesUsed + growthCompressionCyclesUsed + 1;
        stageHistory(end+1,1) = "Hexagonal-seeded dense fallback";
        particleHistory(end+1,1) = numel(rList);
        areaFractionHistory(end+1,1) = currentVF_hex;

        if currentVF_hex >= params.hexagonalFallbackTargetVF
            growthCompressionReachedTarget = true;
            growthCompressionStopReason = "Specified VF reached by hexagonal-seeded dense fallback";
            fprintf('Hexagonal fallback reached VF = %.5f.\n', currentVF_hex);
        else
            fprintf('Hexagonal fallback ended at VF = %.5f before reaching %.5f.\n', ...
                currentVF_hex, params.hexagonalFallbackTargetVF);
        end
    end

    %% Final maximum achieved fraction

    maxAreaFraction = currentArea / domainArea;

    %% Prepare particle table

    particleID = (1:numel(cxList))';

    particleTable = table(particleID, cxList, cyList, rList, ...
        'VariableNames', {'ID', 'x', 'y', 'radius'});

    %% Prepare history table

    historyTable = table(cycleHistory, stageHistory, particleHistory, areaFractionHistory, ...
        'VariableNames', {'Cycle', 'Stage', 'NumberOfParticles', 'AreaFraction'});

    %% If target was not reached

    if ~targetReached
        targetAreaFractionReached = maxAreaFraction;
        warning('Target area fraction was not reached before the stopping criterion.');
    end

    %% Output structure

    result.X = X;
    result.Y = Y;
    result.x = x;
    result.y = y;

    result.occCore = occCore;

    result.particleTable = particleTable;
    result.historyTable = historyTable;

    result.targetAreaFractionReached = targetAreaFractionReached;
    result.maxAreaFraction = maxAreaFraction;

    result.cycles = cycles;
    result.failedCycles = failedCycles;
    result.relaxationCyclesUsed = relaxationCyclesUsed;
    result.relaxationInsertedParticles = relaxationInsertedParticles;

    result.growthCompressionCyclesUsed = growthCompressionCyclesUsed;
    result.growthCompressionAcceptedSteps = growthCompressionAcceptedSteps;
    result.growthCompressionReachedTarget = growthCompressionReachedTarget;
    result.growthCompressionStopReason = growthCompressionStopReason;

    result.hexagonalFallbackUsed = hexagonalFallbackUsed;
    result.hexagonalFallbackStopReason = hexagonalFallbackStopReason;

    result.dxGrid = dxGrid;
    result.dyGrid = dyGrid;
    result.pixelSize = pixelSize;
end

%% ------------------------------------------------------------------------

function [cxList, cyList, rList, occCore, distOcc, nInserted, currentArea] = ...
    attempt_DTPAA_batch_insertion(cxList, cyList, rList, occCore, distOcc, ...
    currentArea, X, Y, pixelSize, params, updateDistanceMap)

    Lx = params.Lx;
    Ly = params.Ly;

    nInserted = 0;

    %% Trial particle radius
    rCurrent = generate_particle_radius(params);

    %% Update distance transform if needed
    if updateDistanceMap
        if any(occCore(:))
            distOcc = periodic_distance_transform(occCore, pixelSize);
        else
            distOcc = inf(size(occCore));
        end
    end

    %% Identify admissible points
    if any(occCore(:))
        admissible = distOcc >= (rCurrent + params.minGap);
    else
        admissible = true(size(occCore));
    end

    admissibleIndices = find(admissible);

    if isempty(admissibleIndices)
        return;
    end

    %% Select candidate points
    nAvailable = numel(admissibleIndices);
    nCandidates = min(nAvailable, params.nCandidates);

    if nCandidates < 1
        return;
    end

    selectedIndices = select_candidate_indices_LV_hybrid( ...
        admissibleIndices, distOcc, nCandidates, params, X, Y);

    candX = X(selectedIndices);
    candY = Y(selectedIndices);

    %% Build one insertion batch
    batchX = [];
    batchY = [];
    batchR = [];

    for k = 1:numel(candX)

        xNew = candX(k);
        yNew = candY(k);
        rNew = rCurrent;

        %% Exact periodic check against existing particles
        if ~isempty(cxList)

            dExisting = periodic_center_distance( ...
                xNew, yNew, cxList, cyList, Lx, Ly);

            requiredDistanceExisting = rList + rNew + params.minGap;

            if any(dExisting < requiredDistanceExisting)
                continue;
            end
        end

        %% Exact periodic check against current batch particles
        if ~isempty(batchX)

            dBatch = periodic_center_distance( ...
                xNew, yNew, batchX, batchY, Lx, Ly);

            requiredDistanceBatch = batchR + rNew + params.minGap;

            if any(dBatch < requiredDistanceBatch)
                continue;
            end
        end

        %% Optional cap at the specified target VF
        if isfield(params,'stopAtSpecifiedTarget') && params.stopAtSpecifiedTarget && ...
                isfield(params,'activeTargetVF') && isfinite(params.activeTargetVF)
            domainAreaLocal = params.Lx * params.Ly;
            targetAreaLocal = params.activeTargetVF * domainAreaLocal;
            batchAreaNow = sum(pi * batchR.^2);
            projectedArea = currentArea + batchAreaNow + pi*rNew^2;
            if projectedArea > targetAreaLocal + params.targetVFTolerance*domainAreaLocal
                continue;
            end
        end

        %% Accept candidate into the batch
        batchX(end+1,1) = xNew;
        batchY(end+1,1) = yNew;
        batchR(end+1,1) = rNew;

        if numel(batchX) >= params.batchSize
            break;
        end
    end

    %% Insert accepted batch into occupancy map
    if isempty(batchX)
        return;
    end

    for p = 1:numel(batchX)

        xp = batchX(p);
        yp = batchY(p);
        rp = batchR(p);

        particleMask = periodic_particle_mask(X, Y, xp, yp, rp, Lx, Ly);
        occCore = occCore | particleMask;

        cxList(end+1,1) = xp;
        cyList(end+1,1) = yp;
        rList(end+1,1)  = rp;

        currentArea = currentArea + pi * rp^2;
        nInserted = nInserted + 1;
    end
end

%% ------------------------------------------------------------------------

function [cxList, cyList] = stochastic_relaxation_PBC(cxList, cyList, rList, params, moveAmplitude, nMoves)

    Lx = params.Lx;
    Ly = params.Ly;
    minGap = params.minGap;

    nParticles = numel(rList);

    if nParticles < 2
        return;
    end

    for moveID = 1:nMoves

        %% Randomly select one particle
        i = randi(nParticles);

        oldX = cxList(i);
        oldY = cyList(i);

        %% Random displacement
        dxMove = moveAmplitude * (2*rand - 1);
        dyMove = moveAmplitude * (2*rand - 1);

        newX = mod(oldX + dxMove, Lx);
        newY = mod(oldY + dyMove, Ly);

        %% Check overlap with all other particles
        otherIdx = true(nParticles,1);
        otherIdx(i) = false;

        xOthers = cxList(otherIdx);
        yOthers = cyList(otherIdx);
        rOthers = rList(otherIdx);

        d = periodic_center_distance(newX, newY, xOthers, yOthers, Lx, Ly);
        requiredDistance = rList(i) + rOthers + minGap;

        if all(d >= requiredDistance)
            cxList(i) = newX;
            cyList(i) = newY;
        end
    end
end

%% ------------------------------------------------------------------------

function occCore = rebuild_periodic_occupancy(X, Y, cxList, cyList, rList, params)

    occCore = false(size(X));

    for i = 1:numel(rList)
        particleMask = periodic_particle_mask(X, Y, cxList(i), cyList(i), ...
            rList(i), params.Lx, params.Ly);
        occCore = occCore | particleMask;
    end
end


%% ------------------------------------------------------------------------

function [cxList, cyList, rList, occCore, currentArea, gcInfo, gcHistory] = ...
    growth_compression_PBC(cxList, cyList, rList, X, Y, params, domainArea)

    % Growth-compression stage for high volume fraction.
    % The algorithm gradually increases all particle radii and then uses
    % collective repulsive relaxation under PBC to remove the overlaps caused
    % by growth. Accepted growth steps increase the VF; failed steps are
    % backtracked by reducing the growth increment.

    if isfield(params,'activeTargetVF') && isfinite(params.activeTargetVF)
        targetVF = params.activeTargetVF;
    else
        targetVF = params.growthCompressionTargetVF;
    end
    maxSteps = params.growthCompressionMaxSteps;
    stepScale = params.growthCompressionStepScale;
    minStepScale = params.growthCompressionMinStepScale;
    backtrackFactor = params.growthCompressionBacktrackFactor;
    maxBacktracks = params.growthCompressionMaxBacktracks;

    gcInfo.cyclesUsed = 0;
    gcInfo.acceptedSteps = 0;
    gcInfo.reachedTarget = false;
    gcInfo.stopReason = "Not started";
    gcHistory.areaFraction = [];

    if isfield(params, 'growthCompressionNoProgressLimit')
        noProgressLimit = params.growthCompressionNoProgressLimit;
    else
        noProgressLimit = 30;
    end

    if isfield(params, 'growthCompressionMinVFImprovement')
        minVFImprovement = params.growthCompressionMinVFImprovement;
    else
        minVFImprovement = 1.0e-5;
    end

    if isfield(params, 'growthCompressionConvergenceWindow')
        convergenceWindow = params.growthCompressionConvergenceWindow;
    else
        convergenceWindow = 15;
    end

    if isfield(params, 'growthCompressionConvergenceTolerance')
        convergenceTolerance = params.growthCompressionConvergenceTolerance;
    else
        convergenceTolerance = 5.0e-4;
    end

    if isfield(params, 'growthCompressionMinVFForConvergenceCheck')
        minVFForConvergenceCheck = params.growthCompressionMinVFForConvergenceCheck;
    else
        minVFForConvergenceCheck = 0.70;
    end

    convergenceWindow = max(2, round(convergenceWindow));

    noProgressCounter = 0;
    bestVF = -inf;

    currentArea = sum(pi * rList.^2);
    currentVF = currentArea / domainArea;

    fprintf('Initial VF before growth-compression: %.5f\n', currentVF);
    fprintf('Requested VF after growth-compression: %.5f\n', targetVF);

    if currentVF >= targetVF
        gcInfo.reachedTarget = true;
        gcInfo.stopReason = "Specified VF already reached before growth-compression";
        occCore = rebuild_periodic_occupancy(X, Y, cxList, cyList, rList, params);
        return;
    end

    for stepID = 1:maxSteps

        gcInfo.cyclesUsed = stepID;

        currentArea = sum(pi * rList.^2);
        currentVF = currentArea / domainArea;

        if currentVF >= targetVF
            gcInfo.reachedTarget = true;
            gcInfo.stopReason = "Specified VF reached";
            break;
        end

        % Use a small growth increment but do not overshoot the target VF.
        requiredScaleToTarget = sqrt(targetVF / currentVF);
        trialScale = min(stepScale, requiredScaleToTarget);

        accepted = false;

        for backtrackID = 0:maxBacktracks

            if trialScale <= 1.0
                trialScale = 1.0 + eps;
            end

            if trialScale < minStepScale
                break;
            end

            rTrial = rList * trialScale;
            cxTrial = cxList;
            cyTrial = cyList;

            [cxTrial, cyTrial, success, finalMaxOverlap] = ...
                resolve_overlaps_repulsive_PBC(cxTrial, cyTrial, rTrial, params);

            if success
                cxList = cxTrial;
                cyList = cyTrial;
                rList  = rTrial;
                accepted = true;
                gcInfo.acceptedSteps = gcInfo.acceptedSteps + 1;
                break;
            else
                % If relaxation fails, reduce the growth increment and retry.
                trialScale = 1.0 + (trialScale - 1.0) * backtrackFactor;
            end
        end

        currentArea = sum(pi * rList.^2);
        currentVF = currentArea / domainArea;
        gcHistory.areaFraction(end+1,1) = currentVF;

        if accepted
            fprintf('GC step %4d: accepted, VF = %.5f, mean r = %.5f\n', ...
                stepID, currentVF, mean(rList));

            if currentVF > bestVF + minVFImprovement
                bestVF = currentVF;
                noProgressCounter = 0;
            else
                noProgressCounter = noProgressCounter + 1;
            end
        else
            gcInfo.stopReason = "No feasible growth increment after backtracking";
            fprintf('GC step %4d: stopped, no feasible growth increment. VF = %.5f\n', ...
                stepID, currentVF);
            break;
        end

        if currentVF >= targetVF
            gcInfo.reachedTarget = true;
            gcInfo.stopReason = "Specified VF reached";
            fprintf('Growth-compression target VF reached: %.5f\n', currentVF);
            break;
        end

        % Convergence/stagnation check based on the recent VF change.
        % This prevents the GC stage from continuing indefinitely when the
        % accepted growth steps become numerically small but physically
        % ineffective. If the target VF is not reached, the main routine will
        % proceed to the next stage, such as the hexagonal fallback.
        if numel(gcHistory.areaFraction) >= convergenceWindow && ...
                currentVF >= minVFForConvergenceCheck

            recentVFChange = gcHistory.areaFraction(end) - ...
                gcHistory.areaFraction(end - convergenceWindow + 1);

            if recentVFChange < convergenceTolerance
                gcInfo.stopReason = "GC converged/stagnated: recent VF change below threshold";
                fprintf('Growth-compression convergence detected.\n');
                fprintf('Recent VF change over last %d GC steps = %.6f < %.6f.\n', ...
                    convergenceWindow, recentVFChange, convergenceTolerance);
                fprintf('Switching to the next stage. Current VF = %.5f.\n', currentVF);
                break;
            end
        end

        if noProgressCounter >= noProgressLimit
            gcInfo.stopReason = "No meaningful VF improvement for the allowed number of GC iterations";
            fprintf('Growth-compression stopped after %d no-progress iterations. VF = %.5f\n', ...
                noProgressLimit, currentVF);
            break;
        end
    end

    currentArea = sum(pi * rList.^2);
    occCore = rebuild_periodic_occupancy(X, Y, cxList, cyList, rList, params);

    if ~gcInfo.reachedTarget
        if gcInfo.stopReason == "Not started"
            gcInfo.stopReason = "Maximum GC iterations reached before specified VF";
        end
        fprintf('Growth-compression ended at VF = %.5f before reaching %.5f.\n', ...
            currentArea/domainArea, targetVF);
        fprintf('Stop reason: %s\n', char(gcInfo.stopReason));
    end
end

%% ------------------------------------------------------------------------

function [cxList, cyList, success, maxOverlap] = ...
    resolve_overlaps_repulsive_PBC(cxList, cyList, rList, params)

    % Collective overlap resolution using pairwise repulsive displacements.
    % This is the core high-VF compaction step. It repeatedly detects all
    % overlapping pairs under periodic boundary conditions and pushes them
    % apart by a fraction of their overlap distance.

    Lx = params.Lx;
    Ly = params.Ly;
    minGap = params.minGap;
    maxIter = params.growthCompressionRelaxIterations;
    relaxFactor = params.growthCompressionRelaxFactor;
    tol = params.growthCompressionOverlapTolerance;

    nParticles = numel(rList);
    success = false;
    maxOverlap = inf;

    if nParticles < 2
        success = true;
        maxOverlap = 0.0;
        return;
    end

    % Pair index list is fixed because the particle number does not change.
    [I, J] = find(triu(true(nParticles), 1));

    for iter = 1:maxIter

        dx = cxList(J) - cxList(I);
        dy = cyList(J) - cyList(I);

        % Minimum-image convention for PBC.
        dx = dx - Lx * round(dx / Lx);
        dy = dy - Ly * round(dy / Ly);

        d = sqrt(dx.^2 + dy.^2);
        requiredDistance = rList(I) + rList(J) + minGap;
        overlap = requiredDistance - d;

        active = overlap > tol;

        if ~any(active)
            success = true;
            maxOverlap = 0.0;
            return;
        end

        maxOverlap = max(overlap(active));

        ia = I(active);
        ja = J(active);
        dxa = dx(active);
        dya = dy(active);
        da  = d(active);
        oa  = overlap(active);

        % Avoid division by zero for coincident centers.
        zeroD = da < eps;
        if any(zeroD)
            randomAngle = 2*pi*rand(sum(zeroD),1);
            dxa(zeroD) = cos(randomAngle);
            dya(zeroD) = sin(randomAngle);
            da(zeroD) = 1.0;
        end

        ux = dxa ./ da;
        uy = dya ./ da;

        % Split the required correction equally between the two particles.
        shift = 0.5 * relaxFactor * oa;

        moveX = zeros(nParticles,1);
        moveY = zeros(nParticles,1);

        moveX = moveX + accumarray(ia, -shift .* ux, [nParticles, 1], @sum, 0);
        moveY = moveY + accumarray(ia, -shift .* uy, [nParticles, 1], @sum, 0);
        moveX = moveX + accumarray(ja,  shift .* ux, [nParticles, 1], @sum, 0);
        moveY = moveY + accumarray(ja,  shift .* uy, [nParticles, 1], @sum, 0);

        % Limit very large moves to avoid oscillation.
        moveMagnitude = sqrt(moveX.^2 + moveY.^2);
        maxMove = 0.25 * mean(rList);
        tooLarge = moveMagnitude > maxMove;

        if any(tooLarge)
            moveX(tooLarge) = moveX(tooLarge) ./ moveMagnitude(tooLarge) * maxMove;
            moveY(tooLarge) = moveY(tooLarge) ./ moveMagnitude(tooLarge) * maxMove;
        end

        cxList = mod(cxList + moveX, Lx);
        cyList = mod(cyList + moveY, Ly);
    end

    % Final check after all iterations.
    dx = cxList(J) - cxList(I);
    dy = cyList(J) - cyList(I);
    dx = dx - Lx * round(dx / Lx);
    dy = dy - Ly * round(dy / Ly);
    d = sqrt(dx.^2 + dy.^2);
    requiredDistance = rList(I) + rList(J) + minGap;
    overlap = requiredDistance - d;
    maxOverlap = max([0; overlap]);
    success = maxOverlap <= tol;
end

%% ------------------------------------------------------------------------


function [cxList, cyList, rList, occCore, hexInfo] = ...
    generate_hexagonal_seeded_dense_PBC(params, X, Y, domainArea)

    % Periodic hexagonal-seeded dense generator with mandatory overlap audit.
    %
    % IMPORTANT:
    %   No theoretical packing-limit criterion is used here. The code accepts
    %   or rejects each candidate only through explicit periodic overlap
    %   checking. Therefore, if the requested VF is geometrically impossible
    %   for the selected radius/domain/gap, the algorithm will automatically
    %   return the highest overlap-free configuration it can find below the
    %   target instead of reporting an invalid overlapping packing.

    Lx = params.Lx;
    Ly = params.Ly;
    minGap = params.minGap;
    targetVF = params.hexagonalFallbackTargetVF;

    if params.particleMode == "monodisperse"
        rBase = params.rMono;
    else
        rBase = params.rMean;
    end

    requiredDistance = 2*rBase + minGap;
    nTarget = ceil(targetVF * domainArea / (pi*rBase^2));

    % Search range for integer periodic grids. This is intentionally broad,
    % but fast geometric screening is used before the expensive all-pair audit.
    maxCells = ceil(sqrt(max(nTarget,1))) * 5 + 40;

    bestBelow = [];
    bestAbove = [];
    bestBelowVF = -inf;
    bestAboveVF = inf;

    fprintf('Hexagonal fallback search started. Requested VF = %.5f, requested N = %d\n', ...
        targetVF, nTarget);
    fprintf('Every candidate grid will be accepted only if the full PBC overlap audit passes.\n');

    for nRow = 2:maxCells
        for nCol = 2:maxCells

            nGrid = nRow * nCol;

            sx = Lx / nCol;
            sy = Ly / nRow;

            % Fast geometric pre-screening only to reduce cost.
            % This does NOT replace the final overlap audit.
            if sx < requiredDistance
                continue;
            end

            dDiag = sqrt((0.5*sx)^2 + sy^2);
            if dDiag < requiredDistance
                continue;
            end

            % Generate the full periodic staggered grid.
            [cxAll, cyAll] = generate_periodic_hex_coordinates( ...
                nRow, nCol, sx, sy, Lx, Ly);

            rAll = rBase * ones(numel(cxAll),1);

            % Mandatory full-grid overlap audit under PBC. This catches:
            %   1. particle-particle overlap inside the box,
            %   2. particle-mirrored-particle overlap across opposite edges,
            %   3. incompatible top-bottom row staggering,
            %   4. incompatible left-right periodic contacts.
            [hasOverlapFull, maxOverlapFull, minClearanceFull, nOverlapPairsFull] = ...
                final_overlap_audit_PBC(cxAll, cyAll, rAll, params);

            if hasOverlapFull
                continue;
            end

            % If this full grid contains more particles than needed, retain a
            % subset of exactly nTarget particles. Any subset of an overlap-free
            % full grid remains overlap-free, but an audit is still performed.
            if nGrid >= nTarget
                nKeep = nTarget;
            else
                nKeep = nGrid;
            end

            if nKeep < nGrid
                % Deterministic, spatially distributed subset selection.
                % The sorted index list avoids random run-to-run changes.
                keep = round(linspace(1, nGrid, nKeep));
                keep = unique(keep, 'stable');

                % If uniqueness reduced the count, fill missing indices.
                if numel(keep) < nKeep
                    missing = setdiff(1:nGrid, keep, 'stable');
                    keep = [keep(:); missing(1:(nKeep-numel(keep)))'];
                    keep = sort(keep);
                end
            else
                keep = 1:nGrid;
            end

            cxTrial = cxAll(keep);
            cyTrial = cyAll(keep);
            rTrial = rBase * ones(numel(cxTrial),1);

            [hasOverlap, maxOverlap, minClearance, nOverlapPairs] = ...
                final_overlap_audit_PBC(cxTrial, cyTrial, rTrial, params);

            if hasOverlap
                % This should rarely happen if the full grid passed, but this
                % second audit makes the final decision fully explicit.
                continue;
            end

            VFtrial = sum(pi*rTrial.^2) / domainArea;

            candidate.nRow = nRow;
            candidate.nCol = nCol;
            candidate.sx = sx;
            candidate.sy = sy;
            candidate.cx = cxTrial;
            candidate.cy = cyTrial;
            candidate.r  = rTrial;
            candidate.VF = VFtrial;
            candidate.nParticles = numel(rTrial);
            candidate.minClearance = minClearance;
            candidate.maxOverlap = maxOverlap;
            candidate.nOverlapPairs = nOverlapPairs;
            candidate.fullGridMinClearance = minClearanceFull;
            candidate.fullGridMaxOverlap = maxOverlapFull;
            candidate.fullGridOverlapPairs = nOverlapPairsFull;

            if VFtrial >= targetVF
                % Choose the closest valid packing at or above target.
                if VFtrial < bestAboveVF
                    bestAboveVF = VFtrial;
                    bestAbove = candidate;
                end
            else
                % Otherwise keep the densest valid packing below target.
                if VFtrial > bestBelowVF
                    bestBelowVF = VFtrial;
                    bestBelow = candidate;
                end
            end
        end
    end

    if ~isempty(bestAbove)
        best = bestAbove;
        reachedTarget = true;
        stopReason = "Specified VF reached by overlap-audited periodic hexagonal fallback";
    elseif ~isempty(bestBelow)
        best = bestBelow;
        reachedTarget = false;
        stopReason = "Requested VF could not be reached without overlap; densest audited fallback below target returned";
    else
        error('Hexagonal fallback failed: no overlap-free periodic grid was found for the current radius, gap, and domain.');
    end

    cxList = best.cx;
    cyList = best.cy;
    rList  = best.r;

    % Optional small jitter. After jitter, the configuration is relaxed and
    % audited again. If the audit fails, the non-jittered audited grid is kept.
    jitterFraction = params.hexagonalFallbackJitterFraction;

    if jitterFraction > 0
        cxBeforeJitter = cxList;
        cyBeforeJitter = cyList;

        jitterAmp = jitterFraction * requiredDistance;
        cxList = mod(cxList + jitterAmp*(2*rand(size(cxList))-1), Lx);
        cyList = mod(cyList + jitterAmp*(2*rand(size(cyList))-1), Ly);

        paramsLocal = params;
        paramsLocal.growthCompressionRelaxIterations = params.hexagonalFallbackRelaxIterations;
        paramsLocal.growthCompressionRelaxFactor = params.hexagonalFallbackRelaxFactor;
        paramsLocal.growthCompressionOverlapTolerance = 1.0e-8;

        [cxList, cyList, success, maxOverlapAfterRelax] = ...
            resolve_overlaps_repulsive_PBC(cxList, cyList, rList, paramsLocal);

        [hasOverlapAfterJitter, ~, ~, ~] = ...
            final_overlap_audit_PBC(cxList, cyList, rList, params);

        if ~success || hasOverlapAfterJitter
            warning('Hexagonal fallback jitter caused overlap after relaxation. Reverting to audited non-jittered lattice. Max overlap = %.3e.', maxOverlapAfterRelax);
            cxList = cxBeforeJitter;
            cyList = cyBeforeJitter;
        end
    end

    % Mandatory final audit of the actual returned configuration.
    [hasOverlapFinal, maxOverlapFinal, minClearanceFinal, nOverlapPairsFinal, worstPairFinal] = ...
        final_overlap_audit_PBC(cxList, cyList, rList, params);

    if hasOverlapFinal
        error(['Hexagonal fallback rejected by final audit. ', ...
               'Overlapping pairs = %d, max overlap = %.6e, worst pair = [%d, %d].'], ...
               nOverlapPairsFinal, maxOverlapFinal, worstPairFinal(1), worstPairFinal(2));
    end

    occCore = rebuild_periodic_occupancy(X, Y, cxList, cyList, rList, params);

    finalVF = sum(pi*rList.^2)/domainArea;

    hexInfo.finalVF = finalVF;
    hexInfo.nParticles = numel(rList);
    hexInfo.nRow = best.nRow;
    hexInfo.nCol = best.nCol;
    hexInfo.hasOverlap = hasOverlapFinal;
    hexInfo.maxOverlap = maxOverlapFinal;
    hexInfo.minClearance = minClearanceFinal;
    hexInfo.nOverlappingPairs = nOverlapPairsFinal;
    hexInfo.worstPair = worstPairFinal;
    hexInfo.reachedTarget = reachedTarget;
    hexInfo.stopReason = stopReason;

    fprintf('Final PBC overlap audit passed.\n');
    fprintf('Minimum clearance = %.6e, overlapping pairs = %d\n', ...
        minClearanceFinal, nOverlapPairsFinal);
    fprintf('Hexagonal fallback grid: nRow = %d, nCol = %d, N = %d, VF = %.5f\n', ...
        best.nRow, best.nCol, numel(rList), finalVF);

    if ~reachedTarget
        fprintf('Requested VF = %.5f was not accepted because every higher candidate failed the overlap audit.\n', targetVF);
        fprintf('Returned highest audited non-overlapping fallback VF = %.5f.\n', finalVF);
    end
end

%% ------------------------------------------------------------------------

function [cx, cy] = generate_periodic_hex_coordinates(nRow, nCol, sx, sy, Lx, Ly)

    % Generates staggered hexagonal coordinates in a periodic rectangular box.

    cx = zeros(nRow*nCol,1);
    cy = zeros(nRow*nCol,1);

    counter = 0;

    for row = 1:nRow
        y = (row - 0.5) * sy;
        xShift = 0.5 * sx * mod(row-1, 2);

        for col = 1:nCol
            x = mod((col - 0.5) * sx + xShift, Lx);
            counter = counter + 1;
            cx(counter) = x;
            cy(counter) = y;
        end
    end
end

%% ------------------------------------------------------------------------

%% ------------------------------------------------------------------------

function distCenter = periodic_distance_transform(occCore, pixelSize)

    [Ny, Nx] = size(occCore);

    occTiled = repmat(occCore, 3, 3);
    distTiled = bwdist(occTiled) * pixelSize;

    rowStart = Ny + 1;
    rowEnd   = 2 * Ny;

    colStart = Nx + 1;
    colEnd   = 2 * Nx;

    distCenter = distTiled(rowStart:rowEnd, colStart:colEnd);
end

%% ------------------------------------------------------------------------

function mask = periodic_particle_mask(X, Y, cx, cy, r, Lx, Ly)

    dx = abs(X - cx);
    dy = abs(Y - cy);

    dx = min(dx, Lx - dx);
    dy = min(dy, Ly - dy);

    mask = dx.^2 + dy.^2 <= r^2;
end

%% ------------------------------------------------------------------------

function d = periodic_center_distance(xNew, yNew, xList, yList, Lx, Ly)

    dx = abs(xList - xNew);
    dx = min(dx, Lx - dx);

    dy = abs(yList - yNew);
    dy = min(dy, Ly - dy);

    d = sqrt(dx.^2 + dy.^2);
end

%% ------------------------------------------------------------------------

function selectedIndices = select_candidate_indices(admissibleIndices, distOcc, ...
    nCandidates, useWeightedSampling)

    nAvailable = numel(admissibleIndices);
    nCandidates = min(nCandidates, nAvailable);

    if nAvailable == 0 || nCandidates == 0
        selectedIndices = [];
        return;
    end

    %% Uniform random selection
    if ~useWeightedSampling
        selectedLocal = randperm(nAvailable, nCandidates);
        selectedIndices = admissibleIndices(selectedLocal);
        return;
    end

    %% Weighted random selection
    weights = distOcc(admissibleIndices);
    weights = weights(:);

    if isempty(weights) || all(isinf(weights)) || all(isnan(weights))
        selectedLocal = randperm(nAvailable, nCandidates);
        selectedIndices = admissibleIndices(selectedLocal);
        return;
    end

    finiteWeights = weights(isfinite(weights));

    if isempty(finiteWeights)
        selectedLocal = randperm(nAvailable, nCandidates);
        selectedIndices = admissibleIndices(selectedLocal);
        return;
    end

    maxFiniteWeight = max(finiteWeights);

    weights(isinf(weights)) = maxFiniteWeight;
    weights(isnan(weights)) = 0;
    weights(weights < 0) = 0;

    if sum(weights) <= 0
        selectedLocal = randperm(nAvailable, nCandidates);
        selectedIndices = admissibleIndices(selectedLocal);
        return;
    end

    selectedLocal = weighted_sample_without_replacement(weights, nCandidates);
    selectedIndices = admissibleIndices(selectedLocal);
end

%% ------------------------------------------------------------------------

function selected = weighted_sample_without_replacement(weights, n)

    weights = weights(:);

    nTotal = numel(weights);
    n = min(n, nTotal);

    weights(isnan(weights)) = 0;

    finiteWeights = weights(isfinite(weights));

    if isempty(finiteWeights)
        selected = randperm(nTotal, n)';
        return;
    end

    maxFiniteWeight = max(finiteWeights);

    weights(isinf(weights)) = maxFiniteWeight;
    weights(weights < 0) = 0;

    if sum(weights) <= 0
        selected = randperm(nTotal, n)';
        return;
    end

    weights = weights ./ sum(weights);

    % Efraimidis-Spirakis weighted sampling without replacement.
    u = rand(nTotal, 1);

    keys = zeros(nTotal, 1);
    positive = weights > 0;

    keys(positive) = u(positive).^(1 ./ weights(positive));
    keys(~positive) = -inf;

    [~, order] = sort(keys, 'descend');

    selected = order(1:n);
end


%% ------------------------------------------------------------------------

function selectedIndices = select_candidate_indices_LV_hybrid( ...
    admissibleIndices, distOcc, nCandidates, params, X, Y)

    % Hybrid largest-void + weighted-random candidate selection.
    %
    % A prescribed fraction of candidate centers is selected from local
    % maxima of the periodic distance-transform field. These points
    % represent the largest remaining admissible voids. The remaining
    % candidates are selected using the weighted random selector to preserve
    % stochastic morphology.

    nAvailable = numel(admissibleIndices);
    nCandidates = min(nCandidates, nAvailable);

    if nAvailable == 0 || nCandidates == 0
        selectedIndices = [];
        return;
    end

    %% If largest-void insertion is disabled, use the original selector

    if ~isfield(params, 'useLargestVoidInsertion') || ~params.useLargestVoidInsertion

        selectedIndices = select_candidate_indices( ...
            admissibleIndices, distOcc, nCandidates, params.useWeightedSampling);

        return;
    end

    %% If distance field is not informative, use original selector

    weights = distOcc(admissibleIndices);

    if isempty(weights) || all(isinf(weights)) || all(isnan(weights))

        selectedIndices = select_candidate_indices( ...
            admissibleIndices, distOcc, nCandidates, params.useWeightedSampling);

        return;
    end

    %% Number of largest-void candidates

    if isfield(params, 'largestVoidFraction')
        largestVoidFraction = params.largestVoidFraction;
    else
        largestVoidFraction = 0.70;
    end

    largestVoidFraction = max(0, min(1, largestVoidFraction));

    nVoid = round(largestVoidFraction * nCandidates);
    nVoid = min(max(nVoid, 1), nCandidates);

    %% Largest-void candidate centers

    voidIndices = find_largest_void_candidates( ...
        admissibleIndices, distOcc, nVoid, params, X, Y);

    %% Weighted-random candidate centers for the remaining quota

    nRemaining = nCandidates - numel(voidIndices);

    if nRemaining > 0

        remainingPool = setdiff(admissibleIndices, voidIndices, 'stable');

        if ~isempty(remainingPool)
            randomIndices = select_candidate_indices( ...
                remainingPool, distOcc, nRemaining, params.useWeightedSampling);
        else
            randomIndices = [];
        end
    else
        randomIndices = [];
    end

    %% Combine candidates

    selectedIndices = unique([voidIndices(:); randomIndices(:)], 'stable');

    %% Fill if unique operation reduced candidate count

    if numel(selectedIndices) < nCandidates

        missingCount = nCandidates - numel(selectedIndices);
        remainingPool = setdiff(admissibleIndices, selectedIndices, 'stable');

        if ~isempty(remainingPool)
            missingCount = min(missingCount, numel(remainingPool));
            extraLocal = randperm(numel(remainingPool), missingCount);
            extraIndices = remainingPool(extraLocal);
            selectedIndices = [selectedIndices(:); extraIndices(:)];
        end
    end
end

%% ------------------------------------------------------------------------

function voidIndices = find_largest_void_candidates( ...
    admissibleIndices, distOcc, nVoid, params, X, Y)

    % Detects and ranks the largest void centers using local maxima of
    % the periodic Euclidean distance-transform field.

    if nVoid <= 0
        voidIndices = [];
        return;
    end

    %% Build admissible mask

    admissibleMask = false(size(distOcc));
    admissibleMask(admissibleIndices) = true;

    %% Clean distance field

    distClean = distOcc;

    finiteValues = distClean(isfinite(distClean));

    if isempty(finiteValues)
        voidIndices = [];
        return;
    end

    maxFinite = max(finiteValues);

    distClean(isinf(distClean)) = maxFinite;
    distClean(isnan(distClean)) = 0;

    %% Detect local maxima

    if exist('imregionalmax', 'file') == 2
        localMaxima = imregionalmax(distClean);
    else
        % Fallback if imregionalmax is unavailable.
        % The sorted admissible distance list still prioritizes largest voids.
        localMaxima = false(size(distClean));
    end

    voidMask = localMaxima & admissibleMask;
    candidateVoidIndices = find(voidMask);

    %% If local maxima are not found, use largest admissible distance values

    if isempty(candidateVoidIndices)
        [~, order] = sort(distClean(admissibleIndices), 'descend');
        sortedIndices = admissibleIndices(order);
    else
        voidScores = distClean(candidateVoidIndices);
        [~, order] = sort(voidScores, 'descend');
        sortedIndices = candidateVoidIndices(order);
    end

    %% Spatially thin largest void candidates to avoid clustering

    if isfield(params, 'largestVoidMinSeparation') && params.largestVoidMinSeparation > 0

        voidIndices = spatially_thin_void_candidates( ...
            sortedIndices, nVoid, params.largestVoidMinSeparation, ...
            X, Y, params.Lx, params.Ly);
    else
        nSelect = min(nVoid, numel(sortedIndices));
        voidIndices = sortedIndices(1:nSelect);
    end
end

%% ------------------------------------------------------------------------

function selectedVoidIndices = spatially_thin_void_candidates( ...
    sortedIndices, nVoid, minSeparation, X, Y, Lx, Ly)

    % Keeps largest void centers while enforcing a minimum separation
    % between selected candidate centers under periodic boundary conditions.

    selectedVoidIndices = [];
    selectedX = [];
    selectedY = [];

    for k = 1:numel(sortedIndices)

        idx = sortedIndices(k);

        xNew = X(idx);
        yNew = Y(idx);

        if isempty(selectedX)

            selectedVoidIndices(end+1,1) = idx;
            selectedX(end+1,1) = xNew;
            selectedY(end+1,1) = yNew;

        else

            d = periodic_center_distance(xNew, yNew, selectedX, selectedY, Lx, Ly);

            if all(d >= minSeparation)

                selectedVoidIndices(end+1,1) = idx;
                selectedX(end+1,1) = xNew;
                selectedY(end+1,1) = yNew;
            end
        end

        if numel(selectedVoidIndices) >= nVoid
            break;
        end
    end

    %% If thinning was too strict, fill from the ranked list

    if numel(selectedVoidIndices) < nVoid

        remaining = setdiff(sortedIndices, selectedVoidIndices, 'stable');
        nMissing = nVoid - numel(selectedVoidIndices);
        nMissing = min(nMissing, numel(remaining));

        selectedVoidIndices = [selectedVoidIndices(:); remaining(1:nMissing)];
    end
end

%% ------------------------------------------------------------------------

function r = generate_particle_radius(params)

    if params.particleMode == "monodisperse"
        r = params.rMono;

    elseif params.particleMode == "polydisperse"
        accepted = false;

        while ~accepted
            rTry = params.rMean + params.rStd * randn;

            if rTry >= params.rMin && rTry <= params.rMax
                r = rTry;
                accepted = true;
            end
        end

    else
        error('particleMode must be either "monodisperse" or "polydisperse".');
    end
end

%% ------------------------------------------------------------------------

function plot_LV_R_DTPAA_PBC_result(result, params)

    particleTable = result.particleTable;

    %% Plot 1: occupancy map
    figure('Color','w');
    imagesc(result.x, result.y, result.occCore);
    axis equal tight;
    set(gca, 'YDir', 'normal');
    xlabel('x');
    ylabel('y');
    title('Periodic hard-particle occupancy map generated by VGPP');
    colorbar;

    %% Plot 2: particle distribution with periodic images
    figure('Color','w');
    hold on;
    box on;

    theta = linspace(0, 2*pi, 120);

    for i = 1:height(particleTable)
        cx = particleTable.x(i);
        cy = particleTable.y(i);
        r  = particleTable.radius(i);
        draw_periodic_particle(cx, cy, r, params.Lx, params.Ly, theta);
    end

    axis equal;
    xlim([0 params.Lx]);
    ylim([0 params.Ly]);
    xlabel('x');
    ylabel('y');
    title(sprintf('VGPP distribution, maximum area fraction = %.4f', ...
        result.maxAreaFraction));

    %% Plot 3: particle radius distribution
    figure('Color','w');
    histogram(particleTable.radius, 20);
    xlabel('Particle radius');
    ylabel('Frequency');
    title('Particle-size distribution');

    %% Plot 4: area fraction evolution by cycle
    if ~isempty(result.historyTable)
        figure('Color','w');
        plot(result.historyTable.Cycle, result.historyTable.AreaFraction, ...
            'LineWidth', 2);
        xlabel('Insertion / relaxation cycle');
        ylabel('Area/volume fraction');
        title('Area/volume fraction evolution during VGPP generation');
        grid off;
        box on;
    end

    %% Plot 5: nearest-neighbor distance under PBC
    if height(particleTable) > 1
        nParticles = height(particleTable);
        nearestDistance = zeros(nParticles,1);

        for i = 1:nParticles
            xi = particleTable.x(i);
            yi = particleTable.y(i);
            xOthers = particleTable.x;
            yOthers = particleTable.y;

            d = periodic_center_distance(xi, yi, xOthers, yOthers, ...
                params.Lx, params.Ly);
            d(i) = inf;
            nearestDistance(i) = min(d);
        end

        figure('Color','w');
        histogram(nearestDistance, 20);
        xlabel('Nearest-neighbor center distance under PBC');
        ylabel('Frequency');
        title('Periodic nearest-neighbor distance distribution');
    end
end

%% ------------------------------------------------------------------------

%% ------------------------------------------------------------------------


function result = run_RSA_PBC(params)

    % Classical Random Sequential Adsorption with periodic boundary condition.
    % Computation time is not recorded or used.

    Lx = params.Lx;
    Ly = params.Ly;
    domainArea = Lx * Ly;

    cxList = [];
    cyList = [];
    rList  = [];

    currentArea = 0.0;
    trialCounter = 0;
    failedTrials = 0;

    trialHistory = [];
    particleHistory = [];
    areaFractionHistory = [];

    while trialCounter < params.maxTrials && failedTrials < params.maxFailedTrials

        trialCounter = trialCounter + 1;

        rNew = generate_particle_radius(params);
        xNew = rand * Lx;
        yNew = rand * Ly;

        accept = true;

        if ~isempty(cxList)
            d = periodic_center_distance(xNew, yNew, cxList, cyList, Lx, Ly);
            requiredDistance = rList + rNew + params.minGap;

            if any(d < requiredDistance)
                accept = false;
            end
        end

        if accept
            cxList(end+1,1) = xNew;
            cyList(end+1,1) = yNew;
            rList(end+1,1)  = rNew;

            currentArea = currentArea + pi * rNew^2;
            failedTrials = 0;

            trialHistory(end+1,1) = trialCounter;
            particleHistory(end+1,1) = numel(rList);
            areaFractionHistory(end+1,1) = currentArea / domainArea;

            % Stop at the exact prescribed particle count when supplied.
            if isfield(params,'targetN') && isfinite(params.targetN) && ...
                    numel(rList) >= params.targetN
                break;
            elseif isfield(params,'targetVF') && isfinite(params.targetVF) && ...
                    currentArea/domainArea >= params.targetVF - 1.0e-12
                break;
            end
        else
            failedTrials = failedTrials + 1;
        end
    end

    maxAreaFraction = currentArea / domainArea;
    particleID = (1:numel(cxList))';

    particleTable = table(particleID, cxList, cyList, rList, ...
        'VariableNames', {'ID', 'x', 'y', 'radius'});

    historyTable = table(trialHistory, particleHistory, areaFractionHistory, ...
        'VariableNames', {'Trial', 'NumberOfParticles', 'AreaFraction'});

    result.particleTable = particleTable;
    result.historyTable = historyTable;
    result.maxAreaFraction = maxAreaFraction;
    result.trials = trialCounter;
    result.failedTrials = failedTrials;
end

%% ------------------------------------------------------------------------

function plot_method_comparison_particles(result_DTPAA, result_RSA, params)

    theta = linspace(0, 2*pi, 120);

    figure('Color','w');

    %% VGPP plot
    subplot(1,2,1);
    hold on;
    box on;

    particleTable = result_DTPAA.particleTable;

    for i = 1:height(particleTable)
        cx = particleTable.x(i);
        cy = particleTable.y(i);
        r  = particleTable.radius(i);
        draw_periodic_particle(cx, cy, r, params.Lx, params.Ly, theta);
    end

    axis equal;
    xlim([0 params.Lx]);
    ylim([0 params.Ly]);
    xlabel('x');
    ylabel('y');
    title(sprintf('Best VGPP, VF = %.4f', ...
        result_DTPAA.maxAreaFraction));

    %% RSA-PBC plot
    subplot(1,2,2);
    hold on;
    box on;

    particleTable = result_RSA.particleTable;

    for i = 1:height(particleTable)
        cx = particleTable.x(i);
        cy = particleTable.y(i);
        r  = particleTable.radius(i);
        draw_periodic_particle(cx, cy, r, params.Lx, params.Ly, theta);
    end

    axis equal;
    xlim([0 params.Lx]);
    ylim([0 params.Ly]);
    xlabel('x');
    ylabel('y');
    title(sprintf('Classical RSA-PBC, VF = %.4f', ...
        result_RSA.maxAreaFraction));
end

%% ------------------------------------------------------------------------

function draw_periodic_particle(cx, cy, r, Lx, Ly, theta)

    shiftX = [-Lx, 0, Lx];
    shiftY = [-Ly, 0, Ly];

    for sx = shiftX
        for sy = shiftY

            cxImage = cx + sx;
            cyImage = cy + sy;

            if cxImage + r < 0
                continue;
            end

            if cxImage - r > Lx
                continue;
            end

            if cyImage + r < 0
                continue;
            end

            if cyImage - r > Ly
                continue;
            end

            xp = cxImage + r * cos(theta);
            yp = cyImage + r * sin(theta);

            fill(xp, yp, [0.75 0.75 0.75], ...
                'EdgeColor', 'k', ...
                'LineWidth', 0.5);
        end
    end
end


%% ------------------------------------------------------------------------

function plot_largest_voids(result, params, nVoidToShow)

    % Visualizes the largest residual void centers in the final configuration.

    if nargin < 3
        nVoidToShow = 30;
    end

    X = result.X;
    Y = result.Y;
    occCore = result.occCore;

    distOcc = periodic_distance_transform(occCore, result.pixelSize);

    if params.particleMode == "monodisperse"
        rCheck = params.rMono;
    else
        rCheck = params.rMin;
    end

    admissible = distOcc >= (rCheck + params.minGap);
    admissibleIndices = find(admissible);

    if isempty(admissibleIndices)
        warning('No admissible void centers found for plotting.');
        return;
    end

    voidIndices = find_largest_void_candidates( ...
        admissibleIndices, distOcc, nVoidToShow, params, X, Y);

    figure('Color','w');
    hold on;
    box on;

    imagesc(result.x, result.y, distOcc);
    set(gca, 'YDir', 'normal');
    axis equal tight;
    colorbar;

    scatter(X(voidIndices), Y(voidIndices), ...
        70, 'r', 'filled', 'MarkerEdgeColor', 'k');

    xlabel('x');
    ylabel('y');
    title('Largest detected residual void centers');
end


%% ------------------------------------------------------------------------

function flag = has_admissible_voids(result, params)

    % Checks whether there is at least one residual void large enough
    % to place another particle.

    occCore = result.occCore;

    distOcc = periodic_distance_transform(occCore, result.pixelSize);

    if params.particleMode == "monodisperse"
        rCheck = params.rMono;
    else
        rCheck = params.rMin;
    end

    admissible = distOcc >= (rCheck + params.minGap);

    flag = any(admissible(:));
end



function plot_final_best_vs_RSA_all_useful(result_DTPAA, result_RSA, params, comparisonTable)

    % Generates a compact set of useful diagnostic plots comparing the final
    % best VGPP result with the classical RSA-PBC result.
    % Left panels always show the final best setting; right panels show RSA.

    theta = linspace(0, 2*pi, 120);
    finalTable = result_DTPAA.particleTable;
    rsaTable = result_RSA.particleTable;

    %% Audit both geometries before plotting
    [hasOverlapFinal, maxOverlapFinal, minClearanceFinal, nOverlapFinal, worstPairFinal] = ...
        final_overlap_audit_PBC(finalTable.x, finalTable.y, finalTable.radius, params);

    [hasOverlapRSA, maxOverlapRSA, minClearanceRSA, nOverlapRSA, worstPairRSA] = ...
        final_overlap_audit_PBC(rsaTable.x, rsaTable.y, rsaTable.radius, params);

    Geometry = ["Final best VGPP"; "Classical RSA-PBC"];
    VF = [result_DTPAA.maxAreaFraction; result_RSA.maxAreaFraction];
    NumberOfParticles = [height(finalTable); height(rsaTable)];
    MeanRadius = [mean(finalTable.radius); mean(rsaTable.radius)];
    MinClearance = [minClearanceFinal; minClearanceRSA];
    MaxOverlap = [maxOverlapFinal; maxOverlapRSA];
    NumberOfOverlappingPairs = [nOverlapFinal; nOverlapRSA];
    HasOverlap = [hasOverlapFinal; hasOverlapRSA];

    auditComparisonTable = table(Geometry, VF, NumberOfParticles, MeanRadius, ...
        HasOverlap, NumberOfOverlappingPairs, MaxOverlap, MinClearance);

    disp('Final best versus RSA overlap-audit comparison:');
    disp(auditComparisonTable);
    writetable(auditComparisonTable, 'Final_best_vs_RSA_overlap_audit_table.csv');

    %% Plot A: side-by-side particle maps, final best on left and RSA on right
    figure('Color','w', 'Name', 'Final best versus RSA particle maps');

    subplot(1,2,1);
    hold on;
    box on;
    for i = 1:height(finalTable)
        draw_periodic_particle(finalTable.x(i), finalTable.y(i), ...
            finalTable.radius(i), params.Lx, params.Ly, theta);
    end
    axis equal;
    xlim([0 params.Lx]);
    ylim([0 params.Ly]);
    xlabel('x');
    ylabel('y');
    title(sprintf('Final best setting, VF = %.4f, N = %d', ...
        result_DTPAA.maxAreaFraction, height(finalTable)));

    subplot(1,2,2);
    hold on;
    box on;
    for i = 1:height(rsaTable)
        draw_periodic_particle(rsaTable.x(i), rsaTable.y(i), ...
            rsaTable.radius(i), params.Lx, params.Ly, theta);
    end
    axis equal;
    xlim([0 params.Lx]);
    ylim([0 params.Ly]);
    xlabel('x');
    ylabel('y');
    title(sprintf('Classical RSA-PBC, VF = %.4f, N = %d', ...
        result_RSA.maxAreaFraction, height(rsaTable)));

    %% Plot B: occupancy maps side-by-side
    occFinal = rebuild_periodic_occupancy(result_DTPAA.X, result_DTPAA.Y, ...
        finalTable.x, finalTable.y, finalTable.radius, params);
    occRSA = rebuild_periodic_occupancy(result_DTPAA.X, result_DTPAA.Y, ...
        rsaTable.x, rsaTable.y, rsaTable.radius, params);

    figure('Color','w', 'Name', 'Final best versus RSA occupancy maps');

    subplot(1,2,1);
    imagesc(result_DTPAA.x, result_DTPAA.y, occFinal);
    axis equal tight;
    set(gca, 'YDir', 'normal');
    xlabel('x');
    ylabel('y');
    title('Final best occupancy map');
    colorbar;

    subplot(1,2,2);
    imagesc(result_DTPAA.x, result_DTPAA.y, occRSA);
    axis equal tight;
    set(gca, 'YDir', 'normal');
    xlabel('x');
    ylabel('y');
    title('RSA occupancy map');
    colorbar;

    %% Plot C: maximum VF and particle count comparison
    figure('Color','w', 'Name', 'VF and particle-count comparison');

    subplot(1,2,1);
    bar(categorical(Geometry), VF);
    ylabel('Area/volume fraction');
    title('Final VF comparison');
    grid off;
    box on;

    subplot(1,2,2);
    bar(categorical(Geometry), NumberOfParticles);
    ylabel('Number of particles');
    title('Particle-count comparison');
    grid off;
    box on;

    %% Plot D: evolution history, final best versus RSA
    figure('Color','w', 'Name', 'VF evolution comparison');
    hold on;
    box on;
    grid off;

    if ~isempty(result_DTPAA.historyTable)
        plot(result_DTPAA.historyTable.Cycle, result_DTPAA.historyTable.AreaFraction, ...
            'LineWidth', 2, 'DisplayName', 'Final best method');
    end

    if ~isempty(result_RSA.historyTable)
        plot(result_RSA.historyTable.Trial, result_RSA.historyTable.AreaFraction, ...
            'LineWidth', 2, 'DisplayName', 'RSA-PBC');
    end

    xlabel('Cycle / trial number');
    ylabel('Area/volume fraction');
    title('VF evolution: final best method versus RSA');
    legend('Location','best');

    %% Plot E: nearest-neighbor center-distance distributions
    nnFinal = nearest_neighbor_distance_PBC(finalTable.x, finalTable.y, params);
    nnRSA = nearest_neighbor_distance_PBC(rsaTable.x, rsaTable.y, params);

    figure('Color','w', 'Name', 'Nearest-neighbor comparison');
    hold on;
    box on;
    histogram(nnFinal, 25, 'DisplayName', 'Final best');
    histogram(nnRSA, 25, 'DisplayName', 'RSA-PBC');
    xlabel('Nearest-neighbor center distance under PBC');
    ylabel('Frequency');
    title('Nearest-neighbor distance comparison');
    legend('Location','best');

    %% Plot F: surface-to-surface clearance distributions
    clearanceFinal = nearest_surface_clearance_PBC(finalTable.x, finalTable.y, finalTable.radius, params);
    clearanceRSA = nearest_surface_clearance_PBC(rsaTable.x, rsaTable.y, rsaTable.radius, params);

    figure('Color','w', 'Name', 'Nearest-clearance comparison');
    hold on;
    box on;
    histogram(clearanceFinal, 25, 'DisplayName', 'Final best');
    histogram(clearanceRSA, 25, 'DisplayName', 'RSA-PBC');
    xlabel('Nearest surface-to-surface clearance under PBC');
    ylabel('Frequency');
    title('Nearest clearance comparison');
    legend('Location','best');

    %% Plot G: radius distributions
    figure('Color','w', 'Name', 'Particle-radius comparison');
    hold on;
    box on;
    histogram(finalTable.radius, 25, 'DisplayName', 'Final best');
    histogram(rsaTable.radius, 25, 'DisplayName', 'RSA-PBC');
    xlabel('Particle radius');
    ylabel('Frequency');
    title('Particle-radius distribution comparison');
    legend('Location','best');

    %% Plot H: overlap-audit metrics
    figure('Color','w', 'Name', 'Overlap audit comparison');

    subplot(1,3,1);
    bar(categorical(Geometry), NumberOfOverlappingPairs);
    ylabel('Number of overlapping pairs');
    title('Overlap-pair count');
    grid off;
    box on;

    subplot(1,3,2);
    bar(categorical(Geometry), MaxOverlap);
    ylabel('Maximum overlap');
    title('Maximum overlap');
    grid off;
    box on;

    subplot(1,3,3);
    bar(categorical(Geometry), MinClearance);
    ylabel('Minimum clearance');
    title('Minimum clearance');
    grid off;
    box on;


    %% Plot I: reconstructed packing evolution from 50% VF to final maximum VF
    plot_packing_evolution_from_final_best(result_DTPAA, params, 0.50, result_DTPAA.maxAreaFraction);

    %% Save compact numerical comparison table as well
    if nargin >= 4 && ~isempty(comparisonTable)
        writetable(comparisonTable, 'Final_best_vs_RSA_VF_particle_count_table.csv');
    end

    fprintf('Useful final-best-versus-RSA plots and audit tables generated.\n');
end

%% ------------------------------------------------------------------------


function plot_packing_evolution_from_final_best(result_DTPAA, params, startVF, endVF)

    % Visualizes how the packing level increases from a selected lower VF
    % to the final best VF. Because the final dense stage may be generated
    % directly rather than saved as many intermediate coordinate states, this
    % figure reconstructs intermediate packing levels by progressively
    % revealing a uniformly distributed subset of particles from the final
    % accepted non-overlapping configuration. The final panel is the actual
    % final best geometry.

    particleTable = result_DTPAA.particleTable;
    domainArea = params.Lx * params.Ly;

    if isempty(particleTable) || height(particleTable) < 2
        warning('Packing-evolution plot skipped because fewer than two particles are available.');
        return;
    end

    finalVF = result_DTPAA.maxAreaFraction;
    endVF = min(endVF, finalVF);

    if endVF <= startVF
        warning('Packing-evolution plot skipped because final VF is not greater than the requested start VF.');
        return;
    end

    % Choose visually useful VF levels. The last level is always the final VF.
    stepVF = 0.05;
    vfLevels = startVF:stepVF:endVF;

    if isempty(vfLevels) || abs(vfLevels(end) - endVF) > 1.0e-8
        vfLevels = [vfLevels, endVF];
    else
        vfLevels(end) = endVF;
    end

    vfLevels = unique(round(vfLevels, 6), 'stable');
    nLevels = numel(vfLevels);

    % Build a spatially balanced randomized reveal order.
    % The previous sorted x-y key may produce diagonal/triangular artificial
    % empty bands. The new order divides the domain into spatial bins,
    % randomly shuffles particles inside each bin, and then reveals particles
    % in round-robin form across all bins. Therefore, every intermediate VF
    % level remains random-looking and spatially homogeneous.
    revealOrder = spatially_balanced_random_order_PBC(particleTable, params);

    nCols = min(4, nLevels);
    nRows = ceil(nLevels / nCols);
    theta = linspace(0, 2*pi, 80);

    figure('Color','w', 'Name', 'Packing evolution from 50 percent to final VF');

    evolutionVF = zeros(nLevels,1);
    evolutionN = zeros(nLevels,1);

    for k = 1:nLevels

        targetVF = vfLevels(k);

        % Select enough particles to approximate the requested VF. For a
        % polydisperse case, accumulate particle areas until the selected
        % area is closest to the requested VF.
        targetArea = targetVF * domainArea;
        selected = [];
        selectedArea = 0.0;

        for q = 1:numel(revealOrder)
            idx = revealOrder(q);
            trialArea = selectedArea + pi * particleTable.radius(idx)^2;

            if trialArea <= targetArea || isempty(selected)
                selected(end+1,1) = idx; %#ok<AGROW>
                selectedArea = trialArea;
            else
                % Decide whether adding this particle or stopping gives a
                % closer approximation to the requested VF.
                errWithout = abs(targetArea - selectedArea);
                errWith = abs(targetArea - trialArea);
                if errWith < errWithout && targetVF < finalVF
                    selected(end+1,1) = idx; %#ok<AGROW>
                    selectedArea = trialArea;
                end
                break;
            end
        end

        % The last panel should show the actual final best setting.
        if k == nLevels
            selected = (1:height(particleTable))';
            selectedArea = sum(pi * particleTable.radius.^2);
        end

        evolutionVF(k) = selectedArea / domainArea;
        evolutionN(k) = numel(selected);

        subplot(nRows, nCols, k);
        hold on;
        box on;

        for iSel = 1:numel(selected)
            i = selected(iSel);
            draw_periodic_particle(particleTable.x(i), particleTable.y(i), ...
                particleTable.radius(i), params.Lx, params.Ly, theta);
        end

        axis equal;
        xlim([0 params.Lx]);
        ylim([0 params.Ly]);
        set(gca, 'XTick', [], 'YTick', []);

        if k == nLevels
            title(sprintf('Final VF = %.4f, N = %d', evolutionVF(k), evolutionN(k)));
        else
            title(sprintf('VF \x2248 %.4f, N = %d', evolutionVF(k), evolutionN(k)));
        end
    end

    sgtitle('Spatially balanced reconstructed packing evolution from 50% VF to final best VF');

    PackingLevel = (1:nLevels)';
    TargetVF = vfLevels(:);
    ReconstructedVF = evolutionVF;
    NumberOfParticles = evolutionN;

    evolutionTable = table(PackingLevel, TargetVF, ReconstructedVF, NumberOfParticles);
    writetable(evolutionTable, 'Packing_evolution_50_to_finalVF_table.csv');

    fprintf('Packing-evolution plot generated from %.2f VF to %.5f VF.\n', startVF, finalVF);
    fprintf('Packing-evolution data saved to Packing_evolution_50_to_finalVF_table.csv.\n');
end

%% ------------------------------------------------------------------------

function revealOrder = spatially_balanced_random_order_PBC(particleTable, params)

    % Produces a randomized but spatially homogeneous particle reveal order.
    % This function is used only for reconstructed packing-evolution
    % visualization. It does not modify the final accepted particle geometry.
    %
    % The previous version used a simple bin round-robin strategy. That was
    % homogeneous, but still looked somewhat artificial for a dense hexagonal
    % final lattice. This revised version uses a hybrid stochastic-blue-noise
    % ordering:
    %   1. the domain is divided into bins for spatial balance;
    %   2. candidate particles are sampled randomly from under-represented
    %      bins;
    %   3. among those candidates, a blue-noise score favors particles far
    %      from already revealed particles;
    %   4. a controlled random term prevents overly regular visual evolution.
    %
    % The result is a more natural, random, and homogeneous apparent packing
    % evolution from lower VF to the final dense state.

    N = height(particleTable);

    if N == 0
        revealOrder = [];
        return;
    end

    %% Visualization-only user controls with safe defaults
    if isfield(params, 'evolutionRandomSeed')
        randomSeed = params.evolutionRandomSeed;
    else
        randomSeed = 2026;
    end

    if isfield(params, 'evolutionRandomnessStrength')
        randomnessStrength = params.evolutionRandomnessStrength;
    else
        randomnessStrength = 0.75;
    end

    if isfield(params, 'evolutionBlueNoiseStrength')
        blueNoiseStrength = params.evolutionBlueNoiseStrength;
    else
        blueNoiseStrength = 0.65;
    end

    if isfield(params, 'evolutionCandidatePoolSize')
        candidatePoolSize = params.evolutionCandidatePoolSize;
    else
        candidatePoolSize = 90;
    end

    if isfield(params, 'evolutionGridBins')
        nBins = params.evolutionGridBins;
    else
        nBins = max(6, ceil(sqrt(N/8)));
    end

    nBins = max(3, round(nBins));
    randomnessStrength = max(0, min(1.5, randomnessStrength));
    blueNoiseStrength = max(0, min(2.0, blueNoiseStrength));
    candidatePoolSize = max(10, round(candidatePoolSize));

    %% Save and restore RNG state so plotting does not alter later operations
    rngState = rng;
    rng(randomSeed, 'twister');

    x = particleTable.x(:);
    y = particleTable.y(:);

    % Randomly shift the bin grid under PBC. This avoids visual alignment
    % between the reveal order and the hexagonal rows/columns.
    phaseX = rand * params.Lx;
    phaseY = rand * params.Ly;

    xShifted = mod(x + phaseX, params.Lx);
    yShifted = mod(y + phaseY, params.Ly);

    binX = floor(xShifted / params.Lx * nBins) + 1;
    binY = floor(yShifted / params.Ly * nBins) + 1;

    binX = min(max(binX,1), nBins);
    binY = min(max(binY,1), nBins);

    binID = (binY - 1) * nBins + binX;
    nTotalBins = nBins * nBins;

    binParticles = cell(nTotalBins, 1);
    binPopulation = zeros(nTotalBins, 1);

    for b = 1:nTotalBins
        ids = find(binID == b);
        if ~isempty(ids)
            ids = ids(randperm(numel(ids)));
        end
        binParticles{b} = ids(:);
        binPopulation(b) = numel(ids);
    end

    %% Hybrid stochastic-blue-noise reveal order
    revealOrder = zeros(N,1);
    selectedX = zeros(N,1);
    selectedY = zeros(N,1);
    selectedCountByBin = zeros(nTotalBins,1);
    isRemaining = true(N,1);

    % Start from a random spatial bin to avoid systematic patterns.
    nonEmptyBins = find(binPopulation > 0);
    firstBin = nonEmptyBins(randi(numel(nonEmptyBins)));
    firstCandidates = binParticles{firstBin};
    firstID = firstCandidates(randi(numel(firstCandidates)));

    revealOrder(1) = firstID;
    selectedX(1) = x(firstID);
    selectedY(1) = y(firstID);
    selectedCountByBin(firstBin) = selectedCountByBin(firstBin) + 1;
    isRemaining(firstID) = false;
    binParticles{firstBin}(binParticles{firstBin} == firstID) = [];

    counter = 1;

    while counter < N

        remainingBins = find(~cellfun(@isempty, binParticles));

        if isempty(remainingBins)
            break;
        end

        % Target selected count in each bin at the current reveal step. Bins
        % with larger deficit are preferred, preserving global homogeneity.
        expectedSelected = (counter / N) * binPopulation;
        deficit = expectedSelected - selectedCountByBin;
        deficit(~ismember(1:nTotalBins, remainingBins)) = -inf;

        % Rank bins by deficit, but include randomization to avoid a rigid
        % round-robin order.
        finiteDeficit = deficit(isfinite(deficit));
        if isempty(finiteDeficit)
            candidateBins = remainingBins;
        else
            noiseAmplitude = randomnessStrength * max(std(finiteDeficit), 1.0e-6);
            noisyDeficit = deficit + noiseAmplitude * randn(size(deficit));
            [~, binOrder] = sort(noisyDeficit, 'descend');
            candidateBins = binOrder(ismember(binOrder, remainingBins));
        end

        % Build a randomized candidate pool from under-represented bins.
        candidateIDs = [];
        maxBinsToProbe = min(numel(candidateBins), max(6, ceil(candidatePoolSize/3)));

        for kk = 1:maxBinsToProbe
            b = candidateBins(kk);
            ids = binParticles{b};

            if isempty(ids)
                continue;
            end

            nTake = min(numel(ids), max(1, ceil(candidatePoolSize / maxBinsToProbe)));
            ids = ids(randperm(numel(ids), nTake));
            candidateIDs = [candidateIDs; ids(:)]; %#ok<AGROW>

            if numel(candidateIDs) >= candidatePoolSize
                break;
            end
        end

        if isempty(candidateIDs)
            remainingIDs = find(isRemaining);
            candidateIDs = remainingIDs(randperm(numel(remainingIDs), ...
                min(candidatePoolSize, numel(remainingIDs))));
        end

        candidateIDs = unique(candidateIDs, 'stable');

        %% Score candidates
        % Blue-noise term: prefer candidates farther away from already
        % revealed particles, preventing clustered appearance at low VF.
        if counter == 0
            minDistScore = ones(numel(candidateIDs),1);
        else
            minDistScore = zeros(numel(candidateIDs),1);
            for c = 1:numel(candidateIDs)
                id = candidateIDs(c);
                d = periodic_center_distance(x(id), y(id), ...
                    selectedX(1:counter), selectedY(1:counter), ...
                    params.Lx, params.Ly);
                minDistScore(c) = min(d);
            end
        end

        if max(minDistScore) > 0
            minDistScore = minDistScore ./ max(minDistScore);
        end

        % Bin-deficit term: keeps the reveal spatially balanced.
        candidateBinID = binID(candidateIDs);
        binScore = deficit(candidateBinID);
        binScore(~isfinite(binScore)) = 0;
        if max(abs(binScore)) > 0
            binScore = binScore ./ max(abs(binScore));
        end

        % Random term: makes the evolution look less deterministic.
        randomScore = rand(numel(candidateIDs),1);

        score = blueNoiseStrength * minDistScore + ...
                0.75 * binScore + ...
                randomnessStrength * randomScore;

        [~, bestLocal] = max(score);
        chosenID = candidateIDs(bestLocal);
        chosenBin = binID(chosenID);

        counter = counter + 1;
        revealOrder(counter) = chosenID;
        selectedX(counter) = x(chosenID);
        selectedY(counter) = y(chosenID);
        selectedCountByBin(chosenBin) = selectedCountByBin(chosenBin) + 1;
        isRemaining(chosenID) = false;

        binParticles{chosenBin}(binParticles{chosenBin} == chosenID) = [];
    end

    revealOrder = revealOrder(1:counter);

    % If any numerical edge case left missing particles, append them randomly.
    if numel(revealOrder) < N
        missing = setdiff((1:N)', revealOrder, 'stable');
        missing = missing(randperm(numel(missing)));
        revealOrder = [revealOrder(:); missing(:)];
    end

    %% Restore RNG state
    rng(rngState);
end

%% ------------------------------------------------------------------------

function nn = nearest_neighbor_distance_PBC(cxList, cyList, params)

    nParticles = numel(cxList);
    nn = NaN(nParticles,1);

    if nParticles < 2
        return;
    end

    for i = 1:nParticles
        d = periodic_center_distance(cxList(i), cyList(i), cxList, cyList, ...
            params.Lx, params.Ly);
        d(i) = inf;
        nn(i) = min(d);
    end
end

%% ------------------------------------------------------------------------

function clearance = nearest_surface_clearance_PBC(cxList, cyList, rList, params)

    nParticles = numel(cxList);
    clearance = NaN(nParticles,1);

    if nParticles < 2
        return;
    end

    for i = 1:nParticles
        d = periodic_center_distance(cxList(i), cyList(i), cxList, cyList, ...
            params.Lx, params.Ly);
        requiredDistance = rList(i) + rList + params.minGap;
        c = d - requiredDistance;
        c(i) = inf;
        clearance(i) = min(c);
    end
end

%% ------------------------------------------------------------------------

function [hasOverlap, maxOverlap, minClearance, nOverlappingPairs, worstPair] = ...
    final_overlap_audit_PBC(cxList, cyList, rList, params)

    % Explicit all-pair overlap audit under periodic boundary conditions.
    %
    % This function is the final authority for geometric validity. It checks
    % every particle pair using the minimum-image convention, so it detects
    % both ordinary in-domain overlaps and mirrored overlaps across opposite
    % periodic edges.

    Lx = params.Lx;
    Ly = params.Ly;
    minGap = params.minGap;

    nParticles = numel(rList);

    hasOverlap = false;
    maxOverlap = 0;
    minClearance = inf;
    nOverlappingPairs = 0;
    worstPair = [NaN, NaN];

    if nParticles < 2
        minClearance = NaN;
        return;
    end

    for i = 1:nParticles-1

        xNew = cxList(i);
        yNew = cyList(i);

        xOthers = cxList(i+1:end);
        yOthers = cyList(i+1:end);
        rOthers = rList(i+1:end);

        d = periodic_center_distance(xNew, yNew, xOthers, yOthers, Lx, Ly);

        requiredDistance = rList(i) + rOthers + minGap;
        clearance = d - requiredDistance;
        overlap = requiredDistance - d;

        minClearance = min(minClearance, min(clearance));

        if isfield(params,'overlapAuditTolerance')
            auditTol = params.overlapAuditTolerance;
        else
            auditTol = 1.0e-8;
        end
        overlappingLocal = overlap > auditTol;
        nOverlappingPairs = nOverlappingPairs + sum(overlappingLocal);

        [localMaxOverlap, localIdx] = max([0; overlap(:)]);

        if localMaxOverlap > maxOverlap
            maxOverlap = localMaxOverlap;

            if localIdx > 1
                worstPair = [i, i + localIdx - 1];
            end
        end
    end

    hasOverlap = nOverlappingPairs > 0;

    if isfield(params,'overlapAuditTolerance')
        auditTol = params.overlapAuditTolerance;
    else
        auditTol = 1.0e-8;
    end
    if maxOverlap <= auditTol
        maxOverlap = 0;
    end

    if isinf(minClearance)
        minClearance = NaN;
    end
end


%% ------------------------------------------------------------------------
% Literature-inspired statistical validation block
% -------------------------------------------------------------------------
% This block adds post-processing plots commonly used in recent RVE-generation
% studies: nearest-neighbor distance/orientation, Ripley's K/L function,
% radial distribution function, local volume-fraction uniformity, and a
% compact statistical validation table. It does not modify the generated
% geometry; it only analyzes the final VGPP and RSA-PBC particle coordinates.

function plot_literature_inspired_VGPP_statistics(result_VGPP, result_RSA, params)

    finalTable = result_VGPP.particleTable;
    rsaTable   = result_RSA.particleTable;

    if height(finalTable) < 2 || height(rsaTable) < 2
        warning('Literature-inspired statistical validation skipped because one geometry has fewer than two particles.');
        return;
    end

    fprintf('\nStarting literature-inspired statistical validation plots...\n');

    statsVGPP = compute_spatial_statistics_PBC(finalTable, params, 'VGPP');
    statsRSA  = compute_spatial_statistics_PBC(rsaTable,   params, 'RSA-PBC');

    %% Plot 1: normalized nearest-neighbor distance PDFs
    figure('Color','w', 'Name', 'Literature plot - nearest-neighbor distance PDF');
    hold on; box on; grid off;
    histogram(statsVGPP.normalizedNND, 30, 'Normalization','pdf', ...
        'DisplayName','VGPP');
    histogram(statsRSA.normalizedNND, 30, 'Normalization','pdf', ...
        'DisplayName','RSA-PBC');
    xline(1.0, '--', 'Contact limit', 'LabelVerticalAlignment','bottom');
    xlabel('Normalized nearest-neighbor distance, d_{NN}/(r_i+r_j+g_{min})');
    ylabel('Probability density');
    title('Nearest-neighbor distance distribution');
    legend('Location','best');

    %% Plot 2: nearest-neighbor orientation CDF versus CSR diagonal
    figure('Color','w', 'Name', 'Literature plot - nearest-neighbor orientation CDF');
    hold on; box on; grid off;
    plot(statsVGPP.orientationAngleDeg, statsVGPP.orientationCDF, 'LineWidth', 2, ...
        'DisplayName','VGPP');
    plot(statsRSA.orientationAngleDeg, statsRSA.orientationCDF, 'LineWidth', 2, ...
        'DisplayName','RSA-PBC');
    plot([0 180], [0 1], 'k--', 'LineWidth', 1.5, 'DisplayName','CSR diagonal');
    xlabel('Nearest-neighbor orientation angle, degree');
    ylabel('Cumulative probability');
    title('Nearest-neighbor orientation randomness');
    legend('Location','best');

    %% Plot 3: Ripley L(r)-r function under PBC
    figure('Color','w', 'Name', 'Literature plot - Ripley L minus r');
    hold on; box on; grid off;
    plot(statsVGPP.rK, statsVGPP.LminusR, 'LineWidth', 2, 'DisplayName','VGPP');
    plot(statsRSA.rK,  statsRSA.LminusR,  'LineWidth', 2, 'DisplayName','RSA-PBC');
    yline(0, 'k--', 'CSR reference', 'LineWidth', 1.5);
    xlabel('Search radius, r');
    ylabel('L(r) - r');
    title('Second-order spatial interaction');
    legend('Location','best');

    %% Plot 4: radial distribution function g(r)
    figure('Color','w', 'Name', 'Literature plot - radial distribution function');
    hold on; box on; grid off;
    plot(statsVGPP.rRDF, statsVGPP.rdf, 'LineWidth', 2, 'DisplayName','VGPP');
    plot(statsRSA.rRDF,  statsRSA.rdf,  'LineWidth', 2, 'DisplayName','RSA-PBC');
    yline(1, 'k--', 'CSR far-field reference', 'LineWidth', 1.5);
    xlabel('Pair distance, r');
    ylabel('g(r)');
    title('Radial distribution function');
    legend('Location','best');

    %% Plot 5: local volume-fraction maps
    figure('Color','w', 'Name', 'Literature plot - local VF uniformity maps');

    subplot(1,2,1);
    imagesc(statsVGPP.localVFMap);
    axis equal tight;
    set(gca,'YDir','normal');
    colorbar;
    title(sprintf('VGPP local VF map, RSD = %.3f', statsVGPP.localVFRSD));
    xlabel('x bin'); ylabel('y bin');

    subplot(1,2,2);
    imagesc(statsRSA.localVFMap);
    axis equal tight;
    set(gca,'YDir','normal');
    colorbar;
    title(sprintf('RSA-PBC local VF map, RSD = %.3f', statsRSA.localVFRSD));
    xlabel('x bin'); ylabel('y bin');

    %% Plot 6: compact morphology-quality metrics
    metricLabels = {'NND CV', 'Orientation KS', '|L-r| mean', 'Local VF RSD', 'RDF first peak'};
    metricNames = categorical(metricLabels);
    metricNames = reordercats(metricNames, metricLabels);

    vgppMetrics = [ ...
        statsVGPP.nndCV, ...
        statsVGPP.orientationKS, ...
        statsVGPP.meanAbsLminusR, ...
        statsVGPP.localVFRSD, ...
        statsVGPP.rdfFirstPeak];

    rsaMetrics = [ ...
        statsRSA.nndCV, ...
        statsRSA.orientationKS, ...
        statsRSA.meanAbsLminusR, ...
        statsRSA.localVFRSD, ...
        statsRSA.rdfFirstPeak];

    figure('Color','w', 'Name', 'Literature plot - morphology quality metrics');
    bar(metricNames, [vgppMetrics(:), rsaMetrics(:)]);
    ylabel('Metric value');
    title('Spatial-statistical morphology descriptors');
    legend({'VGPP','RSA-PBC'}, 'Location','best');
    grid off; box on;

    %% Save statistical summary table
    Method = ["VGPP"; "RSA-PBC"];
    VolumeFraction = [statsVGPP.volumeFraction; statsRSA.volumeFraction];
    NumberOfParticles = [statsVGPP.N; statsRSA.N];
    MeanRadius = [statsVGPP.meanRadius; statsRSA.meanRadius];
    MeanNND = [statsVGPP.meanNND; statsRSA.meanNND];
    NND_CV = [statsVGPP.nndCV; statsRSA.nndCV];
    MinClearance = [statsVGPP.minClearance; statsRSA.minClearance];
    OrientationKS = [statsVGPP.orientationKS; statsRSA.orientationKS];
    MeanAbsLminusR = [statsVGPP.meanAbsLminusR; statsRSA.meanAbsLminusR];
    LocalVFRSD = [statsVGPP.localVFRSD; statsRSA.localVFRSD];
    RDFFirstPeak = [statsVGPP.rdfFirstPeak; statsRSA.rdfFirstPeak];

    literatureStatsTable = table(Method, VolumeFraction, NumberOfParticles, ...
        MeanRadius, MeanNND, NND_CV, MinClearance, OrientationKS, ...
        MeanAbsLminusR, LocalVFRSD, RDFFirstPeak);

    disp('Literature-inspired spatial statistical validation table:');
    disp(literatureStatsTable);
    writetable(literatureStatsTable, 'VGPP_literature_inspired_spatial_statistics.csv');

    fprintf('Literature-inspired statistical validation completed.\n');
end

%% ------------------------------------------------------------------------

function stats = compute_spatial_statistics_PBC(particleTable, params, methodName)

    x = particleTable.x(:);
    y = particleTable.y(:);
    r = particleTable.radius(:);

    Lx = params.Lx;
    Ly = params.Ly;
    A  = Lx * Ly;
    N  = numel(r);

    %% All pair distances and pair indices under PBC
    [pairD, pairI, pairJ] = all_pair_distances_PBC(x, y, Lx, Ly);

    %% Nearest-neighbor distance and orientation
    nearestD = zeros(N,1);
    nearestRequired = zeros(N,1);
    nearestAngle = zeros(N,1);

    for i = 1:N
        dAll = periodic_center_distance(x(i), y(i), x, y, Lx, Ly);
        dAll(i) = inf;
        [nearestD(i), jMin] = min(dAll);
        nearestRequired(i) = r(i) + r(jMin) + params.minGap;

        dx = x(jMin) - x(i);
        dy = y(jMin) - y(i);
        dx = dx - Lx * round(dx / Lx);
        dy = dy - Ly * round(dy / Ly);

        angle = atan2d(dy, dx);
        angle = mod(angle, 180);   % undirected nearest-neighbor orientation
        nearestAngle(i) = angle;
    end

    normalizedNND = nearestD ./ nearestRequired;

    angleGrid = linspace(0, 180, 181)';
    orientationCDF = zeros(size(angleGrid));
    for k = 1:numel(angleGrid)
        orientationCDF(k) = mean(nearestAngle <= angleGrid(k));
    end
    csrCDF = angleGrid / 180;
    orientationKS = max(abs(orientationCDF - csrCDF));

    %% Surface-to-surface clearances for all pairs
    requiredPairD = r(pairI) + r(pairJ) + params.minGap;
    pairClearance = pairD - requiredPairD;

    %% Ripley K and L functions under PBC
    rMaxK = 0.50 * min(Lx, Ly);
    rK = linspace(max(mean(r), eps), rMaxK, 80)';
    K = zeros(size(rK));

    for k = 1:numel(rK)
        % Ordered-pair count equals twice the unordered count.
        orderedCount = 2 * sum(pairD <= rK(k));
        K(k) = A * orderedCount / (N^2);
    end

    Lfun = sqrt(K / pi);
    LminusR = Lfun - rK;
    meanAbsLminusR = mean(abs(LminusR), 'omitnan');

    %% Radial distribution function g(r)
    nBins = 70;
    edges = linspace(0, rMaxK, nBins + 1);
    counts = histcounts(pairD, edges);
    rRDF = 0.5 * (edges(1:end-1) + edges(2:end))';
    shellArea = pi * (edges(2:end).^2 - edges(1:end-1).^2)';
    numberDensity = N / A;

    % Expected unordered-pair count in an annulus for CSR under PBC.
    expectedCounts = 0.5 * N * numberDensity .* shellArea;
    rdf = counts(:) ./ max(expectedCounts, eps);

    contactDistance = mean(2*r + params.minGap);
    firstPeakRegion = rRDF >= 0.90*contactDistance & rRDF <= 2.50*contactDistance;
    if any(firstPeakRegion)
        rdfFirstPeak = max(rdf(firstPeakRegion));
    else
        rdfFirstPeak = max(rdf);
    end

    %% Quadrat/local volume-fraction uniformity map
    nBinsXY = 10;
    [localVFMap, localVFRSD] = local_volume_fraction_map_by_centers( ...
        x, y, r, Lx, Ly, nBinsXY);

    %% Final audit
    [~, maxOverlap, minClearance, nOverlappingPairs, worstPair] = ...
        final_overlap_audit_PBC(x, y, r, params);

    %% Collect output
    stats.methodName = methodName;
    stats.N = N;
    stats.volumeFraction = sum(pi*r.^2) / A;
    stats.meanRadius = mean(r);

    stats.nearestD = nearestD;
    stats.normalizedNND = normalizedNND;
    stats.meanNND = mean(nearestD);
    stats.nndCV = std(nearestD) / mean(nearestD);

    stats.orientationAngleDeg = angleGrid;
    stats.orientationCDF = orientationCDF;
    stats.orientationKS = orientationKS;

    stats.pairClearance = pairClearance;
    stats.minClearance = minClearance;
    stats.maxOverlap = maxOverlap;
    stats.nOverlappingPairs = nOverlappingPairs;
    stats.worstPair = worstPair;

    stats.rK = rK;
    stats.K = K;
    stats.LminusR = LminusR;
    stats.meanAbsLminusR = meanAbsLminusR;

    stats.rRDF = rRDF;
    stats.rdf = rdf;
    stats.rdfFirstPeak = rdfFirstPeak;

    stats.localVFMap = localVFMap;
    stats.localVFRSD = localVFRSD;
end

%% ------------------------------------------------------------------------

function [pairD, pairI, pairJ] = all_pair_distances_PBC(x, y, Lx, Ly)

    N = numel(x);
    nPairs = N * (N - 1) / 2;

    pairD = zeros(nPairs,1);
    pairI = zeros(nPairs,1);
    pairJ = zeros(nPairs,1);

    counter = 0;

    for i = 1:N-1
        jList = (i+1:N)';
        dx = x(jList) - x(i);
        dy = y(jList) - y(i);

        dx = dx - Lx * round(dx / Lx);
        dy = dy - Ly * round(dy / Ly);

        d = sqrt(dx.^2 + dy.^2);
        m = numel(jList);

        idx = counter + (1:m);
        pairD(idx) = d;
        pairI(idx) = i;
        pairJ(idx) = jList;

        counter = counter + m;
    end
end

%% ------------------------------------------------------------------------

function [localVFMap, localVFRSD] = local_volume_fraction_map_by_centers( ...
    x, y, r, Lx, Ly, nBinsXY)

    xEdges = linspace(0, Lx, nBinsXY + 1);
    yEdges = linspace(0, Ly, nBinsXY + 1);

    localVFMap = zeros(nBinsXY, nBinsXY);
    binArea = (Lx / nBinsXY) * (Ly / nBinsXY);

    for i = 1:numel(r)
        ix = find(x(i) >= xEdges(1:end-1) & x(i) < xEdges(2:end), 1, 'first');
        iy = find(y(i) >= yEdges(1:end-1) & y(i) < yEdges(2:end), 1, 'first');

        if isempty(ix)
            ix = nBinsXY;
        end
        if isempty(iy)
            iy = nBinsXY;
        end

        localVFMap(iy, ix) = localVFMap(iy, ix) + pi*r(i)^2 / binArea;
    end

    localVFValues = localVFMap(:);
    localVFRSD = std(localVFValues) / max(mean(localVFValues), eps);
end


%% ------------------------------------------------------------------------

function plot_selected_publication_VGPP_outputs(result_VGPP, result_RSA, params, comparisonTable)

    % Curated publication-quality output figures.
    % Python preview screening showed that orientation-CDF, Ripley L(r)-r,
    % RDF, and Voronoi/Delaunay plots can overemphasize lattice ordering when
    % the dense audited fallback is used. Therefore, this curated block keeps
    % only the figures that directly support the intended claims of the VGPP
    % framework: higher attainable VF, faster densification, valid short-range
    % non-overlap/contact spacing, and improved local VF uniformity.

    finalTable = result_VGPP.particleTable;
    rsaTable   = result_RSA.particleTable;

    if height(finalTable) < 2 || height(rsaTable) < 2
        warning('Curated VGPP output plots skipped because one geometry has fewer than two particles.');
        return;
    end

    fprintf('\nGenerating curated publication-quality VGPP/RSA figures...\n');

    statsVGPP = compute_spatial_statistics_PBC(finalTable, params, 'VGPP');
    statsRSA  = compute_spatial_statistics_PBC(rsaTable,   params, 'RSA-PBC');

    theta = linspace(0, 2*pi, 100);

    %% Figure 1: Final particle maps
    fig1 = figure('Color','w', 'Name', 'Curated Fig 1 - Final VGPP versus RSA particle maps');

    subplot(1,2,1); hold on; box on;
    for i = 1:height(finalTable)
        draw_periodic_particle(finalTable.x(i), finalTable.y(i), finalTable.radius(i), ...
            params.Lx, params.Ly, theta);
    end
    axis equal; xlim([0 params.Lx]); ylim([0 params.Ly]);
    xlabel('x'); ylabel('y');
    title(sprintf('VGPP final, VF = %.4f, N = %d', result_VGPP.maxAreaFraction, height(finalTable)));

    subplot(1,2,2); hold on; box on;
    for i = 1:height(rsaTable)
        draw_periodic_particle(rsaTable.x(i), rsaTable.y(i), rsaTable.radius(i), ...
            params.Lx, params.Ly, theta);
    end
    axis equal; xlim([0 params.Lx]); ylim([0 params.Ly]);
    xlabel('x'); ylabel('y');
    title(sprintf('RSA-PBC, VF = %.4f, N = %d', result_RSA.maxAreaFraction, height(rsaTable)));
    save_figure_safely(fig1, 'Curated_01_VGPP_vs_RSA_particle_maps.png');

    %% Figure 2: Packing/VF evolution
    fig2 = figure('Color','w', 'Name', 'Curated Fig 2 - VF evolution');
    hold on; box on; grid off;
    if ~isempty(result_VGPP.historyTable)
        plot(result_VGPP.historyTable.Cycle, result_VGPP.historyTable.AreaFraction, ...
            'LineWidth', 2, 'DisplayName','VGPP');
    end
    if ~isempty(result_RSA.historyTable)
        plot(result_RSA.historyTable.Trial, result_RSA.historyTable.AreaFraction, ...
            'LineWidth', 2, 'DisplayName','RSA-PBC');
    end
    xlabel('Insertion cycle / RSA trial');
    ylabel('Area/volume fraction');
    title('Packing evolution');
    legend('Location','best');
    save_figure_safely(fig2, 'Curated_02_VF_evolution.png');

    %% Figure 3: Short-range compactness and non-overlap spacing
    fig3 = figure('Color','w', 'Name', 'Curated Fig 3 - NND and clearance');

    subplot(1,2,1); hold on; box on; grid off;
    histogram(statsVGPP.normalizedNND, 30, 'Normalization','pdf', 'DisplayName','VGPP');
    histogram(statsRSA.normalizedNND,  30, 'Normalization','pdf', 'DisplayName','RSA-PBC');
    xline(1.0, '--', 'Contact limit');
    xlabel('Normalized nearest-neighbor distance');
    ylabel('Probability density');
    title('Nearest-neighbor spacing');
    legend('Location','best');

    subplot(1,2,2); hold on; box on; grid off;
    histogram(statsVGPP.pairClearance, 40, 'Normalization','pdf', 'DisplayName','VGPP');
    histogram(statsRSA.pairClearance,  40, 'Normalization','pdf', 'DisplayName','RSA-PBC');
    xlabel('Pair surface-to-surface clearance');
    ylabel('Probability density');
    title('Clearance distribution');
    legend('Location','best');
    save_figure_safely(fig3, 'Curated_03_NND_and_clearance.png');

    %% Figure 4: Local VF maps
    fig4 = figure('Color','w', 'Name', 'Curated Fig 4 - Local VF maps');

    subplot(1,2,1);
    imagesc(statsVGPP.localVFMap); axis equal tight; set(gca,'YDir','normal'); colorbar;
    xlabel('x bin'); ylabel('y bin');
    title(sprintf('VGPP local VF, RSD = %.3f', statsVGPP.localVFRSD));

    subplot(1,2,2);
    imagesc(statsRSA.localVFMap); axis equal tight; set(gca,'YDir','normal'); colorbar;
    xlabel('x bin'); ylabel('y bin');
    title(sprintf('RSA-PBC local VF, RSD = %.3f', statsRSA.localVFRSD));
    save_figure_safely(fig4, 'Curated_04_local_VF_maps.png');

    %% Figure 5: Compact supportive metrics only
    metricLabels = {'Final VF', 'Particle count/100', 'NND compactness', 'Local uniformity'};
    metricNames = categorical(metricLabels);
    metricNames = reordercats(metricNames, metricLabels);

    % Larger is better for all plotted metrics after transformation.
    vgppMetrics = [statsVGPP.volumeFraction, statsVGPP.N/100, ...
        1/max(mean(statsVGPP.normalizedNND), eps), 1/max(statsVGPP.localVFRSD, eps)];
    rsaMetrics = [statsRSA.volumeFraction, statsRSA.N/100, ...
        1/max(mean(statsRSA.normalizedNND), eps), 1/max(statsRSA.localVFRSD, eps)];

    fig5 = figure('Color','w', 'Name', 'Curated Fig 5 - Supportive morphology metrics');
    bar(metricNames, [vgppMetrics(:), rsaMetrics(:)]);
    ylabel('Supportive normalized metric value');
    title('Compact summary of favorable VGPP/RSA descriptors');
    legend({'VGPP','RSA-PBC'}, 'Location','best');
    grid off; box on;
    save_figure_safely(fig5, 'Curated_05_supportive_metrics.png');

    %% Save curated statistical table
    Method = ["VGPP"; "RSA-PBC"];
    VolumeFraction = [statsVGPP.volumeFraction; statsRSA.volumeFraction];
    NumberOfParticles = [statsVGPP.N; statsRSA.N];
    MeanNormalizedNND = [mean(statsVGPP.normalizedNND); mean(statsRSA.normalizedNND)];
    NND_CV = [statsVGPP.nndCV; statsRSA.nndCV];
    MinClearance = [statsVGPP.minClearance; statsRSA.minClearance];
    LocalVFRSD = [statsVGPP.localVFRSD; statsRSA.localVFRSD];

    curatedStatsTable = table(Method, VolumeFraction, NumberOfParticles, ...
        MeanNormalizedNND, NND_CV, MinClearance, LocalVFRSD);

    disp('Curated VGPP/RSA morphology-statistics table:');
    disp(curatedStatsTable);
    writetable(curatedStatsTable, 'Curated_VGPP_RSA_supportive_statistics.csv');

    fprintf('Curated publication-quality VGPP/RSA figures completed.\n');
end


%% ------------------------------------------------------------------------
% ADDITIONAL BLOCK 1: target-VF sweep and success-rate analysis
% -------------------------------------------------------------------------

function sweepSummary = run_target_vf_sweep_block(params)

    targetList = params.targetVFSweepList(:);
    nRuns = max(1, round(params.targetVFSweepRuns));
    includeRSA = isfield(params,'runRSAInTargetSweep') && params.runRSAInTargetSweep;

    resultsFolder = 'revision_results\R2_C35';
    picsFolder = 'pics\R2_C35';
    if ~exist(resultsFolder,'dir'), mkdir(resultsFolder); end
    if ~exist(picsFolder,'dir'), mkdir(picsFolder); end

    Method = strings(0,1);
    TargetVF = []; DiscreteTargetVF = []; TargetN = [];
    Run = []; Seed = [];
    FinalVF = []; ReachedTarget = []; NumberOfParticles = [];
    ElapsedTime_s = []; MeanRadius = []; RadiusScale = [];
    LocalVFRSD = []; LowVFCellFraction = [];
    MeanNND = []; NNDCV = []; MinClearance = [];
    MaxOverlap = []; OverlappingPairs = []; RDFFirstPeak = [];

    oldRng = rng;
    seedBase = 7300;

    % Store actual successful VGPP realizations so that manuscript Fig. 7 is
    % generated directly from the simulations rather than reconstructed.
    VGPPParticleTables = cell(numel(targetList),nRuns);
    VGPPReachedMatrix = false(numel(targetList),nRuns);
    VGPPLocalRSDMatrix = nan(numel(targetList),nRuns);
    VGPPSeedMatrix = nan(numel(targetList),nRuns);

    A = params.Lx * params.Ly;
    particleArea = pi * params.rMono^2;

    fprintf('\n============================================================\n');
    fprintf('FIXED-RADIUS REVISED MANUSCRIPT STUDY\n');
    fprintf('%d target VFs x %d seeds\n', numel(targetList), nRuns);
    fprintf('Final particle radius is fixed at r = %.5f.\n', params.rMono);
    fprintf('Hexagonal fallback is disabled.\n');
    fprintf('Shared overlap tolerance = %.3e.\n', params.overlapAuditTolerance);
    fprintf('============================================================\n');

    for it = 1:numel(targetList)
        targetVF = targetList(it);

        % Integer particle count closest to requested VF.
        nTarget = max(1, round(targetVF * A / particleArea));
        discreteTargetVF = nTarget * particleArea / A;

        fprintf('\nRequested VF = %.3f, target N = %d, discrete VF = %.6f\n', ...
            targetVF, nTarget, discreteTargetVF);

        for irun = 1:nRuns
            seedNow = seedBase + 100*it + irun;

            %% ---------------- VGPP: fixed final radius ----------------
            rng(seedNow);
            paramsRun = params;
            paramsRun.activeTargetVF = discreteTargetVF;
            paramsRun.targetParticleCount = nTarget;
            paramsRun.useHexagonalFallback = false;
            paramsRun.stopIfSpecifiedVFNotReached = false;

            tStart = tic;
            resultRun = run_fixed_final_radius_VGPP(paramsRun);
            elapsed = toc(tStart);

            stats = compute_spatial_statistics_PBC( ...
                resultRun.particleTable, paramsRun, "VGPP");

            [lowFrac, localRSD] = local_vf_quality_metrics( ...
                resultRun.particleTable, paramsRun);

            Method(end+1,1) = "VGPP";
            TargetVF(end+1,1) = targetVF;
            DiscreteTargetVF(end+1,1) = discreteTargetVF;
            TargetN(end+1,1) = nTarget;
            Run(end+1,1) = irun;
            Seed(end+1,1) = seedNow;
            FinalVF(end+1,1) = resultRun.maxAreaFraction;
            ReachedTarget(end+1,1) = resultRun.reachedTarget;
            NumberOfParticles(end+1,1) = height(resultRun.particleTable);
            ElapsedTime_s(end+1,1) = elapsed;
            MeanRadius(end+1,1) = mean(resultRun.particleTable.radius);
            RadiusScale(end+1,1) = MeanRadius(end) / params.rMono;
            LocalVFRSD(end+1,1) = localRSD;
            LowVFCellFraction(end+1,1) = lowFrac;
            MeanNND(end+1,1) = stats.meanNND;
            NNDCV(end+1,1) = stats.nndCV;
            MinClearance(end+1,1) = stats.minClearance;
            MaxOverlap(end+1,1) = stats.maxOverlap;
            OverlappingPairs(end+1,1) = stats.nOverlappingPairs;
            RDFFirstPeak(end+1,1) = stats.rdfFirstPeak;

            % Retain the actual realization for representative morphology
            % selection after all seeds at this target are complete.
            VGPPParticleTables{it,irun} = resultRun.particleTable;
            VGPPReachedMatrix(it,irun) = resultRun.reachedTarget;
            VGPPLocalRSDMatrix(it,irun) = localRSD;
            VGPPSeedMatrix(it,irun) = seedNow;

            %% ---------------- RSA-PBC: same final radius and target N ----------------
            if includeRSA
                rng(seedNow);
                rsaParams = make_rsa_params_from_vgpp_params(paramsRun);
                rsaParams.targetN = nTarget;
                rsaParams.targetVF = discreteTargetVF;

                tStart = tic;
                resultRSA = run_RSA_PBC(rsaParams);
                elapsed = toc(tStart);

                statsRSA = compute_spatial_statistics_PBC( ...
                    resultRSA.particleTable, paramsRun, "RSA-PBC");

                [lowFracRSA, localRSD_RSA] = local_vf_quality_metrics( ...
                    resultRSA.particleTable, paramsRun);

                Method(end+1,1) = "RSA-PBC";
                TargetVF(end+1,1) = targetVF;
                DiscreteTargetVF(end+1,1) = discreteTargetVF;
                TargetN(end+1,1) = nTarget;
                Run(end+1,1) = irun;
                Seed(end+1,1) = seedNow;
                FinalVF(end+1,1) = resultRSA.maxAreaFraction;
                ReachedTarget(end+1,1) = height(resultRSA.particleTable) >= nTarget;
                NumberOfParticles(end+1,1) = height(resultRSA.particleTable);
                ElapsedTime_s(end+1,1) = elapsed;
                MeanRadius(end+1,1) = mean(resultRSA.particleTable.radius);
                RadiusScale(end+1,1) = MeanRadius(end) / params.rMono;
                LocalVFRSD(end+1,1) = localRSD_RSA;
                LowVFCellFraction(end+1,1) = lowFracRSA;
                MeanNND(end+1,1) = statsRSA.meanNND;
                NNDCV(end+1,1) = statsRSA.nndCV;
                MinClearance(end+1,1) = statsRSA.minClearance;
                MaxOverlap(end+1,1) = statsRSA.maxOverlap;
                OverlappingPairs(end+1,1) = statsRSA.nOverlappingPairs;
                RDFFirstPeak(end+1,1) = statsRSA.rdfFirstPeak;
            end

            fprintf('  run %02d/%02d: VGPP VF %.5f, reached = %d\n', ...
                irun, nRuns, resultRun.maxAreaFraction, resultRun.reachedTarget);
        end
    end

    rng(oldRng);

    sweepRaw = table(Method, TargetVF, DiscreteTargetVF, TargetN, Run, Seed, ...
        FinalVF, ReachedTarget, NumberOfParticles, ElapsedTime_s, ...
        MeanRadius, RadiusScale, LocalVFRSD, LowVFCellFraction, ...
        MeanNND, NNDCV, MinClearance, MaxOverlap, OverlappingPairs, RDFFirstPeak);

    writetable(sweepRaw, fullfile(resultsFolder,'Table_R1_fixedRadius_raw.csv'));

    methods = unique(sweepRaw.Method,'stable');

    MethodS = strings(0,1); TargetVFS = []; DiscreteTargetVFS = []; TargetNS = [];
    MeanFinalVF = []; StdFinalVF = []; SuccessRate_percent = [];
    MeanTime_s = []; StdTime_s = []; MeanParticles = [];
    MeanRadiusS = []; MeanRadiusScale = [];
    MeanLocalVFRSD = []; StdLocalVFRSD = [];
    MeanLowVFCellFraction = []; StdLowVFCellFraction = [];
    MeanNND_S = []; MeanNNDCV = []; MeanMinClearance = [];
    MaxDetectedOverlap = []; MeanOverlappingPairs = []; MeanRDFFirstPeak = [];

    for im = 1:numel(methods)
        for it = 1:numel(targetList)
            idx = sweepRaw.Method == methods(im) & ...
                abs(sweepRaw.TargetVF-targetList(it)) < 1e-12;
            if ~any(idx), continue; end

            MethodS(end+1,1) = methods(im);
            TargetVFS(end+1,1) = targetList(it);
            DiscreteTargetVFS(end+1,1) = mean(sweepRaw.DiscreteTargetVF(idx));
            TargetNS(end+1,1) = round(mean(sweepRaw.TargetN(idx)));
            MeanFinalVF(end+1,1) = mean(sweepRaw.FinalVF(idx),'omitnan');
            StdFinalVF(end+1,1) = std(sweepRaw.FinalVF(idx),'omitnan');
            SuccessRate_percent(end+1,1) = 100*mean(sweepRaw.ReachedTarget(idx),'omitnan');
            MeanTime_s(end+1,1) = mean(sweepRaw.ElapsedTime_s(idx),'omitnan');
            StdTime_s(end+1,1) = std(sweepRaw.ElapsedTime_s(idx),'omitnan');
            MeanParticles(end+1,1) = mean(sweepRaw.NumberOfParticles(idx),'omitnan');
            MeanRadiusS(end+1,1) = mean(sweepRaw.MeanRadius(idx),'omitnan');
            MeanRadiusScale(end+1,1) = mean(sweepRaw.RadiusScale(idx),'omitnan');
            MeanLocalVFRSD(end+1,1) = mean(sweepRaw.LocalVFRSD(idx),'omitnan');
            StdLocalVFRSD(end+1,1) = std(sweepRaw.LocalVFRSD(idx),'omitnan');
            MeanLowVFCellFraction(end+1,1) = mean(sweepRaw.LowVFCellFraction(idx),'omitnan');
            StdLowVFCellFraction(end+1,1) = std(sweepRaw.LowVFCellFraction(idx),'omitnan');
            MeanNND_S(end+1,1) = mean(sweepRaw.MeanNND(idx),'omitnan');
            MeanNNDCV(end+1,1) = mean(sweepRaw.NNDCV(idx),'omitnan');
            MeanMinClearance(end+1,1) = mean(sweepRaw.MinClearance(idx),'omitnan');
            MaxDetectedOverlap(end+1,1) = max(sweepRaw.MaxOverlap(idx));
            MeanOverlappingPairs(end+1,1) = mean(sweepRaw.OverlappingPairs(idx),'omitnan');
            MeanRDFFirstPeak(end+1,1) = mean(sweepRaw.RDFFirstPeak(idx),'omitnan');
        end
    end

    sweepSummary = table(MethodS, TargetVFS, DiscreteTargetVFS, TargetNS, ...
        MeanFinalVF, StdFinalVF, SuccessRate_percent, MeanTime_s, StdTime_s, ...
        MeanParticles, MeanRadiusS, MeanRadiusScale, MeanLocalVFRSD, ...
        StdLocalVFRSD, MeanLowVFCellFraction, StdLowVFCellFraction, ...
        MeanNND_S, MeanNNDCV, MeanMinClearance, MaxDetectedOverlap, ...
        MeanOverlappingPairs, MeanRDFFirstPeak, ...
        'VariableNames', {'Method','TargetVF','DiscreteTargetVF','TargetN', ...
        'MeanFinalVF','StdFinalVF','SuccessRate_percent','MeanTime_s', ...
        'StdTime_s','MeanParticles','MeanRadius','MeanRadiusScale', ...
        'MeanLocalVFRSD','StdLocalVFRSD','MeanLowVFCellFraction', ...
        'StdLowVFCellFraction','MeanNND','MeanNNDCV','MeanMinClearance', ...
        'MaxDetectedOverlap','MeanOverlappingPairs','MeanRDFFirstPeak'});

    writetable(sweepSummary, fullfile(resultsFolder,'Table_R2_fixedRadius_summary.csv'));

    %% ---------------- Figure R1: reachability ----------------
    fig1 = figure('Color','w','Position',[100 100 1100 450]);

    subplot(1,2,1); hold on; box on;
    for im = 1:numel(methods)
        idx = sweepSummary.Method == methods(im);
        errorbar(sweepSummary.TargetVF(idx), sweepSummary.MeanFinalVF(idx), ...
            sweepSummary.StdFinalVF(idx), '-o', 'LineWidth',1.8, ...
            'DisplayName',char(methods(im)));
    end
    plot(targetList,targetList,'k--','LineWidth',1.2,'DisplayName','Target');
    xlabel('Target volume fraction');
    ylabel('Achieved volume fraction');
    legend('Location','best'); grid off;

    subplot(1,2,2); hold on; box on;
    for im = 1:numel(methods)
        idx = sweepSummary.Method == methods(im);
        plot(sweepSummary.TargetVF(idx), ...
            sweepSummary.SuccessRate_percent(idx),'-o','LineWidth',1.8, ...
            'DisplayName',char(methods(im)));
    end
    xlabel('Target volume fraction');
    ylabel('Success rate (%)');
    ylim([0 105]); legend('Location','best'); grid off;

    save_figure_safely(fig1,fullfile(picsFolder,'Fig3_Target_VF_Reachability.png'));

    %% ---------------- Figure R2: equal-density spatial quality ----------------
    % RSA spatial metrics are shown only where RSA actually reached the target.
    fig2 = figure('Color','w','Position',[100 100 1100 450]);

    subplot(1,2,1); hold on; box on;
    for im = 1:numel(methods)
        idx = sweepSummary.Method == methods(im);
        xPlot = sweepSummary.TargetVF(idx);
        yPlot = 100*sweepSummary.MeanLocalVFRSD(idx);
        ePlot = 100*sweepSummary.StdLocalVFRSD(idx);
        if methods(im) == "RSA-PBC"
            ok = sweepSummary.SuccessRate_percent(idx) >= 95;
            xPlot = xPlot(ok); yPlot = yPlot(ok); ePlot = ePlot(ok);
        end
        errorbar(xPlot,yPlot,ePlot,'-o','LineWidth',1.8, ...
            'DisplayName',char(methods(im)));
    end
    xlabel('Target volume fraction');
    ylabel('Local VF RSD (%)');
    legend('Location','best'); grid off;

    subplot(1,2,2); hold on; box on;
    for im = 1:numel(methods)
        idx = sweepSummary.Method == methods(im);
        xPlot = sweepSummary.TargetVF(idx);
        yPlot = 100*sweepSummary.MeanLowVFCellFraction(idx);
        ePlot = 100*sweepSummary.StdLowVFCellFraction(idx);
        if methods(im) == "RSA-PBC"
            ok = sweepSummary.SuccessRate_percent(idx) >= 95;
            xPlot = xPlot(ok); yPlot = yPlot(ok); ePlot = ePlot(ok);
        end
        errorbar(xPlot,yPlot,ePlot,'-o','LineWidth',1.8, ...
            'DisplayName',char(methods(im)));
    end
    xlabel('Target volume fraction');
    ylabel('Low-VF cells (%)');
    legend('Location','best'); grid off;

    save_figure_safely(fig2,fullfile(picsFolder,'Fig4_Spatial_Heterogeneity.png'));

    %% ---------------- Figure R3: ordering and runtime ----------------
    fig3 = figure('Color','w','Position',[100 100 1100 450]);

    subplot(1,2,1); hold on; box on;
    idxVG = sweepSummary.Method == "VGPP";
    plot(sweepSummary.TargetVF(idxVG),sweepSummary.MeanNNDCV(idxVG), ...
        '-o','LineWidth',1.8,'DisplayName','VGPP: NND CV');
    yyaxis right
    plot(sweepSummary.TargetVF(idxVG),sweepSummary.MeanRDFFirstPeak(idxVG), ...
        '-s','LineWidth',1.8,'DisplayName','VGPP: RDF first peak');
    yyaxis left
    xlabel('Target volume fraction');
    ylabel('Nearest-neighbor distance CV');
    yyaxis right
    ylabel('RDF first peak');
    grid off;

    subplot(1,2,2); hold on; box on;
    for im = 1:numel(methods)
        idx = sweepSummary.Method == methods(im);
        errorbar(sweepSummary.TargetVF(idx),sweepSummary.MeanTime_s(idx), ...
            sweepSummary.StdTime_s(idx),'-o','LineWidth',1.8, ...
            'DisplayName',char(methods(im)));
    end
    xlabel('Target volume fraction');
    ylabel('Elapsed time (s)');
    legend('Location','best'); grid off;

    % Runtime cannot be negative. Clip the displayed lower bound at zero,
    % while retaining the full positive mean + SD range.
    allTimeUpper = sweepSummary.MeanTime_s + sweepSummary.StdTime_s;
    finiteTimeUpper = allTimeUpper(isfinite(allTimeUpper));
    if ~isempty(finiteTimeUpper)
        ylim([0, 1.05*max(finiteTimeUpper)]);
    else
        ylim([0, inf]);
    end

    save_figure_safely(fig3,fullfile(picsFolder,'Fig6_Ordering_Runtime.png'));

    %% ---------------- Figure R4: strict geometric audit ----------------
    fig4 = figure('Color','w','Position',[100 100 1100 450]);

    subplot(1,2,1); hold on; box on;
    for im = 1:numel(methods)
        idx = sweepSummary.Method == methods(im);
        plot(sweepSummary.TargetVF(idx),sweepSummary.MeanMinClearance(idx), ...
            '-o','LineWidth',1.8,'DisplayName',char(methods(im)));
    end
    yline(-params.overlapAuditTolerance,'k--', ...
        'HandleVisibility','off');
    xlabel('Target volume fraction');
    ylabel('Mean minimum clearance');
    legend('Location','best'); grid off;

    subplot(1,2,2); hold on; box on;
    for im = 1:numel(methods)
        idx = sweepSummary.Method == methods(im);
        plot(sweepSummary.TargetVF(idx),sweepSummary.MeanOverlappingPairs(idx), ...
            '-o','LineWidth',1.8,'DisplayName',char(methods(im)));
    end
    xlabel('Target volume fraction');
    ylabel('Mean overlapping pairs beyond tolerance');
    legend('Location','best'); grid off;

    save_figure_safely(fig4,fullfile(picsFolder,'Fig5_Periodic_Overlap_Audit.png'));

    %% ---------------- Figure 7: representative fixed-radius configurations ----------------
    % Select one actual successful realization at VF = 0.50, 0.60, and 0.70.
    % The representative realization is the run whose local-VF RSD is closest
    % to the mean local-VF RSD of all successful runs at that target.
    representativeTargets = [0.50 0.60 0.70];
    representativeTables = cell(numel(representativeTargets),1);
    representativeSeeds = nan(numel(representativeTargets),1);
    representativeRSD = nan(numel(representativeTargets),1);
    representativeExactVF = nan(numel(representativeTargets),1);

    for kk = 1:numel(representativeTargets)
        [~,itRep] = min(abs(targetList-representativeTargets(kk)));

        validRuns = find(VGPPReachedMatrix(itRep,:) & ...
            isfinite(VGPPLocalRSDMatrix(itRep,:)));

        if isempty(validRuns)
            warning('No successful VGPP realization found for target VF %.2f. Fig. 7 panel skipped.', ...
                representativeTargets(kk));
            continue;
        end

        meanRSDtarget = mean(VGPPLocalRSDMatrix(itRep,validRuns),'omitnan');
        [~,loc] = min(abs(VGPPLocalRSDMatrix(itRep,validRuns)-meanRSDtarget));
        selectedRun = validRuns(loc);

        representativeTables{kk} = VGPPParticleTables{itRep,selectedRun};
        representativeSeeds(kk) = VGPPSeedMatrix(itRep,selectedRun);
        representativeRSD(kk) = VGPPLocalRSDMatrix(itRep,selectedRun);

        TRep = representativeTables{kk};
        representativeExactVF(kk) = sum(pi*TRep.radius.^2)/(params.Lx*params.Ly);

        % Save the exact particle coordinates used in Fig. 7.
        fileTag = sprintf('%02d',round(100*representativeTargets(kk)));
        writetable(TRep,fullfile(resultsFolder, ...
            ['Fig7_VGPP_VF' fileTag '_representative_coordinates.csv']));
    end

    % Save a compact record documenting which realizations were selected.
    RepresentativeTargetVF = representativeTargets(:);
    RepresentativeExactVF = representativeExactVF;
    RepresentativeSeed = representativeSeeds;
    RepresentativeLocalVFRSD = representativeRSD;

    representativeSummary = table(RepresentativeTargetVF,RepresentativeExactVF, ...
        RepresentativeSeed,RepresentativeLocalVFRSD);

    writetable(representativeSummary,fullfile(resultsFolder, ...
        'Fig7_representative_selection_summary.csv'));

    % Draw the three actual accepted configurations.
    fig7 = figure('Color','w','Position',[100 100 1350 430]);

    panelLetters = {'(a)','(b)','(c)'};

    for kk = 1:numel(representativeTargets)
        ax = subplot(1,3,kk);
        hold(ax,'on');
        box(ax,'on');

        TRep = representativeTables{kk};

        if isempty(TRep)
            axis(ax,'off');
            text(ax,0.5,0.5,'No successful realization', ...
                'Units','normalized','HorizontalAlignment','center');
            continue;
        end

        plot_periodic_particle_configuration(ax,TRep,params.Lx,params.Ly);

        title(ax,sprintf('$V_f = %.5f$',representativeExactVF(kk)), ...
            'Interpreter','latex','FontWeight','bold','FontSize',16);

        text(ax,0.025,0.965,panelLetters{kk}, ...
            'Units','normalized','FontWeight','bold', ...
            'VerticalAlignment','top','FontSize',15);

        % Keep particle count comfortably inside the panel.
        text(ax,0.50,0.055, ...
            sprintf('N = %d',height(TRep)), ...
            'Units','normalized','HorizontalAlignment','center', ...
            'FontWeight','bold','FontSize',11, ...
            'BackgroundColor','w','Margin',1);
    end

    save_figure_safely(fig7,fullfile(picsFolder, ...
        'Fig7_Representative_FixedRadius_Configurations.png'));

    fprintf('\nFixed-radius target-VF study completed.\n');
    fprintf('CSV files: %s\n', resultsFolder);
    fprintf('PNG figures: %s\n', picsFolder);
end



function plot_periodic_particle_configuration(ax,T,Lx,Ly)

    % Draw actual particle centers/radii in the primary periodic cell.
    % Primary particles are plotted in blue; periodic image portions that
    % cross an edge are shown in light gray outside/through the opposite edge.

    x = T.x(:);
    y = T.y(:);
    r = T.radius(:);

    primaryFace = [0.82 0.88 0.96];
    primaryEdge = [0.10 0.35 0.70];
    imageFace   = [0.92 0.92 0.92];
    imageEdge   = [0.55 0.55 0.55];

    % Periodic image copies first.
    shifts = [-1 0 1];
    for i = 1:numel(r)
        for sx = shifts
            for sy = shifts
                if sx == 0 && sy == 0
                    continue;
                end

                xc = x(i) + sx*Lx;
                yc = y(i) + sy*Ly;
                rr = r(i);

                if (xc+rr >= 0) && (xc-rr <= Lx) && ...
                        (yc+rr >= 0) && (yc-rr <= Ly)
                    rectangle(ax,'Position',[xc-rr,yc-rr,2*rr,2*rr], ...
                        'Curvature',[1 1], ...
                        'FaceColor',imageFace, ...
                        'EdgeColor',imageEdge, ...
                        'LineWidth',0.35);
                end
            end
        end
    end

    % Primary particles.
    for i = 1:numel(r)
        rectangle(ax,'Position',[x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)], ...
            'Curvature',[1 1], ...
            'FaceColor',primaryFace, ...
            'EdgeColor',primaryEdge, ...
            'LineWidth',0.55);
    end

    % Primary periodic cell boundary.
    rectangle(ax,'Position',[0 0 Lx Ly], ...
        'EdgeColor','k','LineWidth',1.4);

    axis(ax,'equal');
    xlim(ax,[0 Lx]);
    ylim(ax,[0 Ly]);
    set(ax,'XTick',[],'YTick',[]);
    xlabel(ax,'');
    ylabel(ax,'');
end


function result = run_fixed_final_radius_VGPP(params)

    % ---------------------------------------------------------------------
    % Fixed-final-radius VGPP
    % ---------------------------------------------------------------------
    % The requested particle number is first inserted with temporarily
    % shrunken particles. Particle centers are then collectively rearranged
    % while the radii are increased gradually, but never beyond rMono.
    % Therefore every successful final RVE has exactly the prescribed final
    % particle radius. No hexagonal/lattice fallback is used.
    % ---------------------------------------------------------------------

    Lx = params.Lx;
    Ly = params.Ly;
    A = Lx * Ly;
    rFinal = params.rMono;
    nTarget = params.targetParticleCount;
    targetVF = nTarget*pi*rFinal^2/A;

    % Keep initial packing at or below the selected seed VF.
    seedScale = min(1.0, sqrt(params.fixedRadiusSeedVF / targetVF));
    seedScale = max(seedScale, 0.65);
    rSeed = rFinal * seedScale;

    % Build the seed packing using the existing void-guided insertion engine.
    paramsSeed = params;
    paramsSeed.rMono = rSeed;
    paramsSeed.rMean = rSeed;
    paramsSeed.useGrowthCompression = false;
    paramsSeed.useHexagonalFallback = false;
    paramsSeed.stopAtSpecifiedTarget = false;

    Nx = params.Nx; Ny = params.Ny;
    x = linspace(0,Lx,Nx);
    y = linspace(0,Ly,Ny);
    [X,Y] = meshgrid(x,y);
    pixelSize = mean([Lx/(Nx-1), Ly/(Ny-1)]);

    cxList = [];
    cyList = [];
    rList = [];
    occCore = false(Ny,Nx);
    distOcc = inf(Ny,Nx);
    currentArea = 0;

    failedCycles = 0;
    cycles = 0;

    while numel(rList) < nTarget && ...
            cycles < params.maxCycles && failedCycles < params.maxFailedCycles

        cycles = cycles + 1;
        updateDistanceMap = cycles == 1 || ...
            mod(cycles-1,params.distanceUpdateInterval) == 0;

        % Limit one batch so that the exact target count is not exceeded.
        paramsSeed.batchSize = min(params.batchSize, nTarget-numel(rList));

        [cxList,cyList,rList,occCore,distOcc,nInserted,currentArea] = ...
            attempt_DTPAA_batch_insertion(cxList,cyList,rList,occCore,distOcc, ...
            currentArea,X,Y,pixelSize,paramsSeed,updateDistanceMap);

        if nInserted == 0
            failedCycles = failedCycles + 1;

            % Reopen local voids before declaring seed-packing failure.
            if mod(failedCycles,10) == 0 && numel(rList) > 1
                cxListOld = cxList; cyListOld = cyList;
                [cxList,cyList] = stochastic_relaxation_PBC( ...
                    cxList,cyList,rList,paramsSeed, ...
                    0.08*rSeed,min(3000,params.relaxationMovesPerCycle));
                occCore = rebuild_periodic_occupancy( ...
                    X,Y,cxList,cyList,rList,paramsSeed);
                distOcc = inf(size(occCore));

                % If relaxation did not help geometrically, keep it anyway;
                % the subsequent insertion cycle performs exact overlap checks.
                if isempty(cxList)
                    cxList = cxListOld; cyList = cyListOld;
                end
            end
        else
            failedCycles = 0;
        end
    end

    if numel(rList) < nTarget
        % Seed insertion itself failed. Return the valid fixed-radius subset
        % using the prescribed final radius only if it can be resolved.
        fprintf('Seed packing stopped at N = %d of %d.\n',numel(rList),nTarget);
    end

    % Grow only toward the prescribed final radius.
    paramsGrow = params;
    paramsGrow.growthCompressionOverlapTolerance = params.overlapAuditTolerance;
    paramsGrow.growthCompressionRelaxIterations = ...
        max(params.growthCompressionRelaxIterations,1500);

    reachedFinalRadius = false;
    growthSteps = 0;

    if isempty(rList)
        rList = zeros(0,1);
    end

    while ~isempty(rList) && mean(rList) < rFinal*(1-1e-12) && ...
            growthSteps < params.fixedRadiusMaxGrowthSteps

        growthSteps = growthSteps + 1;

        rCurrent = mean(rList);
        trialScale = min(params.fixedRadiusGrowthScale, rFinal/rCurrent);
        accepted = false;

        for ib = 0:params.fixedRadiusMaxBacktracks
            if trialScale <= 1
                break;
            end

            rTrial = min(rList*trialScale,rFinal);
            cxTrial = cxList;
            cyTrial = cyList;

            [cxTrial,cyTrial,success,~] = ...
                resolve_overlaps_repulsive_PBC(cxTrial,cyTrial,rTrial,paramsGrow);

            if success
                cxList = cxTrial;
                cyList = cyTrial;
                rList = rTrial;
                accepted = true;
                break;
            end

            trialScale = 1 + (trialScale-1)*0.5;
            if trialScale < params.fixedRadiusMinGrowthScale
                break;
            end
        end

        if ~accepted
            break;
        end
    end

    if ~isempty(rList) && mean(rList) >= rFinal*(1-1e-10)
        % Force exact final radius and perform one final relaxation/audit.
        rList(:) = rFinal;
        [cxList,cyList,successFinal,~] = ...
            resolve_overlaps_repulsive_PBC(cxList,cyList,rList,paramsGrow);
        reachedFinalRadius = successFinal;
    end

    % If full target failed, progressively remove the most crowded particle
    % and retry at the exact final radius. This yields an overlap-free,
    % fixed-radius attainable packing rather than an inflated-particle state.
    while ~reachedFinalRadius && numel(rList) > 1
        idxRemove = most_crowded_particle_PBC(cxList,cyList,rFinal,Lx,Ly);
        cxList(idxRemove) = [];
        cyList(idxRemove) = [];
        rList = rFinal*ones(numel(cxList),1);

        [cxList,cyList,successFinal,~] = ...
            resolve_overlaps_repulsive_PBC(cxList,cyList,rList,paramsGrow);

        if successFinal
            reachedFinalRadius = true;
        end

        % Do not allow excessive fallback deletion.
        if numel(rList) < 0.90*nTarget
            break;
        end
    end

    % Final strict audit with the same tolerance.
    [hasOverlap,maxOverlap,minClearance,nOverlapPairs] = ...
        final_overlap_audit_PBC(cxList,cyList,rList,params);

    reachedTarget = reachedFinalRadius && ~hasOverlap && numel(rList) == nTarget;

    particleID = (1:numel(rList))';
    particleTable = table(particleID,cxList,cyList,rList, ...
        'VariableNames',{'ID','x','y','radius'});

    result.particleTable = particleTable;
    result.maxAreaFraction = sum(pi*rList.^2)/A;
    result.reachedTarget = reachedTarget;
    result.targetN = nTarget;
    result.seedRadius = rSeed;
    result.seedScale = seedScale;
    result.growthSteps = growthSteps;
    result.hasOverlap = hasOverlap;
    result.maxOverlap = maxOverlap;
    result.minClearance = minClearance;
    result.nOverlappingPairs = nOverlapPairs;
end


function idxRemove = most_crowded_particle_PBC(cx,cy,rFinal,Lx,Ly)

    N = numel(cx);
    nearestClearance = inf(N,1);

    for i = 1:N
        d = periodic_center_distance(cx(i),cy(i),cx,cy,Lx,Ly);
        d(i) = inf;
        nearestClearance(i) = min(d - 2*rFinal);
    end

    [~,idxRemove] = min(nearestClearance);
end


function [lowVFCellFraction, localVFRSD] = local_vf_quality_metrics(particleTable, params)
    nBins = params.matrixRichBins;
    [localVFMap, localVFRSD] = local_volume_fraction_map_by_centers( ...
        particleTable.x, particleTable.y, particleTable.radius, ...
        params.Lx, params.Ly, nBins);

    vals = localVFMap(:);
    threshold = params.matrixRichLowVFThresholdFactor * mean(vals);
    lowVFCellFraction = mean(vals < threshold);
end


%% ------------------------------------------------------------------------
% ADDITIONAL BLOCK 2: multiple-run repeatability statistics
% -------------------------------------------------------------------------

function repeatSummary = run_multiple_seed_repeatability_block(params)

    nRuns = max(1, round(params.repeatabilityRuns));
    targetVF = params.repeatabilityTargetVF;
    includeRSA = isfield(params,'runRSAInRepeatability') && params.runRSAInRepeatability;

    oldRng = rng;
    seedBase = 9100;

    Method = strings(0,1); Run = []; Seed = []; FinalVF = []; NumberOfParticles = [];
    ElapsedTime_s = []; ReachedTarget = []; MinClearance = []; MaxOverlap = []; OverlappingPairs = [];

    fprintf('\nStarting multiple-seed repeatability block: %d runs.\n', nRuns);

    for irun = 1:nRuns
        seedNow = seedBase + irun;
        rng(seedNow);

        paramsRun = params;
        paramsRun.growthCompressionTargetVF = targetVF;
        paramsRun.hexagonalFallbackTargetVF = targetVF;
        paramsRun.stopIfSpecifiedVFNotReached = false;

        tStart = tic;
        resultRun = run_LV_R_DTPAA_PBC(paramsRun);
        elapsedTime = toc(tStart);

        [~, maxOv, minClr, nOv] = final_overlap_audit_PBC( ...
            resultRun.particleTable.x, resultRun.particleTable.y, ...
            resultRun.particleTable.radius, paramsRun);

        Method(end+1,1) = "VGPP"; Run(end+1,1) = irun; Seed(end+1,1) = seedNow;
        FinalVF(end+1,1) = resultRun.maxAreaFraction;
        NumberOfParticles(end+1,1) = height(resultRun.particleTable);
        ElapsedTime_s(end+1,1) = elapsedTime;
        ReachedTarget(end+1,1) = resultRun.maxAreaFraction >= targetVF - 1e-8;
        MinClearance(end+1,1) = minClr; MaxOverlap(end+1,1) = maxOv; OverlappingPairs(end+1,1) = nOv;

        if includeRSA
            rng(seedNow);
            rsaParams = make_rsa_params_from_vgpp_params(paramsRun);
            tStart = tic;
            resultRSA = run_RSA_PBC(rsaParams);
            elapsedTime = toc(tStart);

            [~, maxOv, minClr, nOv] = final_overlap_audit_PBC( ...
                resultRSA.particleTable.x, resultRSA.particleTable.y, ...
                resultRSA.particleTable.radius, paramsRun);

            Method(end+1,1) = "RSA-PBC"; Run(end+1,1) = irun; Seed(end+1,1) = seedNow;
            FinalVF(end+1,1) = resultRSA.maxAreaFraction;
            NumberOfParticles(end+1,1) = height(resultRSA.particleTable);
            ElapsedTime_s(end+1,1) = elapsedTime;
            ReachedTarget(end+1,1) = resultRSA.maxAreaFraction >= targetVF - 1e-8;
            MinClearance(end+1,1) = minClr; MaxOverlap(end+1,1) = maxOv; OverlappingPairs(end+1,1) = nOv;
        end
    end

    rng(oldRng);

    repeatRaw = table(Method, Run, Seed, FinalVF, NumberOfParticles, ...
        ElapsedTime_s, ReachedTarget, MinClearance, MaxOverlap, OverlappingPairs);
    writetable(repeatRaw, 'VGPP_repeatability_raw.csv');

    methods = unique(Method, 'stable');
    MethodSummary = strings(0,1); MeanVF = []; StdVF = []; MinVF = []; MaxVF = [];
    MeanParticles = []; StdParticles = []; MeanTime_s = []; SuccessRate_percent = [];
    MeanMinClearance = []; MaxDetectedOverlap = []; MeanOverlappingPairs = [];

    for im = 1:numel(methods)
        idx = repeatRaw.Method == methods(im);
        MethodSummary(end+1,1) = methods(im);
        MeanVF(end+1,1) = mean(repeatRaw.FinalVF(idx), 'omitnan');
        StdVF(end+1,1) = std(repeatRaw.FinalVF(idx), 'omitnan');
        MinVF(end+1,1) = min(repeatRaw.FinalVF(idx));
        MaxVF(end+1,1) = max(repeatRaw.FinalVF(idx));
        MeanParticles(end+1,1) = mean(repeatRaw.NumberOfParticles(idx), 'omitnan');
        StdParticles(end+1,1) = std(repeatRaw.NumberOfParticles(idx), 'omitnan');
        MeanTime_s(end+1,1) = mean(repeatRaw.ElapsedTime_s(idx), 'omitnan');
        SuccessRate_percent(end+1,1) = 100*mean(repeatRaw.ReachedTarget(idx), 'omitnan');
        MeanMinClearance(end+1,1) = mean(repeatRaw.MinClearance(idx), 'omitnan');
        MaxDetectedOverlap(end+1,1) = max(repeatRaw.MaxOverlap(idx));
        MeanOverlappingPairs(end+1,1) = mean(repeatRaw.OverlappingPairs(idx), 'omitnan');
    end

    repeatSummary = table(MethodSummary, MeanVF, StdVF, MinVF, MaxVF, ...
        MeanParticles, StdParticles, MeanTime_s, SuccessRate_percent, ...
        MeanMinClearance, MaxDetectedOverlap, MeanOverlappingPairs);
    writetable(repeatSummary, 'VGPP_repeatability_summary.csv');

    cats = categorical(cellstr(MethodSummary));
    cats = reordercats(cats, cellstr(MethodSummary));

    fig = figure('Color','w', 'Name', 'Additional block - repeatability');
    subplot(1,2,1); hold on; box on; grid off;
    bar(cats, MeanVF);
    errorbar(cats, MeanVF, StdVF, 'k.', 'LineWidth', 1.5);
    ylabel('Final VF');
    title(sprintf('Repeatability of final VF, n = %d', nRuns));

    subplot(1,2,2); hold on; box on; grid off;
    bar(cats, MeanTime_s);
    ylabel('Mean elapsed time (s)');
    title('Computational repeatability');
    save_figure_safely(fig, 'Additional_02_repeatability_statistics.png');

    fprintf('Multiple-seed repeatability block completed. Results saved as CSV and PNG.\n');
end

%% ------------------------------------------------------------------------
% ADDITIONAL BLOCK 3: Voronoi and Delaunay morphology analysis
% -------------------------------------------------------------------------

function plot_voronoi_delaunay_validation_block(result_VGPP, result_RSA, params)

    statsVGPP = compute_voronoi_delaunay_stats(result_VGPP.particleTable, params);
    statsRSA  = compute_voronoi_delaunay_stats(result_RSA.particleTable, params);

    fig = figure('Color','w', 'Name', 'Additional block - Voronoi and Delaunay validation');

    subplot(1,2,1); hold on; box on; grid off;
    if ~isempty(statsVGPP.voronoiAreas)
        histogram(statsVGPP.voronoiAreas, 25, 'Normalization','probability', 'DisplayName','VGPP');
    end
    if ~isempty(statsRSA.voronoiAreas)
        histogram(statsRSA.voronoiAreas, 25, 'Normalization','probability', 'DisplayName','RSA-PBC');
    end
    xlabel('Interior Voronoi cell area');
    ylabel('Probability');
    title('Voronoi cell-area distribution');
    legend('Location','best');

    subplot(1,2,2); hold on; box on; grid off;
    if ~isempty(statsVGPP.delaunayEdges)
        histogram(statsVGPP.delaunayEdges, 25, 'Normalization','probability', 'DisplayName','VGPP');
    end
    if ~isempty(statsRSA.delaunayEdges)
        histogram(statsRSA.delaunayEdges, 25, 'Normalization','probability', 'DisplayName','RSA-PBC');
    end
    xlabel('Delaunay edge length');
    ylabel('Probability');
    title('Delaunay edge-length distribution');
    legend('Location','best');
    save_figure_safely(fig, 'Additional_03_Voronoi_Delaunay_validation.png');

    Method = ["VGPP"; "RSA-PBC"];
    VoronoiAreaMean = [statsVGPP.voronoiAreaMean; statsRSA.voronoiAreaMean];
    VoronoiAreaCV = [statsVGPP.voronoiAreaCV; statsRSA.voronoiAreaCV];
    DelaunayEdgeMean = [statsVGPP.delaunayEdgeMean; statsRSA.delaunayEdgeMean];
    DelaunayEdgeCV = [statsVGPP.delaunayEdgeCV; statsRSA.delaunayEdgeCV];
    NumberOfInteriorVoronoiCells = [numel(statsVGPP.voronoiAreas); numel(statsRSA.voronoiAreas)];
    NumberOfDelaunayEdges = [numel(statsVGPP.delaunayEdges); numel(statsRSA.delaunayEdges)];

    vorDelTable = table(Method, VoronoiAreaMean, VoronoiAreaCV, ...
        DelaunayEdgeMean, DelaunayEdgeCV, NumberOfInteriorVoronoiCells, NumberOfDelaunayEdges);
    disp('Voronoi/Delaunay morphology validation table:');
    disp(vorDelTable);
    writetable(vorDelTable, 'VGPP_Voronoi_Delaunay_validation.csv');

    fprintf('Voronoi/Delaunay validation block completed.\n');
end

function stats = compute_voronoi_delaunay_stats(particleTable, params)

    x = particleTable.x(:);
    y = particleTable.y(:);
    n = numel(x);

    stats.voronoiAreas = [];
    stats.voronoiAreaMean = NaN;
    stats.voronoiAreaCV = NaN;
    stats.delaunayEdges = [];
    stats.delaunayEdgeMean = NaN;
    stats.delaunayEdgeCV = NaN;

    if n < 4
        return;
    end

    pts = [x y];

    % Interior finite Voronoi-cell areas. Boundary cells are excluded because
    % the plot is used as a compact morphology descriptor, not as the final
    % PBC geometry validator.
    try
        [V, C] = voronoin(pts);
        areas = [];
        for i = 1:numel(C)
            cellIdx = C{i};
            if all(cellIdx ~= 1)
                vx = V(cellIdx,1);
                vy = V(cellIdx,2);
                inside = all(vx >= 0 & vx <= params.Lx & vy >= 0 & vy <= params.Ly);
                if inside
                    areas(end+1,1) = polyarea(vx, vy); %#ok<AGROW>
                end
            end
        end
        stats.voronoiAreas = areas;
        if ~isempty(areas)
            stats.voronoiAreaMean = mean(areas, 'omitnan');
            stats.voronoiAreaCV = std(areas, 'omitnan') / max(mean(areas, 'omitnan'), eps);
        end
    catch
        warning('Voronoi analysis failed. Skipping Voronoi area statistics.');
    end

    % Delaunay edge lengths, evaluated with the PBC distance formula so that
    % edge-length statistics remain consistent with the RVE distance metric.
    try
        DT = delaunayTriangulation(pts);
        E = edges(DT);
        if ~isempty(E)
            edgeLengths = zeros(size(E,1),1);
            for k = 1:size(E,1)
                edgeLengths(k) = periodic_center_distance(x(E(k,1)), y(E(k,1)), ...
                    x(E(k,2)), y(E(k,2)), params.Lx, params.Ly);
            end
            stats.delaunayEdges = edgeLengths;
            stats.delaunayEdgeMean = mean(edgeLengths, 'omitnan');
            stats.delaunayEdgeCV = std(edgeLengths, 'omitnan') / max(mean(edgeLengths, 'omitnan'), eps);
        end
    catch
        warning('Delaunay analysis failed. Skipping Delaunay edge statistics.');
    end
end

%% ------------------------------------------------------------------------
% ADDITIONAL BLOCK 4: local volume-fraction uniformity versus window size
% -------------------------------------------------------------------------

function localWindowTable = plot_local_vf_window_sensitivity_block(result_VGPP, result_RSA, params)

    binsList = params.localVFWindowBinsList(:);

    Method = strings(0,1);
    NumberOfBinsPerDirection = [];
    TotalWindows = [];
    LocalVFMean = [];
    LocalVFStd = [];
    LocalVFRSD = [];
    LocalVFMin = [];
    LocalVFMax = [];

    methods = ["VGPP", "RSA-PBC"];
    tables = {result_VGPP.particleTable, result_RSA.particleTable};

    for im = 1:numel(methods)
        tbl = tables{im};
        for ib = 1:numel(binsList)
            nBins = binsList(ib);
            [localMap, rsd] = local_volume_fraction_map_by_centers( ...
                tbl.x, tbl.y, tbl.radius, params.Lx, params.Ly, nBins);
            vals = localMap(:);

            Method(end+1,1) = methods(im);
            NumberOfBinsPerDirection(end+1,1) = nBins;
            TotalWindows(end+1,1) = nBins^2;
            LocalVFMean(end+1,1) = mean(vals, 'omitnan');
            LocalVFStd(end+1,1) = std(vals, 'omitnan');
            LocalVFRSD(end+1,1) = rsd;
            LocalVFMin(end+1,1) = min(vals);
            LocalVFMax(end+1,1) = max(vals);
        end
    end

    localWindowTable = table(Method, NumberOfBinsPerDirection, TotalWindows, ...
        LocalVFMean, LocalVFStd, LocalVFRSD, LocalVFMin, LocalVFMax);
    writetable(localWindowTable, 'VGPP_local_VF_window_sensitivity.csv');

    fig = figure('Color','w', 'Name', 'Additional block - Local VF window sensitivity');
    hold on; box on; grid off;
    for im = 1:numel(methods)
        idx = localWindowTable.Method == methods(im);
        plot(localWindowTable.NumberOfBinsPerDirection(idx), ...
            localWindowTable.LocalVFRSD(idx), '-o', 'LineWidth', 1.8, ...
            'DisplayName', char(methods(im)));
    end
    xlabel('Number of windows per direction');
    ylabel('RSD of local VF');
    title('Local volume-fraction uniformity versus window size');
    legend('Location','best');
    save_figure_safely(fig, 'Additional_04_local_VF_window_sensitivity.png');

    disp('Local VF window-sensitivity table:');
    disp(localWindowTable);
    fprintf('Local VF window-sensitivity block completed.\n');
end

%% ------------------------------------------------------------------------
% Helper for repeatability and sweep blocks
% -------------------------------------------------------------------------

function rsaParams = make_rsa_params_from_vgpp_params(params)
    rsaParams = struct();
    rsaParams.Lx = params.Lx;
    rsaParams.Ly = params.Ly;
    rsaParams.particleMode = params.particleMode;
    rsaParams.rMono = params.rMono;
    rsaParams.rMean = params.rMean;
    rsaParams.rStd = params.rStd;
    rsaParams.rMin = params.rMin;
    rsaParams.rMax = params.rMax;
    rsaParams.minGap = params.minGap;
    if isfield(params,'rsaMaxTrials')
        rsaParams.maxTrials = params.rsaMaxTrials;
    else
        rsaParams.maxTrials = 5000000;
    end
    if isfield(params,'activeTargetVF') && isfinite(params.activeTargetVF)
        rsaParams.targetVF = params.activeTargetVF;
    else
        rsaParams.targetVF = inf;
    end
    if isfield(params,'rsaMaxFailedTrials')
        rsaParams.maxFailedTrials = params.rsaMaxFailedTrials;
    else
        rsaParams.maxFailedTrials = 200000;
    end
end


%% ------------------------------------------------------------------------

function plot_selected_stage_contribution_block(result_VGPP, params)

    % Selected manuscript block: contribution of each VGPP stage.
    % This plot makes the hybrid mechanism transparent by showing how the
    % area fraction increases from initial void-guided adsorption to final
    % overlap-audited packing.

    if isempty(result_VGPP.historyTable)
        warning('Stage contribution block skipped: history table is empty.');
        return;
    end

    H = result_VGPP.historyTable;
    stageNames = string(H.Stage);

    selectedStages = [ ...
        "Initial adsorption"; ...
        "Relaxation-reinsertion"; ...
        "Growth-compression"; ...
        "Hexagonal-seeded dense fallback"];

    displayStages = [ ...
        "Initial void-guided adsorption"; ...
        "Relaxation-reinsertion"; ...
        "Growth-compression"; ...
        "Overlap-audited dense fallback"];

    stageVF = NaN(numel(selectedStages),1);
    stageParticles = NaN(numel(selectedStages),1);

    for i = 1:numel(selectedStages)
        idx = find(stageNames == selectedStages(i));
        if ~isempty(idx)
            stageVF(i) = H.AreaFraction(idx(end));
            stageParticles(i) = H.NumberOfParticles(idx(end));
        end
    end

    valid = ~isnan(stageVF);
    displayStages = displayStages(valid);
    stageVF = stageVF(valid);
    stageParticles = stageParticles(valid);

    if isempty(stageVF)
        warning('Stage contribution block skipped: no recognized stages were found.');
        return;
    end

    fig = figure('Color','w','Name','Selected block - VGPP stage contribution');
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(1:numel(stageVF), stageVF, '-o', 'LineWidth', 2, 'MarkerSize', 7);
    hold on; box on; grid off;
    set(gca,'XTick',1:numel(stageVF),'XTickLabel',displayStages);
    xtickangle(25);
    ylabel('Area/volume fraction');
    title('Contribution of VGPP stages');
    ylim([0, min(1, max(stageVF)*1.15)]);
    for i = 1:numel(stageVF)
        text(i, stageVF(i), sprintf(' %.3f',stageVF(i)), ...
            'VerticalAlignment','bottom','FontWeight','bold');
    end

    nexttile;
    bar(categorical(displayStages), stageParticles);
    box on; grid off;
    ylabel('Number of particles');
    title('Particle count after each stage');
    xtickangle(25);

    stageTable = table(displayStages(:), stageVF(:), stageParticles(:), ...
        'VariableNames', {'Stage','AreaVolumeFraction','NumberOfParticles'});
    disp('Selected stage-contribution table:');
    disp(stageTable);
    writetable(stageTable, 'Selected_stage_contribution_table.csv');

    save_figure_safely(fig, 'Selected_01_VGPP_stage_contribution.png');
end

%% ------------------------------------------------------------------------

function plot_selected_matrix_rich_zone_block(result_VGPP, result_RSA, params)

    % Selected manuscript block: local matrix-rich zone detection.
    % This block highlights cells with unusually low local particle fraction.
    % It is useful because biased RVEs with unrealistic matrix-rich zones are
    % a known concern in RVE-generation algorithms.

    if isempty(result_VGPP.particleTable) || isempty(result_RSA.particleTable)
        warning('Matrix-rich zone block skipped: one particle table is empty.');
        return;
    end

    bins = params.matrixRichBins;
    thresholdFactor = params.matrixRichLowVFThresholdFactor;

    tblV = result_VGPP.particleTable;
    tblR = result_RSA.particleTable;

    [mapV, rsdV] = local_volume_fraction_map_by_centers( ...
        tblV.x, tblV.y, tblV.radius, params.Lx, params.Ly, bins);
    [mapR, rsdR] = local_volume_fraction_map_by_centers( ...
        tblR.x, tblR.y, tblR.radius, params.Lx, params.Ly, bins);

    meanV = mean(mapV(:));
    meanR = mean(mapR(:));

    lowMaskV = mapV < thresholdFactor * meanV;
    lowMaskR = mapR < thresholdFactor * meanR;

    lowAreaFractionV = sum(lowMaskV(:)) / numel(lowMaskV);
    lowAreaFractionR = sum(lowMaskR(:)) / numel(lowMaskR);

    fig = figure('Color','w','Name','Selected block - matrix-rich zone detection');
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    imagesc(mapV); axis equal tight; set(gca,'YDir','normal');
    title(sprintf('VGPP local VF, RSD = %.3f', rsdV));
    xlabel('Subwindow index'); ylabel('Subwindow index'); colorbar;

    nexttile;
    imagesc(mapR); axis equal tight; set(gca,'YDir','normal');
    title(sprintf('RSA-PBC local VF, RSD = %.3f', rsdR));
    xlabel('Subwindow index'); ylabel('Subwindow index'); colorbar;

    nexttile;
    imagesc(lowMaskV); axis equal tight; set(gca,'YDir','normal');
    title(sprintf('VGPP matrix-rich cells = %.1f%%',100*lowAreaFractionV));
    xlabel('Subwindow index'); ylabel('Subwindow index'); colorbar;

    nexttile;
    imagesc(lowMaskR); axis equal tight; set(gca,'YDir','normal');
    title(sprintf('RSA matrix-rich cells = %.1f%%',100*lowAreaFractionR));
    xlabel('Subwindow index'); ylabel('Subwindow index'); colorbar;

    Method = ["VGPP"; "RSA-PBC"];
    FinalVF = [result_VGPP.maxAreaFraction; result_RSA.maxAreaFraction];
    LocalVFRSD = [rsdV; rsdR];
    MatrixRichCellFraction = [lowAreaFractionV; lowAreaFractionR];
    MeanLocalVF = [meanV; meanR];

    matrixRichTable = table(Method, FinalVF, MeanLocalVF, LocalVFRSD, MatrixRichCellFraction);
    disp('Selected matrix-rich-zone table:');
    disp(matrixRichTable);
    writetable(matrixRichTable, 'Selected_matrix_rich_zone_table.csv');

    save_figure_safely(fig, 'Selected_02_matrix_rich_zone_detection.png');
end

%% ------------------------------------------------------------------------

function plot_selected_clearance_safety_block(result_VGPP, result_RSA, params)

    % Selected manuscript block: periodic non-overlap and near-contact safety.
    % This block is useful because final high-VF configurations must be accepted
    % only if all particle pairs and periodic images pass the overlap audit.

    if isempty(result_VGPP.particleTable) || isempty(result_RSA.particleTable)
        warning('Clearance safety block skipped: one particle table is empty.');
        return;
    end

    tblV = result_VGPP.particleTable;
    tblR = result_RSA.particleTable;

    clearanceV = nearest_surface_clearance_PBC(tblV.x, tblV.y, tblV.radius, params);
    clearanceR = nearest_surface_clearance_PBC(tblR.x, tblR.y, tblR.radius, params);

    pairClearanceV = all_pair_surface_clearance_PBC(tblV.x, tblV.y, tblV.radius, params);
    pairClearanceR = all_pair_surface_clearance_PBC(tblR.x, tblR.y, tblR.radius, params);

    [hasOverlapV, maxOverlapV, minClearanceV, nOverlapV] = ...
        final_overlap_audit_PBC(tblV.x, tblV.y, tblV.radius, params);
    [hasOverlapR, maxOverlapR, minClearanceR, nOverlapR] = ...
        final_overlap_audit_PBC(tblR.x, tblR.y, tblR.radius, params);

    nearThreshold = params.nearContactThresholdFactor * mean(tblV.radius);
    nearContactV = mean(pairClearanceV >= -1e-10 & pairClearanceV <= nearThreshold);
    nearContactR = mean(pairClearanceR >= -1e-10 & pairClearanceR <= nearThreshold);

    fig = figure('Color','w','Name','Selected block - clearance and overlap safety');
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    hold on; box on; grid off;
    histogram(clearanceV, 35, 'Normalization','pdf', 'DisplayName','VGPP');
    histogram(clearanceR, 35, 'Normalization','pdf', 'DisplayName','RSA-PBC');
    xline(0,'--','Overlap limit','HandleVisibility','off');
    xlabel('Nearest surface-to-surface clearance');
    ylabel('Probability density');
    title('Nearest-clearance distribution');
    legend('Location','best');

    nexttile;
    bar(categorical(MethodNamesForTwo()), [nOverlapV nOverlapR; nearContactV nearContactR]');
    ylabel('Count / fraction');
    title('Overlap and near-contact audit');
    legend({'Overlapping pairs','Near-contact pair fraction'},'Location','best');
    grid off; box on;

    Method = ["VGPP"; "RSA-PBC"];
    HasOverlap = [hasOverlapV; hasOverlapR];
    OverlappingPairs = [nOverlapV; nOverlapR];
    MaximumOverlap = [maxOverlapV; maxOverlapR];
    MinimumClearance = [minClearanceV; minClearanceR];
    NearContactPairFraction = [nearContactV; nearContactR];

    safetyTable = table(Method, HasOverlap, OverlappingPairs, MaximumOverlap, ...
        MinimumClearance, NearContactPairFraction);
    disp('Selected clearance-safety table:');
    disp(safetyTable);
    writetable(safetyTable, 'Selected_clearance_safety_table.csv');

    save_figure_safely(fig, 'Selected_03_clearance_overlap_safety.png');
end

%% ------------------------------------------------------------------------

function names = MethodNamesForTwo()
    names = ["VGPP", "RSA-PBC"];
end

%% ------------------------------------------------------------------------

function pairClearance = all_pair_surface_clearance_PBC(cxList, cyList, rList, params)

    nParticles = numel(rList);
    if nParticles < 2
        pairClearance = [];
        return;
    end

    nPairs = nParticles*(nParticles-1)/2;
    pairClearance = zeros(nPairs,1);
    counter = 0;

    for i = 1:nParticles-1
        xOthers = cxList(i+1:end);
        yOthers = cyList(i+1:end);
        rOthers = rList(i+1:end);

        d = periodic_center_distance(cxList(i), cyList(i), xOthers, yOthers, ...
            params.Lx, params.Ly);
        requiredDistance = rList(i) + rOthers + params.minGap;
        c = d - requiredDistance;

        nLocal = numel(c);
        pairClearance(counter+1:counter+nLocal) = c(:);
        counter = counter + nLocal;
    end
end

%% ------------------------------------------------------------------------

function plot_selected_target_reachability_note()
    % Placeholder intentionally unused. Target-VF reachability uses the existing
    % run_target_vf_sweep_block(params) function when the corresponding flag is
    % enabled. It is disabled by default because it requires many full VGPP/RSA
    % reruns and is intended for final manuscript data generation.
end

%% ------------------------------------------------------------------------

function save_figure_safely(figHandle, fileName)
    try
        [parentFolder,~,~] = fileparts(fileName);
        if ~isempty(parentFolder) && ~exist(parentFolder,'dir')
            mkdir(parentFolder);
        end
        apply_unified_publication_style(figHandle);
        exportgraphics(figHandle, fileName, 'Resolution', 600);
    catch
        try
            saveas(figHandle, fileName);
        catch
            warning('Could not save figure: %s', fileName);
        end
    end
end


%% ------------------------------------------------------------------------
function plot_publication_subplots_Fig3_to_Fig7(result_VGPP, result_RSA, params, comparisonTable)

    % Consolidated publication-ready figures for Results and Discussion.
    % These figures combine the most informative outputs into related subplots.

    finalTable = result_VGPP.particleTable;
    rsaTable   = result_RSA.particleTable;

    if ~exist('VGPP_publication_figures','dir')
        mkdir('VGPP_publication_figures');
    end

    %% ========================== Fig. 3 ================================
    % Corrected layout for Fig. 3:
    %   - enlarged 2 x 2 canvas with looser spacing
    %   - subplot labels are drawn as in-axes text annotations so they remain visible
    %   - packing maps occupy the first row, and quantitative comparisons occupy the second row

    figure('Color','w','Name','Fig. 3 VGPP versus RSA packing and VF', ...
        'Units','pixels','Position',[60 60 1550 1050]);
    t3 = tiledlayout(2,2,'TileSpacing','loose','Padding','compact');

    ax31 = nexttile(t3,1);
    draw_particles_publication(finalTable.x, finalTable.y, finalTable.radius, params);
    set(ax31,'FontName','Times New Roman','FontSize',14,'FontWeight','bold','LineWidth',2, ...
        'XTick',[0 params.Lx/2 params.Lx],'YTick',[0 params.Ly/2 params.Ly]);
    xlabel(ax31,'x','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    ylabel(ax31,'y','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    text(ax31,0.04,0.96,'(a)','Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','top','FontName','Times New Roman','FontSize',28,'FontWeight','bold');

    ax32 = nexttile(t3,2);
    draw_particles_publication(rsaTable.x, rsaTable.y, rsaTable.radius, params);
    set(ax32,'FontName','Times New Roman','FontSize',14,'FontWeight','bold','LineWidth',2, ...
        'XTick',[0 params.Lx/2 params.Lx],'YTick',[0 params.Ly/2 params.Ly]);
    xlabel(ax32,'x','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    ylabel(ax32,'y','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    text(ax32,0.04,0.96,'(b)','Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','top','FontName','Times New Roman','FontSize',28,'FontWeight','bold');

    ax33 = nexttile(t3,3);
    hold(ax33,'on'); box(ax33,'on'); grid(ax33,'off');
    if ~isempty(result_VGPP.historyTable)
        xVGPP = 1:height(result_VGPP.historyTable);
        plot(ax33, xVGPP, result_VGPP.historyTable.AreaFraction, 'LineWidth', 2, 'DisplayName','VGPP');
    end
    if ~isempty(result_RSA.historyTable)
        xRSA = linspace(1, max(2,height(result_VGPP.historyTable)), height(result_RSA.historyTable));
        plot(ax33, xRSA, result_RSA.historyTable.AreaFraction, '--', 'LineWidth', 2, 'DisplayName','RSA-PBC');
    end
    xlabel(ax33,'Normalized algorithmic step');
    ylabel(ax33,'Volume fraction');
    ylim(ax33,[0 1.0]);
    legend(ax33,'Location','southeast');
    text(ax33,0.04,0.96,'(c)','Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','top','FontName','Times New Roman','FontSize',28,'FontWeight','bold');

    ax34 = nexttile(t3,4);
    bar(ax34, categorical({'RSA-PBC','VGPP'}), [result_RSA.maxAreaFraction, result_VGPP.maxAreaFraction], 'LineWidth', 2);
    ylabel(ax34,'Final volume fraction');
    ylim(ax34,[0, 1.0]);
    grid(ax34,'off'); box(ax34,'on');
    text(ax34,0.04,0.96,'(d)','Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','top','FontName','Times New Roman','FontSize',28,'FontWeight','bold');

    apply_unified_publication_style(gcf);
    exportgraphics(gcf, fullfile('VGPP_publication_figures','Fig3_VGPP_RSA_packing_and_VF.png'), 'Resolution', 600);

    %% ========================== Fig. 4 ================================
    % Corrected layout for Fig. 4:
    %   - enlarged 2 x 2 canvas with looser spacing
    %   - subplot labels are drawn as in-axes text annotations so they remain visible
    %   - color maps occupy the first row, and quantitative distributions/metrics occupy the second row

    nBins = 8;
    [localVGPP, rsdVGPP] = local_vf_map_by_centers_publication(finalTable.x, finalTable.y, finalTable.radius, params, nBins);
    [localRSA,  rsdRSA]  = local_vf_map_by_centers_publication(rsaTable.x, rsaTable.y, rsaTable.radius, params, nBins);

    meanVGPP = mean(localVGPP(:), 'omitnan');
    meanRSA  = mean(localRSA(:),  'omitnan');

    thresholdVGPP = params.matrixRichLowVFThresholdFactor * meanVGPP;
    thresholdRSA  = params.matrixRichLowVFThresholdFactor * meanRSA;

    lowVGPP = 100 * mean(localVGPP(:) < thresholdVGPP, 'omitnan');
    lowRSA  = 100 * mean(localRSA(:)  < thresholdRSA,  'omitnan');

    rsdVGPP_percent = 100 * rsdVGPP;
    rsdRSA_percent  = 100 * rsdRSA;

    maxColor = max([localVGPP(:); localRSA(:); 1.0]);

    figure('Color','w','Name','Fig. 4 Local VF uniformity', ...
        'Units','pixels','Position',[80 80 1550 1050]);
    t4 = tiledlayout(2,2,'TileSpacing','loose','Padding','compact');

    % ---- (a) VGPP local volume-fraction map
    ax41 = nexttile(t4,1);
    imagesc(ax41, [0 params.Lx], [0 params.Ly], localVGPP);
    axis(ax41,'equal','tight');
    set(ax41,'YDir','normal','FontName','Times New Roman','FontSize',14,'FontWeight','bold','LineWidth',2);
    grid(ax41,'off'); box(ax41,'on');
    caxis(ax41,[0 maxColor]);
    xlabel(ax41,'x','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    ylabel(ax41,'y','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    text(ax41,0.04,0.96,'(a)','Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','top','FontName','Times New Roman','FontSize',26,'FontWeight','bold');
    cb1 = colorbar(ax41);
    cb1.Label.String = 'Local VF';
    cb1.Label.FontName = 'Times New Roman';
    cb1.Label.FontSize = 14;
    cb1.Label.FontWeight = 'bold';
    cb1.FontName = 'Times New Roman';
    cb1.FontSize = 12;
    cb1.FontWeight = 'bold';
    cb1.LineWidth = 1.5;

    % ---- (b) RSA-PBC local volume-fraction map
    ax42 = nexttile(t4,2);
    imagesc(ax42, [0 params.Lx], [0 params.Ly], localRSA);
    axis(ax42,'equal','tight');
    set(ax42,'YDir','normal','FontName','Times New Roman','FontSize',14,'FontWeight','bold','LineWidth',2);
    grid(ax42,'off'); box(ax42,'on');
    caxis(ax42,[0 maxColor]);
    xlabel(ax42,'x','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    ylabel(ax42,'y','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    text(ax42,0.04,0.96,'(b)','Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','top','FontName','Times New Roman','FontSize',26,'FontWeight','bold');
    cb2 = colorbar(ax42);
    cb2.Label.String = 'Local VF';
    cb2.Label.FontName = 'Times New Roman';
    cb2.Label.FontSize = 14;
    cb2.Label.FontWeight = 'bold';
    cb2.FontName = 'Times New Roman';
    cb2.FontSize = 12;
    cb2.FontWeight = 'bold';
    cb2.LineWidth = 1.5;

    % ---- (c) Cell-level local VF distribution
    ax43 = nexttile(t4,3);
    hold(ax43,'on'); box(ax43,'on'); grid(ax43,'off');
    histogram(ax43, localVGPP(:), 'Normalization','probability', 'BinWidth',0.05, ...
        'DisplayName','VGPP');
    histogram(ax43, localRSA(:), 'Normalization','probability', 'BinWidth',0.05, ...
        'DisplayName','RSA-PBC');
    xline(ax43, meanVGPP, '-',  'HandleVisibility','off', 'LineWidth',2);
    xline(ax43, meanRSA,  '--', 'HandleVisibility','off', 'LineWidth',2);
    text(ax43, min(meanVGPP + 0.006, 0.88), 0.76, 'VGPP mean', 'FontName','Times New Roman', ...
        'FontSize',12, 'FontWeight','bold', 'HorizontalAlignment','left');
    text(ax43, meanRSA + 0.008, 0.76, 'RSA mean', 'FontName','Times New Roman', ...
        'FontSize',12, 'FontWeight','bold', 'HorizontalAlignment','left');
    xlabel(ax43,'Local volume fraction');
    ylabel(ax43,'Cell probability');
    xlim(ax43,[max(0,min([localVGPP(:);localRSA(:)])-0.03) min(1,max([localVGPP(:);localRSA(:)])+0.03)]);
    ylim(ax43,[0 0.8]);
    legend(ax43,'Location','southeast','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
    text(ax43,0.04,0.96,'(c)','Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','top','FontName','Times New Roman','FontSize',26,'FontWeight','bold');
    set(ax43,'FontName','Times New Roman','FontSize',14,'FontWeight','bold','LineWidth',2);

    % ---- (d) Matrix-rich and uniformity metrics
    ax44 = nexttile(t4,4);
    metricNames = categorical({'Low-VF cells','Local VF RSD'});
    metricNames = reordercats(metricNames, {'Low-VF cells','Local VF RSD'});
    metricData = [lowVGPP lowRSA; rsdVGPP_percent rsdRSA_percent];
    b = bar(ax44, metricNames, metricData, 'grouped', 'LineWidth', 2);
    ylabel(ax44,'Metric value (%)');
    ylim(ax44,[0 max(14,1.15*max(metricData(:)))]);
    legend(ax44,{'VGPP','RSA-PBC'},'Location','southwest','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
    grid(ax44,'off'); box(ax44,'on');
    text(ax44,0.04,0.96,'(d)','Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','top','FontName','Times New Roman','FontSize',26,'FontWeight','bold');
    set(ax44,'FontName','Times New Roman','FontSize',14,'FontWeight','bold','LineWidth',2);

    for ib = 1:numel(b)
        xtips = b(ib).XEndPoints;
        ytips = b(ib).YEndPoints;
        labels = string(round(b(ib).YData,2));
        text(ax44, xtips, ytips, labels, 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', 'FontName','Times New Roman', 'FontSize',13, 'FontWeight','bold');
    end

    apply_unified_publication_style(gcf);
    exportgraphics(gcf, fullfile('VGPP_publication_figures','Fig4_local_VF_uniformity.png'), 'Resolution', 600);

    %% ========================== Fig. 5 ================================
    % Corrected layout for Fig. 5:
    %   - reduced figure height to obtain a more compact, square-like layout
    %   - explicit legend handles are used so the vertical reference line appears as "Overlap limit"
    %   - overlap-limit line is red, dashed, and thicker for better distinction
    %   - subplot labels are placed inside the axes where they cannot be hidden by the legend

    clearVGPP = nearest_clearance_publication(finalTable.x, finalTable.y, finalTable.radius, params);
    clearRSA  = nearest_clearance_publication(rsaTable.x, rsaTable.y, rsaTable.radius, params);
    nearThreshold = params.nearContactThresholdFactor * params.rMono;

    % Numerical overlap check before plotting. Negative clearance means overlap.
    % Very small negative values can occur only from numerical round-off, so they
    % are set to zero for plotting, but real negative values trigger an error.
    clearanceTol = 1e-10;
    minClearVGPP = min(clearVGPP,[],'omitnan');
    minClearRSA  = min(clearRSA,[],'omitnan');

    if minClearVGPP < -clearanceTol || minClearRSA < -clearanceTol
        error(['Negative surface-to-surface clearance was detected. ', ...
               'This indicates particle overlap. Check the final overlap audit ', ...
               'before using Fig. 5. Min VGPP = %.4e, Min RSA-PBC = %.4e'], ...
               minClearVGPP, minClearRSA);
    end

    clearVGPP_plot = max(clearVGPP,0);
    clearRSA_plot  = max(clearRSA,0);
    allClearPlot   = [clearVGPP_plot(:); clearRSA_plot(:)];
    xMaxClear      = max(allClearPlot,[],'omitnan');
    if isempty(xMaxClear) || ~isfinite(xMaxClear) || xMaxClear <= 0
        xMaxClear = 1;
    end

    % Use non-negative bin edges so Fig. 5a cannot visually suggest overlap.
    nClearBins = 35;
    clearEdges = linspace(0,1.05*xMaxClear,nClearBins+1);

    figure('Color','w','Name','Fig. 5 Clearance and overlap safety', ...
        'Units','pixels','Position',[80 80 1350 620]);
    t5 = tiledlayout(1,2,'TileSpacing','loose','Padding','compact');

    % ---- (a) Nearest-clearance distribution
    ax51 = nexttile(t5,1);
    hold(ax51,'on'); box(ax51,'on'); grid(ax51,'off');

    hVGPP = histogram(ax51, clearVGPP_plot, clearEdges, ...
        'Normalization','pdf', ...
        'DisplayName','VGPP', ...
        'LineWidth',2);

    hRSA = histogram(ax51, clearRSA_plot, clearEdges, ...
        'Normalization','pdf', ...
        'DisplayName','RSA-PBC', ...
        'LineWidth',2);

    hLimit = xline(ax51, 0, '--', ...
        'DisplayName','Overlap limit', ...
        'Color',[1 0 0], ...
        'LineWidth',3);

    xlabel(ax51,'Nearest surface-to-surface clearance', ...
        'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    ylabel(ax51,'Probability density', ...
        'FontName','Times New Roman','FontSize',14,'FontWeight','bold');

    xlim(ax51,[0 1.05*xMaxClear]);
    yMaxClear = max([hVGPP.Values(:); hRSA.Values(:)],[],'omitnan');
    if isempty(yMaxClear) || ~isfinite(yMaxClear) || yMaxClear <= 0
        yMaxClear = 1;
    end
    ylim(ax51,[0 1.18*yMaxClear]);

    legend(ax51,[hVGPP hRSA hLimit], {'VGPP','RSA-PBC','Overlap limit'}, ...
        'Location','northeast', ...
        'FontName','Times New Roman', ...
        'FontSize',13, ...
        'FontWeight','bold', ...
        'Box','on');

    text(ax51,0.06,0.93,'(a)', ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontName','Times New Roman', ...
        'FontSize',22, ...
        'FontWeight','bold');

    % ---- (b) Near-contact particle fraction
    ax52 = nexttile(t5,2);
    nearVGPP = 100 * mean(clearVGPP_plot <= nearThreshold);
    nearRSA  = 100 * mean(clearRSA_plot  <= nearThreshold);

    bar(ax52, categorical({'RSA-PBC','VGPP'}), [nearRSA, nearVGPP], 'LineWidth', 2);
    ylabel(ax52,'Near-contact particles (%)', ...
        'FontName','Times New Roman','FontSize',14,'FontWeight','bold');

    ylim(ax52,[0 110]);
    grid(ax52,'off'); box(ax52,'on');

    text(ax52,0.06,0.93,'(b)', ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontName','Times New Roman', ...
        'FontSize',22, ...
        'FontWeight','bold');

    apply_unified_publication_style(gcf);
    exportgraphics(gcf, fullfile('VGPP_publication_figures','Fig5_clearance_overlap_safety.png'), 'Resolution', 600);

    %% ========================== Fig. 6 ================================
    % Corrected and polished layout for Fig. 6.
    figure('Color','w','Name','Fig. 6 Contribution of VGPP stages', ...
        'Units','pixels','Position',[60 60 1600 700]);
    t6 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    stages = ["Initial adsorption", "Relaxation-reinsertion", ...
        "Growth-compression", "Hexagonal-seeded dense fallback"];
    stageLabels = {'Initial', 'Relaxation', 'Growth comp.', 'Audited fallback'};
    xStage = 1:numel(stages);
    stageVF = nan(numel(stages),1);

    for k = 1:numel(stages)
        idx = result_VGPP.historyTable.Stage == stages(k);
        if any(idx)
            stageVF(k) = max(result_VGPP.historyTable.AreaFraction(idx));
        end
    end
    if isnan(stageVF(1))
        stageVF(1) = min(result_VGPP.historyTable.AreaFraction);
    end
    for k = 2:numel(stageVF)
        if isnan(stageVF(k))
            stageVF(k) = stageVF(k-1);
        end
    end
    increments = diff([0; stageVF(:)]);

    ax61 = nexttile(t6,1);
    plot(ax61, xStage, stageVF, '-o', 'LineWidth', 2, 'MarkerSize', 9, 'MarkerFaceColor','w');
    xlim(ax61,[0.7, numel(stages)+0.3]);
    xticks(ax61,xStage);
    xticklabels(ax61,stageLabels);
    xtickangle(ax61,0);
    xlabel(ax61,'Stage');
    ylabel(ax61,'Volume fraction');
    ylim(ax61,[0, max(0.95, 1.08*max(stageVF))]);
    grid(ax61,'off'); box(ax61,'on');
    set(ax61,'FontName','Times New Roman','FontSize',14,'FontWeight','bold','LineWidth',2);
    text(ax61, 0.04, 0.94, '(a)', 'Units','normalized', 'HorizontalAlignment','left', ...
        'VerticalAlignment','top', 'FontName','Times New Roman', 'FontSize',20, 'FontWeight','bold');
    for k = 1:numel(xStage)
        text(ax61, xStage(k), stageVF(k) + 0.02, sprintf('%.3f', stageVF(k)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontName','Times New Roman', 'FontSize',12, 'FontWeight','bold');
    end

    ax62 = nexttile(t6,2);
    bar(ax62, xStage, increments, 'LineWidth', 2);
    xlim(ax62,[0.5, numel(stages)+0.5]);
    xticks(ax62,xStage);
    xticklabels(ax62,stageLabels);
    xtickangle(ax62,0);
    xlabel(ax62,'Stage');
    ylabel(ax62,'VF increment');
    ylim(ax62,[0, max(0.60, 1.18*max(increments))]);
    grid(ax62,'off'); box(ax62,'on');
    set(ax62,'FontName','Times New Roman','FontSize',14,'FontWeight','bold','LineWidth',2);
    text(ax62, 0.04, 0.94, '(b)', 'Units','normalized', 'HorizontalAlignment','left', ...
        'VerticalAlignment','top', 'FontName','Times New Roman', 'FontSize',20, 'FontWeight','bold');
    for k = 1:numel(xStage)
        text(ax62, xStage(k), increments(k) + 0.015, sprintf('%.3f', increments(k)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontName','Times New Roman', 'FontSize',12, 'FontWeight','bold');
    end

    apply_unified_publication_style(gcf);
    exportgraphics(gcf, fullfile('VGPP_publication_figures','Fig6_stage_contribution.png'), 'Resolution', 600);

    %% ========================== Fig. 7 ================================
    % Corrected and re-centered layout for Fig. 7.
    figure('Color','w','Name','Fig. 7 Reconstructed VGPP packing evolution', ...
        'Units','pixels','Position',[40 40 1850 1100]);

    targetVFList = [0.50 0.60 0.70 0.80 result_VGPP.maxAreaFraction];
    revealOrder = spatial_reveal_order_publication(finalTable.x, finalTable.y, params, 10);

    % Use a 2 x 6 tiled layout so the top row contains 3 centered wide panels
    % and the bottom row contains 2 centered wide panels.
    t7 = tiledlayout(2,6,'TileSpacing','compact','Padding','compact');
    tileStarts = [1 3 5 8 10];

    for k = 1:numel(targetVFList)
        vfTarget = targetVFList(k);
        if k == numel(targetVFList)
            idx = 1:height(finalTable);
        else
            nKeep = min(height(finalTable), max(1, round(vfTarget * params.Lx * params.Ly / ...
                (pi * mean(finalTable.radius)^2))));
            idx = revealOrder(1:nKeep);
        end

        reconstructedVF = sum(pi*finalTable.radius(idx).^2)/(params.Lx*params.Ly);

        ax = nexttile(t7, tileStarts(k), [1 2]);
        draw_particles_publication(finalTable.x(idx), finalTable.y(idx), finalTable.radius(idx), params);
        set(ax,'FontName','Times New Roman','FontSize',14,'FontWeight','bold','LineWidth',2, ...
            'XTick',[0 params.Lx/2 params.Lx], 'YTick',[0 params.Ly/2 params.Ly]);
        xlabel(ax,'x','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
        ylabel(ax,'y','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
        grid(ax,'off'); box(ax,'on');
        text(ax, 0.04, 1.03, sprintf('(%c)  VF = %.3f', char(96+k), reconstructedVF), ...
            'Units','normalized', 'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
            'FontName','Times New Roman', 'FontSize',16, 'FontWeight','bold', 'Interpreter','none', ...
            'Clipping','off');
    end

    apply_unified_publication_style(gcf);
    exportgraphics(gcf, fullfile('VGPP_publication_figures','Fig7_packing_evolution.png'), 'Resolution', 600);

    %% ========================== Fig. 8 ================================
    % Dominance assessment of VGPP against RSA-PBC.
    % This figure summarizes the main quantitative advantages of VGPP using
    % the already calculated metrics from Figs. 3-5.

    vfVGPP_percent  = 100 * result_VGPP.maxAreaFraction;
    vfRSA_percent   = 100 * result_RSA.maxAreaFraction;
    rsdVGPP_percent = 100 * rsdVGPP;
    rsdRSA_percent  = 100 * rsdRSA;

    vfGain_percent = 100 * (result_VGPP.maxAreaFraction - result_RSA.maxAreaFraction) / ...
        max(result_RSA.maxAreaFraction, eps);

    rsdReduction_percent = 100 * (rsdRSA - rsdVGPP) / max(rsdRSA, eps);

    if lowRSA > eps
        matrixRichReduction_percent = 100 * (lowRSA - lowVGPP) / lowRSA;
    else
        matrixRichReduction_percent = 0;
    end

    figure('Color','w','Name','Fig. 8 VGPP dominance assessment', ...
        'Units','pixels','Position',[80 80 1650 720]);
    t8 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % ---- (a) Direct comparison of key RVE-quality metrics
    ax81 = nexttile(t8,1);
    metricNamesDirect = categorical({'Final VF','Local VF RSD','Low-VF cells'});
    metricNamesDirect = reordercats(metricNamesDirect, {'Final VF','Local VF RSD','Low-VF cells'});
    directData = [vfVGPP_percent vfRSA_percent; rsdVGPP_percent rsdRSA_percent; lowVGPP lowRSA];
    b81 = bar(ax81, metricNamesDirect, directData, 'grouped', 'LineWidth', 2);
    ylabel(ax81,'Value (%)','FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    ylim(ax81,[0, max(100,1.18*max(directData(:)))]);
    legend(ax81,{'VGPP','RSA-PBC'},'Location','northoutside','Orientation','horizontal', ...
        'FontName','Times New Roman','FontSize',12,'FontWeight','bold','Box','on');
    grid(ax81,'off'); box(ax81,'on');
    set(ax81,'FontName','Times New Roman','FontSize',14,'FontWeight','bold','LineWidth',2);
    text(ax81,0.04,0.93,'(a)','Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','top','FontName','Times New Roman','FontSize',26,'FontWeight','bold');

    for ib = 1:numel(b81)
        xtips = b81(ib).XEndPoints;
        ytips = b81(ib).YEndPoints;
        labels = string(round(b81(ib).YData,2));
        text(ax81, xtips, ytips + 1.5, labels, 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', 'FontName','Times New Roman', 'FontSize',12, 'FontWeight','bold');
    end

    % ---- (b) VGPP improvement with respect to RSA-PBC
    ax82 = nexttile(t8,2);
    improvementNames = categorical({'Packing capacity','Uniformity','Matrix-rich suppression'});
    improvementNames = reordercats(improvementNames, {'Packing capacity','Uniformity','Matrix-rich suppression'});
    improvementData = [vfGain_percent, rsdReduction_percent, matrixRichReduction_percent];
    b82 = bar(ax82, improvementNames, improvementData, 'LineWidth', 2);
    ylabel(ax82,'VGPP improvement over RSA-PBC (%)', ...
        'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    ylim(ax82,[0, max(110,1.18*max(improvementData))]);
    grid(ax82,'off'); box(ax82,'on');
    set(ax82,'FontName','Times New Roman','FontSize',14,'FontWeight','bold','LineWidth',2);
    text(ax82,0.04,0.93,'(b)','Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','top','FontName','Times New Roman','FontSize',26,'FontWeight','bold');

    xtips = b82.XEndPoints;
    ytips = b82.YEndPoints;
    labels = strcat(string(round(b82.YData,1)), '%');
    text(ax82, xtips, ytips + 2.0, labels, 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontName','Times New Roman', 'FontSize',13, 'FontWeight','bold');

    apply_unified_publication_style(gcf);
    exportgraphics(gcf, fullfile('VGPP_publication_figures','Fig8_VGPP_dominance_assessment.png'), 'Resolution', 600);

    %% Summary metrics
    packingImprovement = 100 * (result_VGPP.maxAreaFraction - result_RSA.maxAreaFraction) / result_RSA.maxAreaFraction;
    summaryTable = table( ...
        ["VGPP final VF"; "RSA final VF"; "Packing improvement (%)"; "VGPP particle number"; "RSA particle number"; "VGPP local VF RSD"; "RSA local VF RSD"; "VGPP min clearance"; "RSA min clearance"], ...
        [result_VGPP.maxAreaFraction; result_RSA.maxAreaFraction; packingImprovement; height(finalTable); height(rsaTable); rsdVGPP; rsdRSA; min(clearVGPP); min(clearRSA)], ...
        'VariableNames', {'Metric','Value'});
    writetable(summaryTable, fullfile('VGPP_publication_figures','Fig3_to_Fig7_summary_metrics.csv'));
    disp('Publication subplot figures Fig. 3 to Fig. 7 saved in VGPP_publication_figures.');
    disp(summaryTable);
end

%% ------------------------------------------------------------------------
function draw_particles_publication(cx, cy, rr, params)
    hold on; box on;
    axis equal;
    xlim([0 params.Lx]);
    ylim([0 params.Ly]);

    % Draw particles together with their periodic images. This makes PBC
    % visually explicit in Fig. 3 and Fig. 7: a particle crossing one edge
    % reappears from the opposite edge.
    shiftX = [-params.Lx, 0, params.Lx];
    shiftY = [-params.Ly, 0, params.Ly];

    for i = 1:numel(rr)
        for sx = shiftX
            for sy = shiftY
                xShifted = cx(i) + sx;
                yShifted = cy(i) + sy;

                if (xShifted + rr(i) >= 0) && (xShifted - rr(i) <= params.Lx) && ...
                   (yShifted + rr(i) >= 0) && (yShifted - rr(i) <= params.Ly)

                    rectangle('Position',[xShifted-rr(i), yShifted-rr(i), 2*rr(i), 2*rr(i)], ...
                        'Curvature',[1 1], ...
                        'FaceColor',[0.78 0.78 0.78], ...
                        'EdgeColor','k', ...
                        'LineWidth',0.35, ...
                        'Clipping','on');
                end
            end
        end
    end

    % Redraw the primary RVE boundary on top of periodic image particles.
    rectangle('Position',[0 0 params.Lx params.Ly], ...
        'EdgeColor','k','LineWidth',1.8,'LineStyle','-','FaceColor','none');

    xlabel('x');
    ylabel('y');
    set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','bold','LineWidth',1.5,'Box','on');
end

%% ------------------------------------------------------------------------
function clearance = nearest_clearance_publication(cx, cy, rr, params)
    n = numel(rr);
    clearance = nan(n,1);
    for i = 1:n
        idx = true(n,1); idx(i) = false;
        d = periodic_center_distance(cx(i), cy(i), cx(idx), cy(idx), params.Lx, params.Ly);
        clearance(i) = min(d - rr(i) - rr(idx) - params.minGap);
    end
end

%% ------------------------------------------------------------------------
function [localMap, rsdValue] = local_vf_map_by_centers_publication(cx, cy, rr, params, nBins)
    localMap = zeros(nBins,nBins);
    cellArea = (params.Lx/nBins) * (params.Ly/nBins);
    for i = 1:numel(rr)
        ix = floor(cx(i)/params.Lx*nBins) + 1;
        iy = floor(cy(i)/params.Ly*nBins) + 1;
        ix = min(max(ix,1),nBins);
        iy = min(max(iy,1),nBins);
        localMap(iy,ix) = localMap(iy,ix) + pi*rr(i)^2/cellArea;
    end
    rsdValue = std(localMap(:))/mean(localMap(:));
end

%% ------------------------------------------------------------------------
function order = spatial_reveal_order_publication(cx, cy, params, nBins)
    rng(params.evolutionRandomSeed);
    n = numel(cx);
    bins = cell(nBins,nBins);
    for i = 1:n
        ix = min(max(floor(cx(i)/params.Lx*nBins)+1,1),nBins);
        iy = min(max(floor(cy(i)/params.Ly*nBins)+1,1),nBins);
        bins{iy,ix}(end+1) = i;
    end
    for iy = 1:nBins
        for ix = 1:nBins
            bins{iy,ix} = bins{iy,ix}(randperm(numel(bins{iy,ix})));
        end
    end
    order = [];
    stillHasParticles = true;
    while stillHasParticles
        stillHasParticles = false;
        binOrder = randperm(nBins*nBins);
        for b = binOrder
            [iy,ix] = ind2sub([nBins,nBins], b);
            if ~isempty(bins{iy,ix})
                order(end+1,1) = bins{iy,ix}(1); %#ok<AGROW>
                bins{iy,ix}(1) = [];
                stillHasParticles = true;
            end
        end
    end
end

%% ------------------------------------------------------------------------
function apply_unified_publication_style(figHandle)
    % Unified post-processing style for all publication figures.
    % Required format:
    %   - Times New Roman
    %   - grid off, box on
    %   - line thickness = 2
    %   - bold x/y labels and bold tick labels, font size 14
    %   - alphabet-only subplot labels positioned above each axes

    if nargin < 1 || isempty(figHandle)
        figHandle = gcf;
    end

    set(figHandle,'Color','w');

    axList = findall(figHandle,'Type','axes');
    for ia = 1:numel(axList)
        ax = axList(ia);
        try
            set(ax, 'FontName','Times New Roman', ...
                'FontSize',14, ...
                'FontWeight','bold', ...
                'LineWidth',2, ...
                'Box','on');
            grid(ax,'off');

            if ~isempty(ax.XLabel)
                set(ax.XLabel, 'FontName','Times New Roman', ...
                    'FontSize',14, 'FontWeight','bold');
            end
            if ~isempty(ax.YLabel)
                set(ax.YLabel, 'FontName','Times New Roman', ...
                    'FontSize',14, 'FontWeight','bold');
            end
            if ~isempty(ax.ZLabel)
                set(ax.ZLabel, 'FontName','Times New Roman', ...
                    'FontSize',14, 'FontWeight','bold');
            end
            if ~isempty(ax.Title) && strlength(string(ax.Title.String)) > 0
                set(ax.Title, 'FontName','Times New Roman', ...
                    'FontSize',18, 'FontWeight','bold', 'Interpreter','none');
                ax.Title.Units = 'normalized';
                ax.Title.Position(2) = 1.08;
            end
        catch
            % Continue styling the remaining axes even if one object type fails.
        end
    end

    % Make all plotted curves and markers visually clear.
    lineList = findall(figHandle,'Type','Line');
    for il = 1:numel(lineList)
        try
            set(lineList(il),'LineWidth',2);
        catch
        end
    end

    % Make bar/histogram edges clear where supported.
    graphicTypes = {'Bar','Histogram','ConstantLine'};
    for it = 1:numel(graphicTypes)
        objList = findall(figHandle,'Type',graphicTypes{it});
        for io = 1:numel(objList)
            try
                set(objList(io),'LineWidth',2);
            catch
            end
        end
    end

    % Style legends and colorbar labels.
    legList = findall(figHandle,'Type','Legend');
    for il = 1:numel(legList)
        try
            set(legList(il),'FontName','Times New Roman', ...
                'FontSize',14,'FontWeight','bold','Box','on');
        catch
        end
    end

    cbList = findall(figHandle,'Type','ColorBar');
    for ic = 1:numel(cbList)
        try
            set(cbList(ic),'FontName','Times New Roman', ...
                'FontSize',14,'FontWeight','bold','LineWidth',2);
            cbList(ic).Label.FontName = 'Times New Roman';
            cbList(ic).Label.FontSize = 14;
            cbList(ic).Label.FontWeight = 'bold';
        catch
        end
    end

    txtList = findall(figHandle,'Type','Text');
    for itxt = 1:numel(txtList)
        try
            set(txtList(itxt),'FontName','Times New Roman','FontWeight','bold');
            if txtList(itxt).FontSize < 14
                txtList(itxt).FontSize = 14;
            end
        catch
        end
    end
end



%% Independent pairwise coordinate audit (not solver pass flag)
function [maxOverlap,minGap,nViolations] = C35_coordinate_audit(T,Lx,Ly,tol)
    n=height(T); maxOverlap=0; minGap=NaN; nViolations=0;
    if n<2, return; end
    assert(all(isfinite(T.x)) && all(isfinite(T.y)) && ...
        all(isfinite(T.radius)) && all(T.radius>0),'Invalid final coordinates');
    minGap=inf;
    for i=1:n-1
        dx=mod(abs(T.x(i)-T.x(i+1:end)),Lx);
        dy=mod(abs(T.y(i)-T.y(i+1:end)),Ly);
        dx=min(dx,Lx-dx); dy=min(dy,Ly-dy);
        gap=hypot(dx,dy)-T.radius(i)-T.radius(i+1:end);
        minGap=min(minGap,min(gap));
        maxOverlap=max(maxOverlap,max(0,-min(gap)));
        nViolations=nViolations+nnz(gap < -tol);
    end
end
