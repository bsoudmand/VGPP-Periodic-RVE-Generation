function [sweepSummary, plotData] = run_target_vf_sweep_block(params)
%RUN_TARGET_VF_SWEEP_BLOCK Execute or resume the fixed-final-radius benchmark.
targetList = params.targetVFSweepList(:);
nRuns = max(1,round(params.targetVFSweepRuns));
includeRSA = params.runRSAInTargetSweep;
resultsFolder = params.resultsFolder;
picsFolder = params.figuresFolder;
if ~isfolder(resultsFolder), mkdir(resultsFolder); end
if ~isfolder(picsFolder), mkdir(picsFolder); end

oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
VGPPParticleTables = cell(numel(targetList),nRuns);
VGPPReachedMatrix = false(numel(targetList),nRuns);
VGPPLocalRSDMatrix = nan(numel(targetList),nRuns);
VGPPSeedMatrix = nan(numel(targetList),nRuns);
records = cell(numel(targetList)*nRuns*(1+double(includeRSA)),1);
recordCount = 0;
A = params.Lx*params.Ly;
particleArea = pi*params.rMono^2;

for it=1:numel(targetList)
    targetVF = targetList(it);
    nTarget = max(1,round(targetVF*A/particleArea));
    fprintf('Target %.2f, %d particles\n',targetVF,nTarget);
    for irun=1:nRuns
        [v,p] = benchmark_case_checkpoint(params,"VGPP",targetVF,nTarget,it,irun);
        recordCount = recordCount+1;
        records{recordCount}=v;
        if v.GeometryValid && v.ReachedTarget
            VGPPParticleTables{it,irun}=p;
            VGPPReachedMatrix(it,irun)=true;
            VGPPLocalRSDMatrix(it,irun)=v.LocalVFRSD;
            VGPPSeedMatrix(it,irun)=v.Seed;
        end
        if includeRSA
            [a,~] = benchmark_case_checkpoint(params,"RSA-PBC",targetVF,nTarget,it,irun);
            recordCount = recordCount+1;
            records{recordCount}=a;
        end
        % Each completed realization is already saved separately.
        fprintf('  completed run %d/%d\n',irun,nRuns);
    end
end
sweepRaw = vertcat(records{1:recordCount});
writetable(sweepRaw,fullfile(resultsFolder,'Table_R1_fixedRadius_raw.csv'));

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

    plotData = struct('targetList',targetList,'sweepSummary',sweepSummary, ...
        'VGPPParticleTables',{VGPPParticleTables}, ...
        'VGPPReachedMatrix',VGPPReachedMatrix, ...
        'VGPPLocalRSDMatrix',VGPPLocalRSDMatrix, ...
        'VGPPSeedMatrix',VGPPSeedMatrix);
end
