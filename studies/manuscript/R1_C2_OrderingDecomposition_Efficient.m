%% R1_C2_OrderingDecomposition.m
% Scientific Reports, Reviewer 1, Comment 2.
% PURPOSE: separate density-driven local packing order from any measurable
% additional order associated with largest-void-guided insertion.
%
% Three MATCHED-PHYSICS generation arms are compared at EXACTLY the same
% integer particle count N and final monodisperse radius:
%  (1) Full_VGPP        = largest-void priority + weighted admissible mask;
%  (2) WeightedMask_Only= same admissible mask without local-maximum priority;
%  (3) UniformDomain_Growth = full-domain proposals + identical densification.
% A fourth arm, HardDiskMC_fromFull, measures how the Full_VGPP structure
% evolves under periodic nonoverlap-preserving Monte Carlo moves at FIXED N
% and r. It is a decorrelation check, NOT a proved equilibrium ensemble.
%
% Key measurements on successful, exactly density-matched configurations:
%  * |psi6|_local: mean magnitude of the sixfold complex order parameter
%    defined using PERIODIC DELAUNAY bond neighbors; range [0,1].
%  * |psi6|_global: magnitude of the particle-average complex psi6;
%  * fraction of particles with |psi6_i| > 0.7;
%  * g(r) based on exact minimum-image distances with finite-N correction;
%  * near-contact coordination Z(1.1 d) and first-shell g(r) maximum.
% Failure to reach the EXACT target is recorded as failure and not used
% as a density-matched ordering observation.
%
% Original packing kernels copied from uploaded code_refined(5).m via the
% earlier standalone R1_C1 comparator. No external source files required.
% Corrected stable weighted sampler is shared by all weighted variants.
% MATLAB R2024b with Image Processing Toolbox. No Statistics Toolbox needed.
%
% Output folder: revision_results/R1_C2_efficient (CSVs + cached coords)
% Picture folder: pics/R1_C2_efficient
% Prior Comment 1 raw CSVs do NOT hold particle coordinates; cannot reuse
% them to reconstruct psi6. Saved-coordinate files can be reused if present.
% Numerical interpretations MUST wait for actual MATLAB execution.

clc; clear; close all;
resultsDir = fullfile('revision_results','R1_C2_efficient'); picsDir = fullfile('pics','R1_C2_efficient');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end
if ~exist(picsDir,'dir'), mkdir(picsDir); end
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',13);

%% Physical and algorithmic inputs inherited from the uploaded base code
params = struct();
params.Lx = 10.0; params.Ly = 10.0;
params.Nx = 250; params.Ny = 250;
params.particleMode = "monodisperse";
params.rMono = 0.18; params.rMean = 0.18;
params.rStd = 0.035; params.rMin = 0.10; params.rMax = 0.28;
params.minGap = 0;
params.batchSize = 80; params.nCandidates = 4000;
params.distanceUpdateInterval = 5;
params.useLargestVoidInsertion = true;
params.useWeightedSampling = true;
params.largestVoidFraction = 0.70;
params.largestVoidMinSeparation = 0.75*params.rMono;
params.maxCycles = 10000; params.maxFailedCycles = 500;
params.relaxationMovesPerCycle = 5000;
params.growthCompressionRelaxIterations = 800;
params.growthCompressionRelaxFactor = 0.45;
params.overlapAuditTolerance = 1.0e-8;
params.growthCompressionOverlapTolerance = params.overlapAuditTolerance;
params.fixedRadiusSeedVF = 0.50;
params.fixedRadiusGrowthScale = 1.0025;
params.fixedRadiusMaxGrowthSteps = 300;
params.fixedRadiusMaxBacktracks = 14;
params.fixedRadiusMinGrowthScale = 1.00002;
params.targetVFTolerance = 1e-6;
params.useHexagonalFallback = false;
params.useGrowthCompression = false;
params.stopAtSpecifiedTarget = false;


%% Experiment configuration
% Two informative densities only: moderate and high. Increase after inspection.
fractions = [0.60 0.65 0.70];
nRuns = 2;                         % only 12 generator runs, not 360
quickTestMode = false;              % smoke: 3 arms x 1 seed at phi=0.60
if quickTestMode
    fractions = 0.60; nRuns = 1;
end
% IMPORTANT: The uploaded R1_C1_raw.csv has aggregate descriptors but NO
% saved coordinates. Sixfold order and g(r) CANNOT be reconstructed from it.
% If you have saved particle-center CSVs, put them into cacheDir with the
% naming pattern shown below. Otherwise generateMissing = true is needed.
generateMissing = true;           % set false for cached-coordinate analysis only
cacheDir = fullfile(resultsDir,'coordinate_cache');
if ~exist(cacheDir,'dir'), mkdir(cacheDir); end
seedBase = 9200;
methodNames = ["Full_VGPP","WeightedMask_Only","UniformDomain_Growth"];
mcName = "HardDiskMC_fromFull";

% Hard disk Monte Carlo is a finite-chain decorrelation reference only.
% Short MC chains cannot substantiate equilibrium hard-disk statistics.
runHardDiskMC = false; % optional; expensive and not an equilibrium reference
mcBurnSweeps = 20;                   % burn-in, each = N attempts
mcProductionSweeps = 20;             % total further sweeps
mcSamples = 2;                       % two production snapshots
mcDisplacement = 0.30*params.rMono;  % |delta_x|, |delta_y| bound
mcMinTravelDiameters = 0.50;         % descriptive mixing warning only
saveFirstExampleEachTarget = true;

% Radial bins in units of particle DIAMETER d=2*r. All below half-box.
d = 2*params.rMono;
rEdges = linspace(0,5*d,101);
if rEdges(end)>min(params.Lx,params.Ly)/2
    error('RDF radius must not exceed half the shorter box side.');
end

%% Main study: one row per arm, area fraction and independent seed
raw = table(); rdfRaw = table();
for iv = 1:numel(fractions)
    nominalVF = fractions(iv);
    Ntarget = round(nominalVF*params.Lx*params.Ly/(pi*params.rMono^2));
    discreteVF = Ntarget*pi*params.rMono^2/(params.Lx*params.Ly);
    fprintf('\n========== R1 C2: nominal phi = %.3f, exact phi = %.7f, N = %d ==========\n', ...
        nominalVF,discreteVF,Ntarget);

    for iRun = 1:nRuns
        thisSeed = seedBase+100*iv+iRun;
        fullCoordinates = [];

        for im = 1:numel(methodNames)
            method = methodNames(im);
            p = params;
            p.targetParticleCount = Ntarget;
            p.activeTargetVF = discreteVF;
            p.seedPlacementMode = 'full_vgpp';
            if method=="WeightedMask_Only"
                p.seedPlacementMode = 'mask_only';
                p.useLargestVoidInsertion = false;
                p.useWeightedSampling = true;
            elseif method=="UniformDomain_Growth"
                p.seedPlacementMode = 'uniform_domain';
                p.useLargestVoidInsertion = false;
                p.useWeightedSampling = false;
            end
            rng(thisSeed,'twister');
            status = ""; valid = false;
            nFinal=NaN; finalVF=NaN; elapsed=NaN; minClearance=NaN;
            psiLocal=NaN; psiGlobal=NaN; fracOrdered=NaN;
            gPeak=NaN; nearZ=NaN; rdf=[];
            completed = false;
            tRun = tic;
            try
                cacheFile = fullfile(cacheDir,sprintf('R1_C2_phi%02d_seed%d_%s.csv', ...
                    round(100*nominalVF),thisSeed,char(method)));
                if isfile(cacheFile)
                    T=readtable(cacheFile);
                    elapsed=NaN; % generation time unavailable for cached runs
                    status="Loaded cached particle coordinates";
                    completed=true;
                elseif generateMissing
                    result=run_fixed_final_radius_VGPP(p);
                    elapsed=toc(tRun);
                    T=result.particleTable;
                    completed=true;
                    if result.reachedTarget && height(T)==Ntarget
                        writetable(T,cacheFile); % reusable across later analyses
                    end
                else
                    error('R1C2:MissingCoordinates', ...
                        'No coordinate cache for %s, phi %.2f, seed %d',char(method),nominalVF,thisSeed);
                end
                requiredFields={'x','y','radius'};
                if ~all(ismember(requiredFields,T.Properties.VariableNames))
                    error('R1C2:BadCoordinateFile','CSV needs x, y, radius columns.');
                end
                nFinal=height(T); finalVF=sum(pi*T.radius.^2)/(p.Lx*p.Ly);
                [hasOverlap,~,minClearance,nBad]=final_overlap_audit_PBC(T.x,T.y,T.radius,p);
                valid = nFinal==Ntarget && ~hasOverlap && nBad==0 && ...
                    all(abs(T.radius-params.rMono)<1e-12) && ...
                    abs(finalVF-discreteVF)<1e-10;
                if valid
                    M=R1C2_measure_order(T,p,rEdges);
                    psiLocal=M.psi6Local; psiGlobal=M.psi6Global;
                    fracOrdered=M.fracStrongLocal; gPeak=M.gFirstPeak;
                    nearZ=M.nearCoord; rdf=M.rdf;
                    if method=="Full_VGPP"
                        fullCoordinates=T;
                    end
                    if saveFirstExampleEachTarget && iRun==1
                        name=sprintf('R1_C2_%.0f_%s_positions.csv',100*nominalVF,char(method));
                        writetable(T,fullfile(resultsDir,name));
                    end
                else
                    status="Did not pass exact-N / exact-radius / nonoverlap audit";
                end
            catch ME
                elapsed=toc(tRun);
                status=string(ME.identifier)+": "+string(ME.message);
                fprintf(2,'ERROR method %s phi=%.3f seed=%d: %s\n', ...
                    method, nominalVF,thisSeed,ME.message);
            end
            mAccept = NaN; mcTravel=NaN; mcLimitedMix=false;
            raw=[raw; R1C2_make_row(method,nominalVF,discreteVF,Ntarget,...
                iRun,thisSeed,completed,valid,nFinal,finalVF,elapsed,...
                psiLocal,psiGlobal,fracOrdered,gPeak,nearZ,minClearance,...
                mAccept,mcTravel,mcLimitedMix,status)]; %#ok<AGROW>
            if valid
                rdfRaw=[rdfRaw;R1C2_rdf_rows(method,nominalVF,discreteVF,...
                    iRun,thisSeed,rEdges,rdf,p.rMono)]; %#ok<AGROW>
            end
            fprintf('%-23s phi %.3f seed %d, valid %d, psi6 %.4f, runtime %.1fs\n', ...
                method,nominalVF,thisSeed,valid,psiLocal,elapsed);
        end

        % Only exact-N, fully valid VGPP packings may initialize the MC arm.
        mcValid=false; mcStatus="VGPP seed not valid";
        mcPsiLocal=NaN; mcPsiGlobal=NaN; mcFracOrdered=NaN;
        mcPeak=NaN; mcCoord=NaN; mcTime=NaN; mcAccept=NaN;
        mcTravel=NaN; mcLimitedMix=false; mcMinClearance=NaN;
        if runHardDiskMC && ~isempty(fullCoordinates)
            rng(thisSeed+200000,'twister');
            tMC=tic;
            try
                [mcInfo,mcSnapshots]=R1C2_run_mc(fullCoordinates,params,...
                    mcBurnSweeps,mcProductionSweeps,mcSamples,mcDisplacement);
                mcTime=toc(tMC);
                mcAccept=mcInfo.acceptanceRate;
                mcTravel=mcInfo.netTravelDiameters;
                mcLimitedMix=mcTravel<mcMinTravelDiameters;
                allPsiLocal=nan(numel(mcSnapshots),1);
                allPsiGlobal=nan(numel(mcSnapshots),1);
                allFrac=nan(numel(mcSnapshots),1);
                allPeak=nan(numel(mcSnapshots),1);
                allCoord=nan(numel(mcSnapshots),1);
                allRdf=nan(numel(mcSnapshots),numel(rEdges)-1);
                mcValid=true;
                for k=1:numel(mcSnapshots)
                    Ts=mcSnapshots{k};
                    [bad,~,mcMinClearance,nOverlaps]=final_overlap_audit_PBC(...
                        Ts.x,Ts.y,Ts.radius,params);
                    if bad || nOverlaps~=0
                        mcValid=false; mcStatus="Post-MC overlap audit failed";
                        break;
                    end
                    M=R1C2_measure_order(Ts,params,rEdges);
                    allPsiLocal(k)=M.psi6Local;
                    allPsiGlobal(k)=M.psi6Global;
                    allFrac(k)=M.fracStrongLocal;
                    allPeak(k)=M.gFirstPeak;
                    allCoord(k)=M.nearCoord;
                    allRdf(k,:)=M.rdf(:)';
                    rdfRaw=[rdfRaw;R1C2_rdf_rows(mcName,nominalVF,...
                        discreteVF,iRun,thisSeed,rEdges,M.rdf,params.rMono)]; %#ok<AGROW>
                end
                if mcValid
                    mcPsiLocal=mean(allPsiLocal,'omitnan');
                    mcPsiGlobal=mean(allPsiGlobal,'omitnan');
                    mcFracOrdered=mean(allFrac,'omitnan');
                    mcPeak=mean(allPeak,'omitnan');
                    mcCoord=mean(allCoord,'omitnan');
                    mcStatus="";
                    if mcLimitedMix
                        mcStatus="Limited MC displacement; NOT equilibrium-validated";
                    end
                    if saveFirstExampleEachTarget && iRun==1
                        name=sprintf('R1_C2_%.0f_MC_positions.csv',100*nominalVF);
                        writetable(mcSnapshots{end},fullfile(resultsDir,name));
                    end
                end
            catch ME
                mcTime=toc(tMC);
                mcValid=false;
                mcStatus=string(ME.identifier)+": "+string(ME.message);
                fprintf(2,'MC ERROR phi %.3f seed %d: %s\n',... 
                    nominalVF,thisSeed,ME.message);
            end
        end
        mcN=NaN; mcVF=NaN;
        if ~isempty(fullCoordinates),mcN=Ntarget;mcVF=discreteVF;end
        if runHardDiskMC
        raw=[raw;R1C2_make_row(mcName,nominalVF,discreteVF,Ntarget,...
            iRun,thisSeed,runHardDiskMC&&~isempty(fullCoordinates),mcValid,...
            mcN,mcVF,mcTime,...
            mcPsiLocal,mcPsiGlobal,mcFracOrdered,...
            mcPeak,mcCoord,mcMinClearance,mcAccept,mcTravel,mcLimitedMix,...
            mcStatus)]; %#ok<AGROW>
        end
        % Incremental checkpoint even during long runs.
        writetable(raw,fullfile(resultsDir,'R1_C2_raw.csv'));
        if ~isempty(rdfRaw)
            writetable(rdfRaw,fullfile(resultsDir,'R1_C2_rdf_individual.csv'));
        end
    end
end

%% Summary of FULLY VALID exact-density results (failures never silently dropped)
summary=table(); allNames=methodNames;
if runHardDiskMC, allNames=[methodNames,mcName]; end
for iv=1:numel(fractions)
    for im=1:numel(allNames)
        m=allNames(im);
        x=raw(raw.Method==m & abs(raw.NominalVF-fractions(iv))<1e-12,:);
        y=x(x.ValidAtTarget,:); nValid=height(y);
        row=table(m,fractions(iv),height(x),nValid,...
            R1C2_nanmean(y.Psi6Local),R1C2_nanstd(y.Psi6Local),...
            R1C2_nanmean(y.Psi6Global),R1C2_nanstd(y.Psi6Global),...
            R1C2_nanmean(y.FractionStrongLocal),R1C2_nanmean(y.GFirstPeak),...
            R1C2_nanmean(y.NearCoord),R1C2_nanmean(y.Runtime_s),...
            R1C2_nanmean(y.MCAcceptance),R1C2_nanmean(y.MCTravelDiameters),...
            nnz(y.MCLimitedMix),...
            'VariableNames',{'Method','NominalVF','NAttempted','NValid',...
            'Psi6LocalMean','Psi6LocalSD','Psi6GlobalMean','Psi6GlobalSD',...
            'FractionStrongMean','FirstPeakMean','NearCoordMean',...
            'RuntimeMean_s','MCAcceptMean','MCTravelMean','NMCWithLimitedTravel'});
        summary=[summary;row]; %#ok<AGROW>
    end
end
writetable(summary,fullfile(resultsDir,'R1_C2_summary.csv'));

%% Matched-seed differences, only pairs BOTH valid at SAME N and radius
paired=table();
rng(450001,'twister'); % fixed bootstrap stream, separate from generators
for iv=1:numel(fractions)
    full=raw(raw.Method=="Full_VGPP" & ...
             abs(raw.NominalVF-fractions(iv))<1e-12,:);
    comparatorNames=["WeightedMask_Only","UniformDomain_Growth"];
    if runHardDiskMC, comparatorNames=[comparatorNames,mcName]; end
    for comparison=comparatorNames
        alt=raw(raw.Method==comparison & ...
            abs(raw.NominalVF-fractions(iv))<1e-12,:);
        [~,ia,ib]=intersect(full.Seed,alt.Seed,'stable');
        both=full.ValidAtTarget(ia) & alt.ValidAtTarget(ib);
        ia=ia(both);ib=ib(both);
        delta=full.Psi6Local(ia)-alt.Psi6Local(ib);
        [lowerCI,upperCI]=R1C2_bootstrap_ci(delta,1000);
        paired=[paired;table(comparison,fractions(iv),numel(delta),...
            R1C2_nanmean(delta),lowerCI,upperCI,...
            R1C2_nanmean(full.Psi6Global(ia)-alt.Psi6Global(ib)),...
            R1C2_nanmean(full.GFirstPeak(ia)-alt.GFirstPeak(ib)),...
            'VariableNames',{'Comparator','NominalVF','NValidPairs',...
            'DeltaPsi6LocalMean','DeltaPsi6Local_CI95Lower',...
            'DeltaPsi6Local_CI95Upper','DeltaPsi6GlobalMean',...
            'DeltaFirstPeakMean'})]; %#ok<AGROW>
    end
end
writetable(paired,fullfile(resultsDir,'R1_C2_matched_effects.csv'));

%% Ensemble RDF summary. MC has several snapshots/seed: first average per seed.
rdfSummary=table();
if ~isempty(rdfRaw)
    for iv=1:numel(fractions)
        for im=1:numel(allNames)
            m=allNames(im);
            sub=rdfRaw(rdfRaw.Method==m & ...
                abs(rdfRaw.NominalVF-fractions(iv))<1e-12,:);
            if isempty(sub),continue;end
            % Equal weight to each independent generation seed.
            uniqSeeds=unique(sub.Seed);
            centers=unique(sub.RadiusOverDiameter);
            for k=1:numel(centers)
                part=sub(abs(sub.RadiusOverDiameter-centers(k))<1e-10,:);
                perSeed=nan(numel(uniqSeeds),1);
                for is=1:numel(uniqSeeds)
                    perSeed(is)=R1C2_nanmean(part.G(part.Seed==uniqSeeds(is)));
                end
                rdfSummary=[rdfSummary;table(m,fractions(iv),centers(k),...
                    numel(uniqSeeds),R1C2_nanmean(perSeed),...
                    R1C2_nanstd(perSeed),...
                    'VariableNames',{'Method','NominalVF',...
                    'RadiusOverDiameter','NSeeds','MeanG','SDG'})]; %#ok<AGROW>
            end
        end
    end
end
writetable(rdfSummary,fullfile(resultsDir,'R1_C2_rdf_summary.csv'));

%% Figures: ordering, paired differences, and RDF
% MATLAB R2024b
% Uses existing summary, paired, rdfSummary, raw, and picsDir variables.
% No simulations are repeated.

methodIDs = ["Full_VGPP", ...
             "WeightedMask_Only", ...
             "UniformDomain_Growth"];

methodLabels = ["Full VGPP", ...
                "Weighted sampling", ...
                "Uniform sampling"];

colors = [0.0000 0.4470 0.7410;
          0.8500 0.3250 0.0980;
          0.9290 0.6940 0.1250];

fontName = 'Times New Roman';
fontSize = 18;

if ~exist(picsDir,'dir')
    mkdir(picsDir);
end

%% Figure 1: Local sixfold ordering versus area fraction

fig = figure('Color','w','Position',[100 100 1100 620]);
hold on; box on;

for im = 1:numel(methodIDs)

    sub = summary(summary.Method == methodIDs(im), :);

    if isempty(sub)
        continue;
    end

    sub = sortrows(sub,'NominalVF');

    plotSD = sub.Psi6LocalSD;
    plotSD(~isfinite(plotSD)) = 0;

    errorbar(sub.NominalVF,sub.Psi6LocalMean,plotSD,...
        '-o',...
        'Color',colors(im,:),...
        'MarkerFaceColor','w',...
        'MarkerSize',8,...
        'LineWidth',2,...
        'CapSize',10,...
        'DisplayName',char(methodLabels(im)));

end

xlabel('Nominal inclusion area fraction');
ylabel('Mean local sixfold order, \langle|\psi_6|\rangle',...
       'Interpreter','tex');

legend('Location','northwest','Box','on');

set(gca,...
    'FontName',fontName,...
    'FontSize',fontSize,...
    'LineWidth',2,...
    'TickDir','in');

xlim([0.595 0.705]);
ylim([0.40 0.52]);

exportgraphics(fig,...
    fullfile(picsDir,'R1_C2_LocalPsi6_vs_Density.png'),...
    'Resolution',300);

%% Figure 2: Paired differences in local sixfold ordering
% Individual seed differences are shown because n=2 is insufficient
% for reliable inference from bootstrap confidence intervals.

fig = figure('Color','w','Position',[100 100 1100 620]);
hold on; box on;

comparisonIDs = ["WeightedMask_Only","UniformDomain_Growth"];
comparisonLabels = ["Weighted sampling","Uniform sampling"];

for im = 1:numel(comparisonIDs)

    sub = paired(paired.Comparator == comparisonIDs(im),:);

    if isempty(sub)
        continue;
    end

    sub = sortrows(sub,'NominalVF');

    plot(sub.NominalVF,sub.DeltaPsi6LocalMean,...
        '-o',...
        'Color',colors(im+1,:),...
        'MarkerSize',9,...
        'MarkerFaceColor','w',...
        'LineWidth',2,...
        'DisplayName',char(comparisonLabels(im)));

end

yline(0,'k--','LineWidth',1.2,...
    'HandleVisibility','off');

xlabel('Nominal inclusion area fraction');
ylabel('\Delta\langle|\psi_6|\rangle: Full VGPP - comparator',...
       'Interpreter','tex');

legend('Location','northwest','Box','on');

set(gca,...
    'FontName',fontName,...
    'FontSize',fontSize,...
    'LineWidth',2,...
    'TickDir','in');

xlim([0.595 0.705]);

exportgraphics(fig,...
    fullfile(picsDir,'R1_C2_PairedPsi6_Effects.png'),...
    'Resolution',300);

%% Figure 3: RDF at each investigated area fraction

rdfFractions = [0.60 0.65 0.70];

for iv = 1:numel(rdfFractions)

    vfToShow = rdfFractions(iv);

    fig = figure('Color','w','Position',[100 100 1100 620]);
    hold on; box on;

    nPlotted = 0;

    for im = 1:numel(methodIDs)

        sub = rdfSummary( ...
            rdfSummary.Method == methodIDs(im) & ...
            abs(rdfSummary.NominalVF-vfToShow)<1e-12,:);

        if isempty(sub)
            continue;
        end

        sub = sortrows(sub,'RadiusOverDiameter');

        plot(sub.RadiusOverDiameter,sub.MeanG,...
            'Color',colors(im,:),...
            'LineWidth',2,...
            'DisplayName',char(methodLabels(im)));

        nPlotted = nPlotted+1;
    end

    if nPlotted == 0
        warning('No RDF results for phi_A = %.2f',vfToShow);
        close(fig);
        continue;
    end

    xline(1,'k--','LineWidth',1.2,...
        'HandleVisibility','off');

    xlabel('Normalized separation, r/d');
    ylabel('Radial distribution function, g(r)');

    legend('Location','northeast','Box','on');

    set(gca,...
        'FontName',fontName,...
        'FontSize',fontSize,...
        'LineWidth',2,...
        'TickDir','in');

    xlim([0.8 5]);
    ylim([0 inf]);

    exportgraphics(fig,...
        fullfile(picsDir,...
        sprintf('R1_C2_RDF_%.0fPercent.png',100*vfToShow)),...
        'Resolution',300);

end

fprintf('\nR1 C2 figures regenerated from existing data.\n');

%% Helpers specific to this comment
function row=R1C2_make_row(method,nominalVF,discreteVF,Ntarget,run,seed,...
    completed,valid,nFinal,finalVF,elapsed,psiLocal,psiGlobal,fracStrong,...
    gPeak,nearZ,minClearance,mcAccept,mcTravel,mcLimitedMix,message)
row=table(method,nominalVF,discreteVF,Ntarget,run,seed,completed,valid,...
    nFinal,finalVF,elapsed,psiLocal,psiGlobal,fracStrong,gPeak,nearZ,...
    minClearance,mcAccept,mcTravel,mcLimitedMix,message,...
    'VariableNames',{'Method','NominalVF','ExactVF','NTarget','Run','Seed',...
    'Completed','ValidAtTarget','FinalN','FinalVF','Runtime_s','Psi6Local',...
    'Psi6Global','FractionStrongLocal','GFirstPeak','NearCoord',...
    'MinClearance','MCAcceptance','MCTravelDiameters','MCLimitedMix','Message'});
end

function row=R1C2_rdf_rows(method,nominalVF,exactVF,run,seed,edges,rdf,rMono)
    rCenter=(edges(1:end-1)+edges(2:end))/2;
    rCenter=rCenter/(2*rMono);
    n=numel(rCenter);
    row=table(repmat(method,n,1),repmat(nominalVF,n,1),...
        repmat(exactVF,n,1),repmat(run,n,1),repmat(seed,n,1),...
        rCenter(:),rdf(:),'VariableNames',...
        {'Method','NominalVF','ExactVF','Run','Seed','RadiusOverDiameter','G'});
end

function M=R1C2_measure_order(T,p,rEdges)
    n=height(T);
    if n<4, error('R1C2:TooFewParticles','Need >= 4 particles.');end
    x=T.x(:); y=T.y(:); Lx=p.Lx; Ly=p.Ly;
    % 3-by-3 periodic tiling ensures that central particle bonds through
    % the boundary have the correct neighbors. Delaunay triangulation does
    % NOT itself assume periodicity: it is applied to this explicit tiling.
    [SX,SY]=meshgrid(-1:1,-1:1);
    xx=repmat(x,9,1)+kron(SX(:),ones(n,1))*Lx;
    yy=repmat(y,9,1)+kron(SY(:),ones(n,1))*Ly;
    dt=delaunayTriangulation(xx,yy);
    EE=edges(dt);
    isCentral=@(ids) ids>4*n & ids<=5*n;
    ij=[EE(isCentral(EE(:,1)),:);fliplr(EE(isCentral(EE(:,2)),:))];
    ii=ij(:,1)-4*n;
    jj=mod(ij(:,2)-1,n)+1;
    good=ii~=jj;ijUnique=unique([ii(good),jj(good)],'rows');
    ii=ijUnique(:,1);jj=ijUnique(:,2);
    dx=x(jj)-x(ii);dy=y(jj)-y(ii);
    dx=dx-Lx*round(dx/Lx);dy=dy-Ly*round(dy/Ly);
    zz=exp(1i*6*atan2(dy,dx));
    nn=accumarray(ii,1,[n,1],@sum,0);
    if any(nn<3)
        error('R1C2:BadNeighbors','Periodic Delaunay: <3 bonds for a particle.');
    end
    ps=(accumarray(ii,real(zz),[n 1],@sum,0)+1i*...
        accumarray(ii,imag(zz),[n 1],@sum,0))./nn;
    M.psi6Local=mean(abs(ps));
    M.psi6Global=abs(mean(ps));
    M.fracStrongLocal=mean(abs(ps)>0.7);

    % Every unordered pair exactly once: PBC min-image metric.
    [I,J]=find(triu(true(n),1));
    px=x(J)-x(I);py=y(J)-y(I);
    px=px-Lx*round(px/Lx);py=py-Ly*round(py/Ly);
    rr=hypot(px,py);
    c=histcounts(rr,rEdges);
    areaAnn=pi*(rEdges(2:end).^2-rEdges(1:end-1).^2);
    expected=0.5*n*(n-1)*areaAnn/(Lx*Ly); % finite-N correction
    M.rdf=c./expected;
    mids=(rEdges(1:end-1)+rEdges(2:end))/2;
    shell=mids>=2*p.rMono & mids<=3*p.rMono;
    tmp=M.rdf(shell);
    M.gFirstPeak=max(tmp);
    M.nearCoord=2*nnz(rr<=2.2*p.rMono)/n;
end

function [info,snapshots]=R1C2_run_mc(T,p,nBurn,nProduction,nSamples,step)
    if nProduction<nSamples || mod(nProduction,nSamples)~=0
        error('R1C2:MCSettings','Production sweeps must divide equally among samples.');
    end
    x=T.x(:);y=T.y(:);x0=x;y0=y;
    n=numel(x);r=p.rMono;minD=2*r+p.minGap;
    thresh=max(minD-p.overlapAuditTolerance,0)^2;
    snapshots=cell(nSamples,1);
    accepted=0;attempted=0;snapNum=0;
    for sweep=1:(nBurn+nProduction)
        order=randperm(n);
        for idx=1:n
            i=order(idx);
            nx=mod(x(i)+(2*rand-1)*step,p.Lx);
            ny=mod(y(i)+(2*rand-1)*step,p.Ly);
            dx=x-nx; dy=y-ny;
            dx=dx-p.Lx*round(dx/p.Lx);
            dy=dy-p.Ly*round(dy/p.Ly);
            dist2=dx.^2+dy.^2;
            dist2(i)=Inf;
            attempted=attempted+1;
            if all(dist2>=thresh)
                x(i)=nx; y(i)=ny; accepted=accepted+1;
            end
        end
        if sweep>nBurn && mod(sweep-nBurn,nProduction/nSamples)==0
            snapNum=snapNum+1;
            snapshots{snapNum}=table(T.ID,x,y,T.radius,...
                'VariableNames',{'ID','x','y','radius'});
        end
    end
    dx=x-x0;dy=y-y0;
    dx=dx-p.Lx*round(dx/p.Lx);dy=dy-p.Ly*round(dy/p.Ly);
    info.acceptanceRate=accepted/max(attempted,1);
    info.netTravelDiameters=mean(hypot(dx,dy))/(2*r);
end

function m=R1C2_nanmean(x)
    x=x(isfinite(x));
    if isempty(x),m=NaN;else,m=mean(x);end
end
function s=R1C2_nanstd(x)
    x=x(isfinite(x));
    if numel(x)<2,s=NaN;else,s=std(x);end
end
function [lo,hi]=R1C2_bootstrap_ci(delta,nBoot)
    delta=delta(isfinite(delta));n=numel(delta);
    if n<2,lo=NaN;hi=NaN;return;end
    means=nan(nBoot,1);
    for ib=1:nBoot,means(ib)=mean(delta(randi(n,n,1)));end
    means=sort(means);
    lo=R1C2_percentile_sorted(means,0.025);
    hi=R1C2_percentile_sorted(means,0.975);
end
function q=R1C2_percentile_sorted(sorted,p)
    idx=1+(numel(sorted)-1)*p;
    low=floor(idx);high=ceil(idx);alpha=idx-low;
    q=(1-alpha)*sorted(low)+alpha*sorted(high);
end

%% Copied base-generation routines (source noted above)

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

    reachedTarget = reachedFinalRadius && ~hasOverlap && numel(rList) >= nTarget;

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

function [cxList, cyList, rList, occCore, distOcc, nInserted, currentArea] = ...
    attempt_DTPAA_batch_insertion(cxList, cyList, rList, occCore, distOcc, ...
    currentArea, X, Y, pixelSize, params, updateDistanceMap)

    Lx = params.Lx;
    Ly = params.Ly;

    nInserted = 0;

    %% Trial particle radius
    rCurrent = generate_particle_radius(params);

    %% Update distance transform if needed
    % Uniform-domain proposals do not require a distance transform. Avoid
    % charging that comparator for a map it does not use.
    if updateDistanceMap && ~strcmp(params.seedPlacementMode,'uniform_domain')
        if any(occCore(:))
            distOcc = periodic_distance_transform(occCore, pixelSize);
        else
            distOcc = inf(size(occCore));
        end
    end

    %% Identify eligible proposals. Uniform arm intentionally uses the whole
    % periodic Cartesian grid (not the occupancy mask); every proposed
    % center must still pass the identical exact periodic non-overlap test.
    if strcmp(params.seedPlacementMode, 'uniform_domain')
        admissible = true(size(occCore));
    elseif any(occCore(:))
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

    if strcmp(params.seedPlacementMode, 'uniform_domain')
        % Full-domain trial/rejection arm, on the common candidate grid.
        selectedIndices = admissibleIndices(randperm(nAvailable, nCandidates));
    else
        selectedIndices = select_candidate_indices_LV_hybrid( ...
            admissibleIndices, distOcc, nCandidates, params, X, Y);
    end

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

    % Stable Efraimidis-Spirakis weighted sampling without replacement.
    % The power form U^(1/w) in the uploaded baseline underflows to zero
    % for the roughly 6e4 grid points used here. The algebraically
    % equivalent log-key form avoids this artifact. Apply identically to
    % all arms that use weighted proposals, including Full_VGPP.
    weights = weights ./ sum(weights);
    positive = weights > 0;
    keys = -inf(nTotal,1);
    u = max(rand(nnz(positive),1), realmin('double'));
    keys(positive) = log(u)./weights(positive);
    [~, order] = sort(keys, 'descend');
    selected = order(1:n);
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
