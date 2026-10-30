%% R1_C1_FocusedAblation.m
% MATLAB R2024b | Image Processing Toolbox needed for bwdist/imregionalmax
% Focused complementary study: no-radius-growth and no-seed-relaxation.
% Previously completed 360-run benchmark is reused, NOT rerun.
% CPU budget is cooperative, checked between/within solver operations.
% A timed-out run is classified CENSORED, NOT a packing failure.
clc; clear; close all;
quickTestMode = false;
outDir=fullfile('revision_results','R1_C1_focused_ablation');
if ~exist(outDir,'dir'), mkdir(outDir); end
baselineFile='R1_C1_raw.csv';
if ~isfile(baselineFile), baselineFile=fullfile('revision_results','R1_C1_raw.csv'); end
if isfile(baselineFile)
    baseline=readtable(baselineFile,'TextType','string');
    copyfile(baselineFile,fullfile(outDir,'R1_C1_existing_baseline_copy.csv'));
    fprintf('Existing benchmark reused: %d recorded runs.\n',height(baseline));
else
    warning('Previous 360-run CSV not found. New experiments will still run.');
end
params=struct();
params.Lx=10; params.Ly=10; params.Nx=250; params.Ny=250;
params.particleMode="monodisperse";
params.rMono=.18; params.rMean=.18; params.rStd=.035;
params.rMin=.10; params.rMax=.28; params.minGap=0;
params.batchSize=80; params.nCandidates=4000;
params.distanceUpdateInterval=5;
params.useLargestVoidInsertion=true; params.useWeightedSampling=true;
params.largestVoidFraction=.70;
params.largestVoidMinSeparation=.75*params.rMono;
params.maxCycles=10000; params.maxFailedCycles=500;
params.relaxationMovesPerCycle=5000;
params.growthCompressionRelaxIterations=800;
params.growthCompressionRelaxFactor=.45;
params.overlapAuditTolerance=1e-8;
params.growthCompressionOverlapTolerance=1e-8;
params.fixedRadiusSeedVF=.50;
params.fixedRadiusGrowthScale=1.0025;
params.fixedRadiusMaxGrowthSteps=300;
params.fixedRadiusMaxBacktracks=14;
params.fixedRadiusMinGrowthScale=1.00002;
params.targetVFTolerance=1e-6;
params.useHexagonalFallback=false;
params.useGrowthCompression=false;
params.stopAtSpecifiedTarget=false;
params.stopIfSpecifiedVFNotReached=false;
params.ablationDisableGrowth=false;
params.ablationDisableSeedRelax=false;
params.seedPlacementMode='full_vgpp';
if quickTestMode
    targetFractions=[.50 .60]; seeds=1; perRunLimit_s=20;
    params.maxCycles=500; params.maxFailedCycles=70;
else
    targetFractions=[.50 .60 .70]; seeds=1:3; perRunLimit_s=75;
    params.maxCycles=3000; params.maxFailedCycles=150;
end
methods=["NoRadiusGrowth","NoSeedStochasticRelax"];
raw=table();
rawPath=fullfile(outDir,'R1_C1_focused_raw.csv');
for it=1:numel(targetFractions)
  fraction=targetFractions(it);
  targetN=round(fraction*params.Lx*params.Ly/(pi*params.rMono^2));
  exactAF=targetN*pi*params.rMono^2/(params.Lx*params.Ly);
  for ir=seeds
    seed=7400+ir;  % matched with original benchmark seeds 7401--7420
    for im=1:numel(methods)
      method=methods(im); p=params;
      p.targetParticleCount=targetN; p.activeTargetVF=exactAF;
      p.ablationDisableGrowth=(method=="NoRadiusGrowth");
      p.ablationDisableSeedRelax=(method=="NoSeedStochasticRelax");
      rng(seed,'twister');
      global R1C1_START_TIC R1C1_MAX_SECONDS
      R1C1_START_TIC=tic; R1C1_MAX_SECONDS=perRunLimit_s;
      status="ERROR"; completed=false; reached=false; censored=false;
      seedN=NaN; finalN=NaN; removedN=NaN; growthSteps=NaN;
      minClearance=NaN; maxOverlap=NaN; localRSD=NaN; err="";
      fprintf('Starting %s, AF %.2f, seed %d (budget %.0fs)\n',method,fraction,seed,perRunLimit_s);
      try
        result=run_fixed_final_radius_VGPP(p);
        completed=true;
        seedN=result.seedN; finalN=height(result.particleTable);
        removedN=result.removedCount; growthSteps=result.growthSteps;
        minClearance=result.minClearance; maxOverlap=result.maxOverlap;
        reached=logical(result.reachedTarget) && finalN==targetN && ...
            result.nOverlappingPairs==0 && ...
            all(abs(result.particleTable.radius-p.rMono)<1e-12);
        if finalN>1, localRSD=R1C1_localVF_RSD(result.particleTable,p,8); end
        if reached, status="SUCCESS"; else, status="COMPLETED_BELOW_TARGET"; end
      catch ME
        if strcmp(ME.identifier,'R1C1:BudgetExceeded')
           censored=true; status="CENSORED_RUNTIME";
        else
           err=string(ME.message); status="ERROR";
        end
      end
      elapsed=toc(R1C1_START_TIC);
      row=table(method,fraction,exactAF,targetN,ir,seed,string(status),completed,reached,censored,...
          seedN,finalN,removedN,growthSteps,elapsed,minClearance,maxOverlap,localRSD,err,...
          'VariableNames',{'Method','TargetAF','DiscreteTargetAF','TargetN','Run','Seed','Status',...
          'Completed','ReachedTarget','Censored','SeedN','FinalN','RemovedN','GrowthSteps',...
          'ElapsedTime_s','MinClearance','MaxOverlap','LocalAFRSD','ErrorMessage'});
      raw=[raw;row]; %#ok<AGROW>
      writetable(raw,rawPath);
      fprintf('  %s, elapsed %.1fs, final N %.0f / %d\n',status,elapsed,finalN,targetN);
    end
  end
end
summary=table();
for im=1:numel(methods)
 for it=1:numel(targetFractions)
  rows=raw(raw.Method==methods(im) & abs(raw.TargetAF-targetFractions(it))<1e-10,:);
  n=height(rows); nSuccess=sum(rows.ReachedTarget); nCensored=sum(rows.Censored);
  nCompleted=sum(rows.Completed); nError=sum(rows.Status=="ERROR");
  % SuccessPctAll is the observed share, not an estimate of failure probability
  % because timeout runs are censored.
  obsPct=100*nSuccess/n;
  meanTime=mean(rows.ElapsedTime_s,'omitnan');
  medTime=median(rows.ElapsedTime_s,'omitnan');
  meanSeedN=mean(rows.SeedN,'omitnan');
  meanRemovedN=mean(rows.RemovedN,'omitnan');
  row=table(methods(im),targetFractions(it),n,nSuccess,nCompleted,nCensored,nError,...
    obsPct,meanTime,medTime,meanSeedN,meanRemovedN,...
    'VariableNames',{'Method','TargetAF','Runs','Success','Completed','Censored','Errors',...
    'ObservedSuccessPct','MeanElapsed_s','MedianElapsed_s','MeanSeedN','MeanRemovedN'});
  summary=[summary;row]; %#ok<AGROW>
 end
end
writetable(summary,fullfile(outDir,'R1_C1_focused_summary.csv'));
if exist('baseline','var')
  targetAF=unique(baseline.TargetVF);
  compact=table();
  for i=1:numel(targetAF)
    rows=baseline(abs(baseline.TargetVF-targetAF(i))<1e-10,:);
    for im=1:numel(unique(rows.Method))
      names=unique(rows.Method);
      sub=rows(rows.Method==names(im),:);
      compact=[compact;table(string(names(im)),targetAF(i),height(sub),sum(sub.ReachedTarget),...
        median(sub.ElapsedTime_s),'VariableNames',{'Method','TargetAF','Runs','Success','MedianElapsed_s'})]; %#ok<AGROW>
    end
  end
  writetable(compact,fullfile(outDir,'R1_C1_existing_baseline_summary.csv'));
end
fprintf('Finished focused ablation. Outputs: %s\n',outDir);

function R1C1_budgetCheck(~)
  global R1C1_START_TIC R1C1_MAX_SECONDS
  if ~isempty(R1C1_START_TIC) && toc(R1C1_START_TIC)>=R1C1_MAX_SECONDS
    error('R1C1:BudgetExceeded','Cooperative per-run CPU-time budget exceeded');
  end
end

function rsd = R1C1_localVF_RSD(T,p,bins)
    counts = histcounts2(T.y,T.x,linspace(0,p.Ly,bins+1),linspace(0,p.Lx,bins+1));
    localVF = counts*pi*p.rMono^2 / ((p.Lx/bins)*(p.Ly/bins));
    rsd=std(localVF(:))/max(mean(localVF(:)),eps);
end

function cv = R1C1_nearestCV(T,p)
    n=height(T); if n<2, cv=NaN;return;end
    nearest = nan(n,1);
    for i=1:n
        dx=abs(T.x-T.x(i)); dx=min(dx,p.Lx-dx);
        dy=abs(T.y-T.y(i)); dy=min(dy,p.Ly-dy);
        dd=hypot(dx,dy);dd(i)=Inf;nearest(i)=min(dd);
    end
    cv=std(nearest)/max(mean(nearest),eps);
end

function out=R1C1_wilson95(k,n)
    z=1.959963984540054;
    p=k/n;den=1+z*z/n;
    ctr=(p+z*z/(2*n))/den;
    half=z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den;
    out=[max(0,ctr-half),min(1,ctr+half)];
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

    R1C1_budgetCheck(params);
    Lx = params.Lx;
    Ly = params.Ly;
    A = Lx * Ly;
    rFinal = params.rMono;
    nTarget = params.targetParticleCount;
    targetVF = nTarget*pi*rFinal^2/A;

    % Keep initial packing at or below the selected seed VF.
    if isfield(params,'ablationDisableGrowth') && params.ablationDisableGrowth
        seedScale = 1.0;
    else
        seedScale = min(1.0, sqrt(params.fixedRadiusSeedVF / targetVF));
        seedScale = max(seedScale, 0.65);
    end
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

        R1C1_budgetCheck(params);
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
            if (~isfield(params,'ablationDisableSeedRelax') || ~params.ablationDisableSeedRelax) && ...
                    mod(failedCycles,10) == 0 && numel(rList) > 1
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

    seedN = numel(rList);
    seedCycles = cycles;
    seedFailedCycles = failedCycles;
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

        R1C1_budgetCheck(params);
        growthSteps = growthSteps + 1;

        rCurrent = mean(rList);
        trialScale = min(params.fixedRadiusGrowthScale, rFinal/rCurrent);
        accepted = false;

        for ib = 0:params.fixedRadiusMaxBacktracks
            R1C1_budgetCheck(params);
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
        R1C1_budgetCheck(params);
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

    result.seedN = seedN;
    result.seedCycles = seedCycles;
    result.seedFailedCycles = seedFailedCycles;
    result.removedCount = max(0,seedN-numel(rList));
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
        if mod(k,128)==1, R1C1_budgetCheck(params); end

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
        if mod(p,16)==1, R1C1_budgetCheck(params); end

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
        if mod(moveID,100)==1, R1C1_budgetCheck(params); end

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
        if mod(iter,5)==1, R1C1_budgetCheck(params); end

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
