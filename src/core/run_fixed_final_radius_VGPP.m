function result = run_fixed_final_radius_VGPP(params)

    % Fixed-final-radius VGPP
    % The requested particle number is first inserted with temporarily
    % shrunken particles. Particle centers are then collectively rearranged
    % while the radii are increased gradually, but never beyond rMono.
    % Therefore every successful final RVE has exactly the prescribed final
    % particle radius. No hexagonal/lattice fallback is used.

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

    % Preserve upper-boundary storage for the original cell statistics.
    % These values do not feed insertion, relaxation, or geometric audits.
    historicalX = cxList;
    historicalY = cyList;

    % The candidate grid includes the periodic upper endpoints. Represent
    % those exact endpoints by their identical zero-boundary locations.
    % Do not wrap arbitrary invalid or nonfinite coordinates: the audit
    % below must continue to reject them.
    cxList(cxList == Lx) = 0;
    cyList(cyList == Ly) = 0;

    % Final strict audit with the same tolerance.
    [hasOverlap,maxOverlap,minClearance,nOverlapPairs] = ...
        final_overlap_audit_PBC(cxList,cyList,rList,params);

    reachedTarget = reachedFinalRadius && ~hasOverlap && numel(rList) >= nTarget;

    particleID = (1:numel(rList))';
    particleTable = table(particleID,cxList,cyList,rList, ...
        'VariableNames',{'ID','x','y','radius'});

    result.particleTable = particleTable;
    result.historicalParticleTable = table(particleID,historicalX,historicalY,rList, ...
        'VariableNames',{'ID','x','y','radius'});
    result.maxAreaFraction = sum(pi*rList.^2)/A;
    result.reachedTarget = reachedTarget;
    result.geometryValid = reachedFinalRadius && ~hasOverlap && ...
        ~isempty(rList) && all(abs(rList-rFinal) <= 1e-12);
    result.targetN = nTarget;
    result.seedRadius = rSeed;
    result.seedScale = seedScale;
    result.growthSteps = growthSteps;
    result.hasOverlap = hasOverlap;
    result.maxOverlap = maxOverlap;
    result.minClearance = minClearance;
    result.nOverlappingPairs = nOverlapPairs;
end
