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
