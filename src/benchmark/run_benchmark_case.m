function [row, particleTable, historicalParticleTable] = run_benchmark_case(params, method, targetVF, nTarget, it, irun)
%RUN_BENCHMARK_CASE Run one seeded benchmark realization and audit geometry.
seed = 7300 + 100*it + irun;
rng(seed);
A = params.Lx*params.Ly;
paramsRun = params;
paramsRun.activeTargetVF = nTarget*pi*params.rMono^2/A;
paramsRun.targetParticleCount = nTarget;
paramsRun.useHexagonalFallback = false;
paramsRun.stopIfSpecifiedVFNotReached = false;

startTime = tic;
if method == "VGPP"
    result = run_fixed_final_radius_VGPP(paramsRun);
    particleTable = result.particleTable;
    historicalParticleTable = result.historicalParticleTable;
    reached = logical(result.reachedTarget);
else
    rsaParams = make_rsa_params_from_vgpp_params(paramsRun);
    rsaParams.targetN = nTarget;
    rsaParams.targetVF = paramsRun.activeTargetVF;
    result = run_RSA_PBC(rsaParams);
    particleTable = result.particleTable;
    historicalParticleTable = particleTable;
    reached = height(particleTable) >= nTarget;
end
elapsed = toc(startTime);
validate_historical_coordinates(particleTable,historicalParticleTable,paramsRun);
% The audit is independent of the method's completion flag.
[invalid, maxOverlap, minClearance, violatingPairs] = ...
    final_overlap_audit_PBC(particleTable.x, particleTable.y, ...
                            particleTable.radius, paramsRun);
finalRadii = ~isempty(particleTable) && ...
             all(abs(particleTable.radius-params.rMono)<=1e-12);
geometryValid = ~invalid && finalRadii;
if method == "VGPP"
    geometryValid = geometryValid && logical(result.geometryValid);
end
reached = reached && geometryValid;

if ~geometryValid
    completion = "invalid_geometry";
elseif reached
    completion = "target_attained";
else
    completion = "valid_below_target";
end

% Invalid geometries remain visible in the raw audit, but cannot enter
% reported spatial statistics or the successful-target comparison.
meanRadius = NaN; localRSD=NaN; lowFrac=NaN;
meanNND=NaN; nndCV=NaN; rdfFirstPeak=NaN;
finalVF = NaN; radiusScale=NaN; canonicalRSD=NaN; canonicalLow=NaN;
if geometryValid
    stats = compute_spatial_statistics_PBC(particleTable, paramsRun, method);
    % Historical local cells use pre-canonicalization boundary-side storage.
    % Pair distances and all geometry continue to use particleTable.
    [lowFrac,localRSD] = local_vf_quality_metrics(historicalParticleTable,paramsRun);
    [canonicalLow,canonicalRSD] = local_vf_quality_metrics(particleTable,paramsRun);
    meanRadius=mean(particleTable.radius);
    radiusScale=meanRadius/params.rMono;
    finalVF=sum(pi*particleTable.radius.^2)/A;
    meanNND=stats.meanNND;
    nndCV=stats.nndCV;
    rdfFirstPeak=stats.rdfFirstPeak;
end
row = table(method,targetVF,paramsRun.activeTargetVF,nTarget,irun,seed, ...
    finalVF,double(reached),height(particleTable),elapsed,meanRadius, ...
    radiusScale,localRSD,lowFrac,meanNND,nndCV,minClearance, ...
    maxOverlap,violatingPairs,rdfFirstPeak,geometryValid,completion,canonicalRSD,canonicalLow, ...
    'VariableNames',{'Method','TargetVF','DiscreteTargetVF','TargetN', ...
    'Run','Seed','FinalVF','ReachedTarget','NumberOfParticles', ...
    'ElapsedTime_s','MeanRadius','RadiusScale','LocalVFRSD', ...
    'LowVFCellFraction','MeanNND','NNDCV','MinClearance', ...
    'MaxOverlap','OverlappingPairs','RDFFirstPeak','GeometryValid', ...
    'CompletionStatus','CanonicalLocalVFRSD','CanonicalLowVFCellFraction'});
end
