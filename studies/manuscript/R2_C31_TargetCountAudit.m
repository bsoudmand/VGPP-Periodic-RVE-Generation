%% R2_C31_TargetCountAudit.m
% Reviewer 2 Comment 31: verify compatibility of particle-removal fallback
% with reported 100% exact-target attainment.
% MATLAB R2024b; requires the raw CSV produced by code_refined(5).m.
% Standalone audit script: no functions from the base source are required.
clc; clear; close all;

inputFile = fullfile('revision_results','Table_R1_fixedRadius_raw.csv');
outFolder = fullfile('revision_results','R2_C31');
if ~exist(outFolder,'dir'), mkdir(outFolder); end
if ~isfile(inputFile)
    error(['Required simulation data not found: %s\n' ...
       'First run the target-VF sweep in code_refined(5).m and copy its raw CSV here.'],inputFile);
end
T = readtable(inputFile,'TextType','string');
required = ["Method","TargetVF","DiscreteTargetVF","TargetN", ...
    "Run","Seed","FinalVF","ReachedTarget","NumberOfParticles", ...
    "MeanRadius","MaxOverlap","OverlappingPairs"];
missing = setdiff(required,string(T.Properties.VariableNames));
if ~isempty(missing)
    error('Missing CSV columns: %s',strjoin(missing,', '));
end
V = T(T.Method == "VGPP",:);
if isempty(V), error('No VGPP rows found in raw table.'); end

% Existing source uses a fixed prescribed radius and a 10x10 domain.
rFinal = 0.18; Lx = 10; Ly = 10;
tol = 1e-8; targetTol = 1e-9;
expectedDiscrete = V.TargetN .* pi*rFinal^2/(Lx*Ly);
countExact = V.NumberOfParticles == V.TargetN;
countBelow = V.NumberOfParticles < V.TargetN;
countAbove = V.NumberOfParticles > V.TargetN;
radiusExact = abs(V.MeanRadius-rFinal) <= 1e-10;
areaConsistency = abs(V.FinalVF - V.NumberOfParticles.*pi.*V.MeanRadius.^2/(Lx*Ly)) <= targetTol;
targetConsistency = abs(V.DiscreteTargetVF-expectedDiscrete) <= targetTol;
noRecordedOverlaps = (V.OverlappingPairs == 0) & (V.MaxOverlap <= tol);
% A record is eligible only if its target count, radius and audit all pass.
verifiedTarget = countExact & radiusExact & areaConsistency & ...
    targetConsistency & noRecordedOverlaps;
reported = logical(V.ReachedTarget);
disagreement = reported ~= verifiedTarget;

A = table(V.TargetVF,V.Run,V.Seed,V.TargetN,V.NumberOfParticles, ...
    V.DiscreteTargetVF,V.FinalVF,V.MeanRadius,V.MaxOverlap, ...
    V.OverlappingPairs,reported,countExact,countBelow,countAbove, ...
    radiusExact,areaConsistency,targetConsistency,noRecordedOverlaps, ...
    verifiedTarget,disagreement, ...
    'VariableNames',{'TargetAreaFraction','Run','Seed','TargetN', ...
    'FinalN','DiscreteTargetAreaFraction','FinalAreaFraction', ...
    'MeanFinalRadius','MaxOverlap','OverlappingPairs','ReportedReached', ...
    'CountExact','CountBelowTarget','CountAboveTarget','RadiusExact', ...
    'AreaFractionConsistent','TargetAreaFractionConsistent', ...
    'NoRecordedOverlaps','IndependentlyVerifiedTarget','FlagDisagreement'});
writetable(A,fullfile(outFolder,'R2_C31_run_audit.csv'));

u = unique(V.TargetVF);
S = table('Size',[numel(u) 9], 'VariableTypes', ...
    {'double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames',{'TargetAreaFraction','Runs','ReportedSuccesses', ...
    'VerifiedSuccesses','VerifiedRatePercent','CountBelowTarget', ...
    'CountAboveTarget','AuditDisagreements','MeanParticleShortfall'});
for k=1:numel(u)
    ii = V.TargetVF == u(k);
    S.TargetAreaFraction(k) = u(k);
    S.Runs(k) = nnz(ii);
    S.ReportedSuccesses(k) = nnz(reported(ii));
    S.VerifiedSuccesses(k) = nnz(verifiedTarget(ii));
    S.VerifiedRatePercent(k) = 100*mean(verifiedTarget(ii));
    S.CountBelowTarget(k) = nnz(countBelow(ii));
    S.CountAboveTarget(k) = nnz(countAbove(ii));
    S.AuditDisagreements(k) = nnz(disagreement(ii));
    S.MeanParticleShortfall(k) = mean(V.TargetN(ii)-V.NumberOfParticles(ii));
end
writetable(S,fullfile(outFolder,'R2_C31_target_summary.csv'));

set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
f = figure('Color','w','Position',[100 100 1000 520]);
bar(u,[S.VerifiedSuccesses S.Runs-S.VerifiedSuccesses],'stacked');
xlabel('Target particle area fraction'); ylabel('Realizations');
legend({'Exact target verified','Not verified'},'Location','best');
title('Independent exact-target audit'); grid on;
exportgraphics(f,fullfile(outFolder,'R2_C31_exact_target_counts.png'),'Resolution',300);

fprintf('\nReviewer 2 Comment 31: exact target audit\n');
disp(S);
fprintf('Rows with discrepancy between original flag and independent audit: %d\n',nnz(disagreement));
fprintf('Rows below target N (could include failed insertion or removal): %d\n',nnz(countBelow));
fprintf(['Note: the raw CSV does NOT record whether a particle was actually removed.\n' ...
 'A shortfall can arise from initial insertion failure or the removal fallback.\n' ...
 'Therefore this audit verifies the success claim, but not removal frequency.\n']);
if any(disagreement)
    warning('Reported reachability differs from independent audit in %d runs. Inspect CSV.',nnz(disagreement));
end
