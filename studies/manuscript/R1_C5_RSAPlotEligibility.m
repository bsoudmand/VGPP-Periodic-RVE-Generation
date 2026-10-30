%% R1_C5_RSAPlotEligibility.m
% Scientific Reports - VGPP first revision
% Reviewer 1, Comment 5:
% "Clarify in figure captions that RSA PBC data are only plotted for
% volume fractions that RSA successfully achieves."
%
% PURPOSE
%   Rebuild manuscript Figures 1-4 from the ORIGINAL individual-realization data.
%   The distinction is important:
%     Fig. 1: All prescribed targets, including failed RSA-PBC targets,
%             because the figure explicitly demonstrates saturation/failure.
%     Figs. 2-4: RSA-PBC is shown ONLY where every RSA-PBC realization
%                reached its prescribed particle count (strict 100% rule).
%   VGPP statistics in Figs. 2-4 likewise require 100% target reachability.
%   The script DOES NOT simulate new microstructures or infer missing data.
%
% STANDALONE REQUIREMENT
%   No custom MATLAB functions or toolboxes are needed.
%   INPUT: revision_results/Table_R1_fixedRadius_raw.csv from the
%          multi-target sweep in code_refined(5).m.
%   If the original sweep is unfinished, run this script after that sweep
%   creates the CSV. Never substitute unverified or fictional values.
%
% OUTPUTS
%   revision_results/R1_C5_plot_eligibility.csv
%   revision_results/R1_C5_plot_summary.csv
%   pics/R1_C5_Fig1_reachability.png
%   pics/R1_C5_Fig2_spatial_heterogeneity.png
%   pics/R1_C5_Fig3_periodic_overlap_audit.png
%   pics/R1_C5_Fig4_ordering_runtime.png
%
% The files use a new R1_C5 prefix and do not overwrite previous figures.
% Figures/captions in the manuscript must NOT be replaced until the outputs
% have been inspected and verified.

clc;
clear;
close all;

%% --------------------------- USER INPUTS -------------------------------
rawCSV = fullfile('revision_results', 'Table_R1_fixedRadius_raw.csv');
resultsFolder = 'revision_results';
picsFolder = 'pics';
figureDPI = 300;
visibleFigures = 'off';       % 'on' for interactive visual inspection
strictAllSeeds = true;      % mandatory for main analysis
expectedNRuns = 20;        % original manuscript has 20 independent runs per target
expectedMethods = ["VGPP", "RSA-PBC"];

if ~exist(resultsFolder,'dir'), mkdir(resultsFolder); end
if ~exist(picsFolder,'dir'), mkdir(picsFolder); end

if ~isfile(rawCSV)
    alternativeCSV = 'Table_R1_fixedRadius_raw.csv';
    if isfile(alternativeCSV)
        rawCSV = alternativeCSV;
    else
        error(['Cannot find Table_R1_fixedRadius_raw.csv. ', ...
            'First run the target-volume-fraction sweep in code_refined(5).m ', ...
            '(runTargetVFSweepBlock = true). Place the resulting CSV in ', ...
            'revision_results or next to this script.']);
    end
end

set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',16);
set(groot,'defaultAxesLineWidth',2);
set(groot,'defaultLineLineWidth',2);

%% --------------------- LOAD AND VALIDATE RAW DATA ----------------------
T = readtable(rawCSV, 'TextType','string');
requiredColumns = { ...
    'Method','TargetVF','FinalVF','ReachedTarget', ...
    'ElapsedTime_s','LocalVFRSD','LowVFCellFraction', ...
    'NNDCV','RDFFirstPeak','MinClearance','OverlappingPairs'};
missingColumns = setdiff(requiredColumns, T.Properties.VariableNames);
if ~isempty(missingColumns)
    error('Raw CSV is missing the required columns: %s', ...
        strjoin(missingColumns,', '));
end

T.Method = upper(strtrim(string(T.Method)));
T.Method(contains(T.Method,'RSA')) = "RSA-PBC";
T.Method(contains(T.Method,'VGPP')) = "VGPP";

% Defensively handle common CSV type variations.
for j = 1:numel(requiredColumns)
    nm = requiredColumns{j};
    if strcmp(nm,'Method') || strcmp(nm,'ReachedTarget'), continue; end
    if ~isnumeric(T.(nm))
        T.(nm) = str2double(string(T.(nm)));
    end
end
if islogical(T.ReachedTarget)
    % Already correct.
elseif isnumeric(T.ReachedTarget)
    T.ReachedTarget = isfinite(T.ReachedTarget) & T.ReachedTarget > 0.5;
else
    v = lower(strtrim(string(T.ReachedTarget)));
    T.ReachedTarget = ismember(v, ["1","true","yes"]);
end

assert(height(T)>0,'Raw CSV has zero data rows.');
assert(all(isfinite(T.TargetVF)), 'Nonfinite target area fraction in the raw CSV.');
assert(all(isfinite(T.FinalVF)), 'Nonfinite achieved area fraction in the raw CSV.');
assert(all(isfinite(T.ElapsedTime_s) & T.ElapsedTime_s>=0), ...
    'The raw CSV contains missing or negative elapsed times.');
assert(all(ismember(T.Method,expectedMethods)), ...
    'Expected method labels VGPP and RSA-PBC; inspect Method column.');
assert(all(ismember(expectedMethods,unique(T.Method))), ...
    'One or both of VGPP / RSA-PBC are absent from the raw CSV.');

vfList = unique(T.TargetVF,'sorted');
nVF = numel(vfList);

% The default is intentionally immutable: do NOT use a 95% cutoff and
% treat a failed target as fully reached for geometric/statistical plots.
if ~strictAllSeeds
    error('For Comment 5, strictAllSeeds must remain true.');
end

%% -------------------- SUMMARIZE BY METHOD/TARGET ----------------------
Method = strings(0,1); TargetVF = zeros(0,1);
NRuns = zeros(0,1); NSuccess = zeros(0,1);
SuccessRate_percent = zeros(0,1); EligibleForFigs234 = false(0,1);
MeanFinalVF_All = zeros(0,1); StdFinalVF_All = zeros(0,1);
MeanRuntime_All_s = zeros(0,1);
MeanRuntime_Success_s = zeros(0,1); StdRuntime_Success_s = zeros(0,1);
MeanLocalVFRSD_Success = zeros(0,1); StdLocalVFRSD_Success = zeros(0,1);
MeanLowVFCellFraction_Success = zeros(0,1); StdLowVFCellFraction_Success = zeros(0,1);
MeanNNDCV_Success = zeros(0,1); MeanRDFFirstPeak_Success = zeros(0,1);
MeanMinClearance_Success = zeros(0,1); StdMinClearance_Success = zeros(0,1);
MeanOverlappingPairs_Success = zeros(0,1);

for im=1:numel(expectedMethods)
    methodNow = expectedMethods(im);
    for it=1:nVF
        vf = vfList(it);
        isRow = (T.Method==methodNow) & abs(T.TargetVF-vf)<1e-12;
        nTotal = sum(isRow);
        if nTotal==0
            warning('No %s data at target %.4f; no point will be displayed.', ...
                char(methodNow),vf);
            continue;
        end
        isSuccess = isRow & T.ReachedTarget;
        nGood = sum(isSuccess);
        eligibility = (nGood == nTotal) && (nTotal == expectedNRuns); % enforce complete original dataset

        Method(end+1,1) = methodNow;
        TargetVF(end+1,1) = vf;
        NRuns(end+1,1) = nTotal;
        NSuccess(end+1,1) = nGood;
        SuccessRate_percent(end+1,1) = 100*nGood/nTotal;
        EligibleForFigs234(end+1,1) = eligibility;
        MeanFinalVF_All(end+1,1) = mean(T.FinalVF(isRow),'omitnan');
        StdFinalVF_All(end+1,1) = std(T.FinalVF(isRow),'omitnan');
        MeanRuntime_All_s(end+1,1) = mean(T.ElapsedTime_s(isRow),'omitnan');
        MeanRuntime_Success_s(end+1,1) = mean(T.ElapsedTime_s(isSuccess),'omitnan');
        StdRuntime_Success_s(end+1,1) = std(T.ElapsedTime_s(isSuccess),'omitnan');
        MeanLocalVFRSD_Success(end+1,1) = mean(T.LocalVFRSD(isSuccess),'omitnan');
        StdLocalVFRSD_Success(end+1,1) = std(T.LocalVFRSD(isSuccess),'omitnan');
        MeanLowVFCellFraction_Success(end+1,1) = ...
            mean(T.LowVFCellFraction(isSuccess),'omitnan');
        StdLowVFCellFraction_Success(end+1,1) = ...
            std(T.LowVFCellFraction(isSuccess),'omitnan');
        MeanNNDCV_Success(end+1,1) = mean(T.NNDCV(isSuccess),'omitnan');
        MeanRDFFirstPeak_Success(end+1,1) = mean(T.RDFFirstPeak(isSuccess),'omitnan');
        MeanMinClearance_Success(end+1,1) = mean(T.MinClearance(isSuccess),'omitnan');
        StdMinClearance_Success(end+1,1) = std(T.MinClearance(isSuccess),'omitnan');
        MeanOverlappingPairs_Success(end+1,1) = ...
            mean(T.OverlappingPairs(isSuccess),'omitnan');
    end
end

S = table(Method, TargetVF, NRuns, NSuccess, SuccessRate_percent, ...
    EligibleForFigs234, MeanFinalVF_All, StdFinalVF_All, ...
    MeanRuntime_All_s, MeanRuntime_Success_s, StdRuntime_Success_s, ...
    MeanLocalVFRSD_Success, StdLocalVFRSD_Success, ...
    MeanLowVFCellFraction_Success, StdLowVFCellFraction_Success, ...
    MeanNNDCV_Success, MeanRDFFirstPeak_Success, ...
    MeanMinClearance_Success, StdMinClearance_Success, ...
    MeanOverlappingPairs_Success);

writetable(S,fullfile(resultsFolder,'R1_C5_plot_summary.csv'));
RSA = S(S.Method=="RSA-PBC", ...
    {'TargetVF','NRuns','NSuccess','SuccessRate_percent','EligibleForFigs234'});
writetable(RSA,fullfile(resultsFolder,'R1_C5_plot_eligibility.csv'));

disp('RSA-PBC plot-eligibility audit (100% of runs must reach target):');
disp(RSA);

if ~any(RSA.EligibleForFigs234)
    warning(['RSA-PBC did not successfully reach any target in all ', ...
        'realizations. Figures 2-4 will contain only VGPP points.']);
end

for k=1:height(RSA)
    if ~RSA.EligibleForFigs234(k)
        fprintf(['RSA-PBC target %.3f: EXCLUDE from geometric/statistical ', ...
            'and elapsed-time comparisons (successful %d/%d). ', ...
            'RETAIN achieved-VF and success-rate results in Fig. 1.\n'], ...
            RSA.TargetVF(k),RSA.NSuccess(k),RSA.NRuns(k));
    end
end

%% ----------------------- COMMON PLOT STYLES ---------------------------
methodStyle = {'-o','-s'};
methodColor = [0.10 0.33 0.60; 0.80 0.36 0.15];
methodNames = cellstr(expectedMethods);

%% Fig. 1 -- achieveability: intentionally includes FAILED RSA targets
f1=figure('Color','w','Visible',visibleFigures,'Position',[100 100 1120 435]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on; box on;
for im=1:numel(expectedMethods)
    part=S(S.Method==expectedMethods(im),:);
    errorbar(part.TargetVF,part.MeanFinalVF_All,part.StdFinalVF_All, ...
        methodStyle{im},'Color',methodColor(im,:),'LineWidth',1.7, ...
        'MarkerSize',6,'MarkerFaceColor',methodColor(im,:), ...
        'DisplayName',methodNames{im});
end
plot(vfList,vfList,'k--','LineWidth',1.1,'DisplayName','Target');
xlabel('Prescribed target area fraction');
ylabel('Achieved area fraction');
legend('Location','northwest'); grid off;
nexttile; hold on; box on;
for im=1:numel(expectedMethods)
    part=S(S.Method==expectedMethods(im),:);
    plot(part.TargetVF,part.SuccessRate_percent,methodStyle{im}, ...
        'Color',methodColor(im,:),'LineWidth',1.7, ...
        'MarkerSize',6,'MarkerFaceColor',methodColor(im,:), ...
        'DisplayName',methodNames{im});
end
xlabel('Prescribed target area fraction');
ylabel('Target-reaching realizations (%)');
ylim([0 105]); legend('Location','best'); grid off;
exportgraphics(f1,fullfile(picsFolder,'R1_C5_Fig1_reachability.png'), ...
    'Resolution',figureDPI);

%% Fig. 2 -- spatial statistics only where the prescribed target was met
f2=figure('Color','w','Visible',visibleFigures,'Position',[100 100 1120 435]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on; box on;
for im=1:numel(expectedMethods)
    part=S(S.Method==expectedMethods(im) & S.EligibleForFigs234,:);
    if isempty(part),continue;end
    errorbar(part.TargetVF,100*part.MeanLocalVFRSD_Success, ...
        100*part.StdLocalVFRSD_Success,methodStyle{im}, ...
        'Color',methodColor(im,:),'LineWidth',1.7, ...
        'MarkerSize',6,'MarkerFaceColor',methodColor(im,:), ...
        'DisplayName',methodNames{im});
end
xlabel('Prescribed target area fraction');
ylabel('Local area-fraction RSD (%)');
legend('Location','best'); grid off;
nexttile; hold on; box on;
for im=1:numel(expectedMethods)
    part=S(S.Method==expectedMethods(im) & S.EligibleForFigs234,:);
    if isempty(part),continue;end
    errorbar(part.TargetVF,100*part.MeanLowVFCellFraction_Success, ...
        100*part.StdLowVFCellFraction_Success,methodStyle{im}, ...
        'Color',methodColor(im,:),'LineWidth',1.7, ...
        'MarkerSize',6,'MarkerFaceColor',methodColor(im,:), ...
        'DisplayName',methodNames{im});
end
xlabel('Prescribed target area fraction');
ylabel('Low-local-area-fraction cells (%)');
legend('Location','best'); grid off;
exportgraphics(f2,fullfile(picsFolder,'R1_C5_Fig2_spatial_heterogeneity.png'), ...
    'Resolution',figureDPI);

%% Fig. 3 -- geometric validation at *achieved targets* only
f3=figure('Color','w','Visible',visibleFigures,'Position',[100 100 1120 435]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on; box on;
for im=1:numel(expectedMethods)
    part=S(S.Method==expectedMethods(im) & S.EligibleForFigs234,:);
    if isempty(part),continue;end
    errorbar(part.TargetVF,part.MeanMinClearance_Success, ...
        part.StdMinClearance_Success,methodStyle{im}, ...
        'Color',methodColor(im,:),'LineWidth',1.7, ...
        'MarkerSize',6,'MarkerFaceColor',methodColor(im,:), ...
        'DisplayName',methodNames{im});
end
xlabel('Prescribed target area fraction');
ylabel('Minimum periodic surface clearance');
legend('Location','best'); grid off;
nexttile; hold on; box on;
for im=1:numel(expectedMethods)
    part=S(S.Method==expectedMethods(im) & S.EligibleForFigs234,:);
    if isempty(part),continue;end
    plot(part.TargetVF,part.MeanOverlappingPairs_Success,methodStyle{im}, ...
        'Color',methodColor(im,:),'LineWidth',1.7, ...
        'MarkerSize',6,'MarkerFaceColor',methodColor(im,:), ...
        'DisplayName',methodNames{im});
end
xlabel('Prescribed target area fraction');
ylabel('Mean overlapping pairs beyond tolerance');
legend('Location','best'); grid off;
exportgraphics(f3,fullfile(picsFolder,'R1_C5_Fig3_periodic_overlap_audit.png'), ...
    'Resolution',figureDPI);

%% Fig. 4 -- VGPP ordering; RSA runtime only at fully achieved target fractions
f4=figure('Color','w','Visible',visibleFigures,'Position',[100 100 1120 435]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; box on;
VG=S(S.Method=="VGPP" & S.EligibleForFigs234,:);
if ~isempty(VG)
    yyaxis left;
    p1=plot(VG.TargetVF, VG.MeanNNDCV_Success, '-o', ...
        'Color',methodColor(1,:),'LineWidth',1.7,'MarkerSize',6);
    ylabel('Nearest-neighbor distance CV');
    yyaxis right;
    p2=plot(VG.TargetVF, VG.MeanRDFFirstPeak_Success, '-s', ...
        'Color',methodColor(2,:),'LineWidth',1.7,'MarkerSize',6);
    ylabel('First RDF peak');
    legend([p1 p2], {'VGPP: NND CV','VGPP: RDF peak'},'Location','best');
end
xlabel('Prescribed target area fraction'); grid off;

nexttile; hold on; box on;
maxUpper=0;
for im=1:numel(expectedMethods)
    part=S(S.Method==expectedMethods(im) & S.EligibleForFigs234,:);
    if isempty(part),continue;end
    errorbar(part.TargetVF,part.MeanRuntime_Success_s, ...
        part.StdRuntime_Success_s,methodStyle{im}, ...
        'Color',methodColor(im,:),'LineWidth',1.7, ...
        'MarkerSize',6,'MarkerFaceColor',methodColor(im,:), ...
        'DisplayName',methodNames{im});
    maxUpper=max(maxUpper,max(part.MeanRuntime_Success_s+ ...
        part.StdRuntime_Success_s));
end
xlabel('Prescribed target area fraction');
ylabel('Elapsed runtime for reached target (s)');
legend('Location','best'); grid off;
if maxUpper>0 && isfinite(maxUpper)
    ylim([0 1.07*maxUpper]);
else
    ylim([0 1]);
end
exportgraphics(f4,fullfile(picsFolder,'R1_C5_Fig4_ordering_runtime.png'), ...
    'Resolution',figureDPI);

%% --------------------------- FINAL AUDIT -------------------------------
% Ensure excluded RSA points cannot leak into Figs. 2-4.
rsaExclude=S.Method=="RSA-PBC" & ~S.EligibleForFigs234;
rsaInclude=S.Method=="RSA-PBC" & S.EligibleForFigs234;
fprintf('\n===== R1 C5 plot-eligibility summary =====\n');
fprintf('All prescribed target fractions shown in Fig. 1: %s\n', ...
    mat2str(vfList',3));
fprintf('RSA-PBC allowed target fractions in Figs. 2-4: %s\n', ...
    mat2str(S.TargetVF(rsaInclude)',3));
fprintf('RSA-PBC excluded target fractions in Figs. 2-4: %s\n', ...
    mat2str(S.TargetVF(rsaExclude)',3));
fprintf('Figures saved in: %s\n',picsFolder);
fprintf('Audit CSVs saved in: %s\n',resultsFolder);
fprintf('No changes made to the manuscript.\n');
