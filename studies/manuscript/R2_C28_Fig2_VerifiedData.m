%% R2_C28_Fig2_VerifiedData.m
% Scientific Reports, VGPP revision: standalone reproduction of revised Fig. 2
% Two horizontal panels, matching the original Fig. 2 layout.
% Input: Table_R2_fixedRadius_summary.csv (existing 20-realization sweep).
% No particle generation and no modification of the input CSV.
% Only successful-target RSA-PBC values are displayed (0.40 and 0.50).
% Numerical RSD and low-cell fractions in CSV are FRACTIONS, converted to %.
% MATLAB R2024b; base MATLAB only.
clc; clear; close all;

%% Settings
inputCSV = 'revision_results\Table_R2_fixedRadius_summary.csv';  % Optional exact path. Empty = auto-discover or select file.
outputName = 'Fig2_spatial_heterogeneity_R2_C28_verified.png';
figureDPI = 600;
visibleFigures = 'on';  % Set 'off' for batch export.

%% Locate the original summary CSV (no data regeneration)
scriptDir = fileparts(mfilename('fullpath'));
searchCandidates = { ...
    fullfile(pwd,'revision_results','Table_R2_fixedRadius_summary.csv'), ...
    fullfile(pwd,'revision_results','Table_R2_fixedRadius_summary(1).csv'), ...
    fullfile(pwd,'Table_R2_fixedRadius_summary.csv'), ...
    fullfile(pwd,'Table_R2_fixedRadius_summary(1).csv'), ...
    fullfile(scriptDir,'revision_results','Table_R2_fixedRadius_summary.csv'), ...
    fullfile(scriptDir,'Table_R2_fixedRadius_summary.csv'), ...
    fullfile(scriptDir,'Table_R2_fixedRadius_summary(1).csv')};
if isempty(inputCSV)
    for k = 1:numel(searchCandidates)
        if isfile(searchCandidates{k})
            inputCSV = searchCandidates{k};
            break;
        end
    end
end
if isempty(inputCSV) || ~isfile(inputCSV)
    [f,p] = uigetfile('*.csv','Select the ORIGINAL Table_R2_fixedRadius_summary.csv');
    if isequal(f,0), error('No CSV selected; figure not generated.'); end
    inputCSV = fullfile(p,f);
end
fprintf('Fig. 2 data source: %s\n', inputCSV);
T = readtable(inputCSV,'TextType','string','VariableNamingRule','preserve');
required = {'Method','TargetVF','SuccessRate_percent', ...
    'MeanLocalVFRSD','StdLocalVFRSD', ...
    'MeanLowVFCellFraction','StdLowVFCellFraction'};
missing = setdiff(required,T.Properties.VariableNames);
assert(isempty(missing),'Missing CSV columns: %s',strjoin(missing,', '));

methods = strtrim(upper(string(T.Method)));
methods(contains(methods,'VGPP')) = "VGPP";
methods(contains(methods,'RSA')) = "RSA-PBC";
phiTarget = [0.40 0.50 0.55 0.60 0.65 0.70]';
metrics = {'MeanLocalVFRSD','StdLocalVFRSD','MeanLowVFCellFraction','StdLowVFCellFraction'};
YVGPP = nan(6,4); YRSA = nan(6,4);
for j = 1:6
    a = abs(T.TargetVF-phiTarget(j)) < 1e-7;
    iv = find(a & methods=="VGPP");
    assert(numel(iv)==1 && T.SuccessRate_percent(iv)==100, ...
        'VGPP missing/duplicated/not fully successful at %.2f.',phiTarget(j));
    for q = 1:4, YVGPP(j,q)=100*T.(metrics{q})(iv); end
    ir = find(a & methods=="RSA-PBC");
    assert(numel(ir)==1,'RSA-PBC missing or duplicated at %.2f.',phiTarget(j));
    if T.SuccessRate_percent(ir)==100
        for q = 1:4, YRSA(j,q)=100*T.(metrics{q})(ir); end
    end
end
assert(all(isfinite(YVGPP),'all'),'Invalid VGPP summary values.');
assert(all(isfinite(YRSA(1:2,:)),'all') && all(isnan(YRSA(3:end,:)),'all'), ...
    'RSA eligibility mismatch. Only 0.40 and 0.50 should be plotted.');
assert(all(YVGPP(:,2)>=0) && all(YRSA(1:2,2)>=0) && ...
       all(YVGPP(:,4)>=0) && all(YRSA(1:2,4)>=0),'Negative SD detected.');

%% Match original figure: 1 x 2 panels, Times New Roman, blue/orange curves
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
fig = figure('Color','w','Visible',visibleFigures, ...
    'Units','pixels','Position',[70 70 1500 630]);
blue = [0 0.4470 0.7410]; orange = [0.8500 0.3250 0.0980];
ax1 = axes(fig,'Units','normalized','Position',[0.084 0.15 0.395 0.76]);
ax2 = axes(fig,'Units','normalized','Position',[0.572 0.15 0.395 0.76]);
axesAll = [ax1 ax2];
for j = 1:2
    ax = axesAll(j); hold(ax,'on'); box(ax,'on'); grid(ax,'off');
    set(ax,'FontName','Times New Roman','FontSize',18,'LineWidth',3, ...
        'TickDir','in','XMinorTick','off','YMinorTick','off');
    xlim(ax,[0.39 0.71]); xticks(ax,phiTarget);
    xlabel(ax,'Prescribed target area fraction','FontSize',22);
    if j==1
        vCol=1; sCol=2;
        ylabel(ax,'Local area-fraction RSD (%)','FontSize',22);
        ylim(ax,[8 23]);
    else
        vCol=3; sCol=4;
        ylabel(ax,'Low-local-area-fraction cells (%)','FontSize',22);
        ylim(ax,[0 13]);
    end
    % Plot only two RSA points rather than all attempted (censored) fractions.
    hV = errorbar(ax,phiTarget,YVGPP(:,vCol),YVGPP(:,sCol),'-o', ...
        'Color',blue,'LineWidth',3,'MarkerSize',6,'MarkerFaceColor',blue, ...
        'CapSize',6);
    hR = errorbar(ax,phiTarget(1:2),YRSA(1:2,vCol),YRSA(1:2,sCol),'-s', ...
        'Color',orange,'LineWidth',3,'MarkerSize',6,'MarkerFaceColor',orange, ...
        'CapSize',6);
    legend(ax,[hV hR],{'VGPP','RSA-PBC'},'Location','northeast', ...
        'FontName','Times New Roman','FontSize',22,'Box','on');
end

%% Save output only under pics (leave benchmark CSV and manuscript figures intact)
projectDir = scriptDir;
if isfolder(fullfile(pwd,'revision_results'))
    projectDir = pwd;
elseif isfolder(fullfile(scriptDir,'revision_results'))
    projectDir = scriptDir;
else
    parentInput = fileparts(inputCSV);
    if strcmpi(string(lastFolderName(parentInput)),"revision_results")
        projectDir = fileparts(parentInput);
    end
end
picDir = fullfile(projectDir,'pics');
if ~isfolder(picDir), mkdir(picDir); end
pngPath = fullfile(picDir,outputName);
exportgraphics(fig,pngPath,'Resolution',figureDPI,'BackgroundColor','white');
fprintf('Saved revised Fig. 2: %s\n',pngPath);
fprintf('At 0.40: VGPP RSD %.3f%%, RSA-PBC %.3f%%\n',YVGPP(1,1),YRSA(1,1));
fprintf('At 0.50: VGPP RSD %.3f%%, RSA-PBC %.3f%%\n',YVGPP(2,1),YRSA(2,1));

function name = lastFolderName(pathStr)
    [~,name] = fileparts(pathStr);
end
