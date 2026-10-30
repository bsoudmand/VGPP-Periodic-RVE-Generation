clc;
clear;
close all;

%% ============================================================
%  Standalone script to generate Fig. 9 for Reviewer #1 Comment 2
%  Fig. 9(a): mean local sixfold order vs area fraction
%  Fig. 9(b): RDF at phi_A = 0.60
%  Fig. 9(c): RDF at phi_A = 0.65
%  Fig. 9(d): RDF at phi_A = 0.70
%
%  Required input files:
%     R1_C2_summary.csv
%     R1_C2_rdf_summary.csv
%
%  Recommended output:
%     4 square individual panels + 1 combined preview figure
%% ============================================================

%% ---------------- USER SETTINGS ----------------
summaryFile = fullfile('revision_results', 'R1_C2_efficient','R1_C2_summary.csv');
rdfFile     = fullfile('revision_results', 'R1_C2_efficient','R1_C2_rdf_summary.csv');

% If your files are in another folder, modify the lines above, e.g.:
% summaryFile = fullfile('results','R1_C2_summary.csv');
% rdfFile     = fullfile('results','R1_C2_rdf_summary.csv');

outDir = fullfile('pics', 'R1_C2_Fig');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

fontName = 'Times New Roman';
fontSize = 16;
axisLineWidth = 3;
curveLineWidth = 3;
markerSize = 8;

figSize = 800;  % square canvas in pixels for each individual panel

methodIDs = ["Full_VGPP","WeightedMask_Only","UniformDomain_Growth"];
methodLabels = ["Full VGPP","Weighted sampling","Uniform sampling"];

colors = [0.0000 0.4470 0.7410;   % blue
          0.8500 0.3250 0.0980;   % orange
          0.9290 0.6940 0.1250];  % yellow

targetVF = [0.60 0.65 0.70];

%% ---------------- CHECK FILES ----------------
if ~isfile(summaryFile)
    error('Cannot find summary file: %s', summaryFile);
end

if ~isfile(rdfFile)
    error('Cannot find RDF summary file: %s', rdfFile);
end

%% ---------------- READ TABLES ----------------
summaryTbl = readtable(summaryFile,'TextType','string');
rdfTbl     = readtable(rdfFile,'TextType','string');

% Force method column to string if needed
summaryTbl.Method = string(summaryTbl.Method);
rdfTbl.Method     = string(rdfTbl.Method);

%% ---------------- BASIC COLUMN CHECK ----------------
requiredSummaryCols = ["Method","NominalVF","Psi6LocalMean"];
requiredRdfCols     = ["Method","NominalVF","RadiusOverDiameter","MeanG"];

for k = 1:numel(requiredSummaryCols)
    if ~ismember(requiredSummaryCols(k), string(summaryTbl.Properties.VariableNames))
        error('Missing required column "%s" in %s', requiredSummaryCols(k), summaryFile);
    end
end

for k = 1:numel(requiredRdfCols)
    if ~ismember(requiredRdfCols(k), string(rdfTbl.Properties.VariableNames))
        error('Missing required column "%s" in %s', requiredRdfCols(k), rdfFile);
    end
end

% If SD column does not exist, create zeros
if ~ismember("Psi6LocalSD", string(summaryTbl.Properties.VariableNames))
    summaryTbl.Psi6LocalSD = zeros(height(summaryTbl),1);
end

%% ============================================================
%% Fig. 9(a): Mean local sixfold order vs area fraction
%% ============================================================
figA = figure('Color','w','Units','pixels',...
    'Position',[100 100 figSize figSize]);

axA = axes(figA);
hold(axA,'on');
box(axA,'on');

for i = 1:numel(methodIDs)

    sub = summaryTbl(summaryTbl.Method == methodIDs(i), :);
    if isempty(sub)
        continue;
    end

    sub = sortrows(sub,'NominalVF');

    x = sub.NominalVF;
    y = sub.Psi6LocalMean;
    e = sub.Psi6LocalSD;
    e(~isfinite(e)) = 0;

    errorbar(axA, x, y, e, '-o', ...
        'Color', colors(i,:), ...
        'LineWidth', curveLineWidth, ...
        'MarkerSize', markerSize, ...
        'MarkerFaceColor', 'w', ...
        'CapSize', 8, ...
        'DisplayName', methodLabels(i));
end

xlabel(axA,'Nominal inclusion area fraction','FontName',fontName);
ylabel(axA,'Mean local sixfold order, \langle|\psi_6|\rangle',...
    'Interpreter','tex','FontName',fontName);

xlim(axA,[0.595 0.705]);

% automatic y-limits with small margin
allY = summaryTbl.Psi6LocalMean(isfinite(summaryTbl.Psi6LocalMean));
allE = summaryTbl.Psi6LocalSD(isfinite(summaryTbl.Psi6LocalSD));
if isempty(allE), allE = 0; end
yMin = min(allY - allE);
yMax = max(allY + allE);
pad = 0.01;
ylim(axA,[floor((yMin-pad)*100)/100, ceil((yMax+pad)*100)/100]);

set(axA,...
    'FontName',fontName,...
    'FontSize',fontSize,...
    'LineWidth',axisLineWidth,...
    'TickDir','in',...
    'Box','on',...
    'PlotBoxAspectRatio',[1 1 1]);

legend(axA,'Location','northwest','Box','on',...
    'FontName',fontName,'FontSize',13);

axA.Units = 'normalized';
axA.Position = [0.17 0.17 0.75 0.75];

exportgraphics(figA, fullfile(outDir,'Fig9a_MeanLocalOrder.png'),'Resolution',300);

%% ============================================================
%% Fig. 9(b), 9(c), 9(d): RDF at 0.60, 0.65, 0.70
%% ============================================================
panelLetters = {'b','c','d'};

for iv = 1:numel(targetVF)

    vf = targetVF(iv);

    figR = figure('Color','w','Units','pixels',...
        'Position',[100 100 figSize figSize]);

    axR = axes(figR);
    hold(axR,'on');
    box(axR,'on');

    ymax = 0;

    for i = 1:numel(methodIDs)

        sub = rdfTbl(rdfTbl.Method == methodIDs(i) & abs(rdfTbl.NominalVF - vf) < 1e-12, :);
        if isempty(sub)
            continue;
        end

        sub = sortrows(sub,'RadiusOverDiameter');

        plot(axR, sub.RadiusOverDiameter, sub.MeanG, ...
            'LineWidth', curveLineWidth, ...
            'Color', colors(i,:), ...
            'DisplayName', methodLabels(i));

        ymax = max(ymax, max(sub.MeanG));
    end

    xline(axR,1,'k--','LineWidth',1.0,'HandleVisibility','off');

    xlabel(axR,'Normalized separation, r/d','FontName',fontName);
    ylabel(axR,'Radial distribution function, g(r)',...
        'Interpreter','tex','FontName',fontName);

    xlim(axR,[0.8 5.0]);
    ylim(axR,[0, 1.05*ymax]);

    set(axR,...
        'FontName',fontName,...
        'FontSize',fontSize,...
        'LineWidth',axisLineWidth,...
        'TickDir','in',...
        'Box','on',...
        'PlotBoxAspectRatio',[1 1 1]);

    legend(axR,'Location','north','Box','on',...
        'FontName',fontName,'FontSize',13);

    axR.Units = 'normalized';
    axR.Position = [0.17 0.17 0.75 0.75];

    outName = sprintf('Fig9%s_RDF_%02d.png', panelLetters{iv}, round(100*vf));
    exportgraphics(figR, fullfile(outDir,outName),'Resolution',300);
end

%% ============================================================
%% Combined 2x2 figure with reduced middle margins
%% ============================================================

figC = figure('Color','w', ...
    'Units','pixels', ...
    'Position',[50 50 1400 1050]);

%% Compact 2x2 layout with 2:1 plotting frames
% [left bottom width height]

positions = [
    0.150 0.625 0.365 0.260;   % (a) top-left
    0.550 0.625 0.365 0.260;   % (b) top-right
    0.150 0.175 0.365 0.260;   % (c) bottom-left
    0.550 0.175 0.365 0.260];  % (d) bottom-right

for ip = 1:4

    ax = axes('Parent',figC, ...
        'Units','normalized', ...
        'Position',positions(ip,:));

    hold(ax,'on');
    box(ax,'on');

    if ip == 1

        for im = 1:numel(methodIDs)

            sub = summaryTbl( ...
                summaryTbl.Method == methodIDs(im),:);

            if isempty(sub), continue; end

            sub = sortrows(sub,'NominalVF');

            e = sub.Psi6LocalSD;
            e(~isfinite(e)) = 0;

            errorbar(ax, ...
                sub.NominalVF, ...
                sub.Psi6LocalMean,e, ...
                '-o', ...
                'Color',colors(im,:), ...
                'LineWidth',1.6, ...
                'MarkerSize',6, ...
                'MarkerFaceColor','w', ...
                'CapSize',6, ...
                'DisplayName',char(methodLabels(im)));
        end

        xlim(ax,[0.595 0.705]);
        ylim(ax,[0.40 0.52]);

        xlabel(ax,'Nominal inclusion area fraction');
        ylabel(ax,'Mean local sixfold order, \langle|\psi_6|\rangle', ...
            'Interpreter','tex');

        title(ax,'(a)','FontWeight','bold');

    else

        vf = targetVF(ip-1);
        ymax = 0;

        for im = 1:numel(methodIDs)

            sub = rdfTbl( ...
                rdfTbl.Method == methodIDs(im) & ...
                abs(rdfTbl.NominalVF-vf)<1e-12,:);

            if isempty(sub),continue;end

            sub = sortrows(sub,'RadiusOverDiameter');

            plot(ax,sub.RadiusOverDiameter,sub.MeanG, ...
                'LineWidth',1.6, ...
                'Color',colors(im,:), ...
                'DisplayName',char(methodLabels(im)));

            ymax = max(ymax,max(sub.MeanG));
        end

        xline(ax,1,'k--', ...
            'HandleVisibility','off');

        xlim(ax,[0.8 5]);
        ylim(ax,[0 1.05*ymax]);

        xlabel(ax,'Normalized separation, r/d');
        ylabel(ax,'Radial distribution function, g(r)');

        title(ax,sprintf('(%s)  \\phi_A = %.2f', ...
            panelLetters{ip-1},vf), ...
            'Interpreter','tex', ...
            'FontWeight','bold');
    end

     set(ax, ...
        'FontName','Times New Roman', ...
        'FontSize',16, ...
        'LineWidth',2, ...
        'TickDir','in', ...
        'PlotBoxAspectRatio',[2 1 1]);
    
    legend(ax,'Location','north', ...
        'FontSize',16, ...
        'Box','on');

end

exportgraphics(figC, ...
    fullfile(outDir,'R1_C2_Fig_CombinedPreview.png'), ...
    'Resolution',600);

fprintf('Compact Figure 9 preview generated.\n');