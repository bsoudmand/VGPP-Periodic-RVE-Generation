%% R5-C1: Periodically consistent conceptual RSA-PBC / VGPP schematic
% -------------------------------------------------------------------------
% Four-panel schematic: RSA-PBC is an independent baseline (a); VGPP is (b-d).
%
% Improvements in this version:
% - 14 pt bold subplot titles
% - 14 pt bold axis labels and tick labels
% - 14 pt bold explanation-box text
% - Increased vertical distance between subplot titles and axes
% - True minimum-image periodic free-space map and periodic candidate separation
% - Actual independent pairwise PBC audit of schematic final configuration
% - Explicit separation of baseline panel (a) from VGPP panels (b-d)
% -------------------------------------------------------------------------

clc; clear; close all;

%% ------------------------------------------------------------------------
% General settings
% -------------------------------------------------------------------------

rng(7);

Lx = 10;
Ly = 10;
r  = 0.55;

theta = linspace(0, 2*pi, 180);

fontName = 'Times New Roman';

axisFontSize  = 14;
titleFontSize = 14;
boxFontSize   = 14;
mainTitleSize = 17;

% This controls how far above the axes each panel title is placed
titleGapY = 1.10;

fig = figure('Color','w', ...
    'Position',[60 60 1900 1020], ...
    'Toolbar','none', ...
    'MenuBar','none');

%% ------------------------------------------------------------------------
% Manual layout
% -------------------------------------------------------------------------

boxW  = 0.165;
plotW = 0.225;

gapBoxPlot = 0.045;
gapGroups  = 0.070;

rowH = 0.310;

topY = 0.535;
botY = 0.095;

xBox1  = 0.035;
xPlot1 = xBox1 + boxW + gapBoxPlot;

xBox2  = xPlot1 + plotW + gapGroups;
xPlot2 = xBox2 + boxW + gapBoxPlot;

%% ------------------------------------------------------------------------
% Base particle configuration
% -------------------------------------------------------------------------

particles = [
    2.0 2.0 r
    4.2 2.3 r
    6.6 2.0 r
    8.2 3.2 r
    2.6 5.0 r
    5.2 5.1 r
    7.5 5.6 r
    3.8 7.8 r
    6.5 8.0 r
];

cx = particles(:,1);
cy = particles(:,2);
rr = particles(:,3);

%% ------------------------------------------------------------------------
% Distance-transform field for panels (b) and (c)
% -------------------------------------------------------------------------

Nx = 350;
Ny = 350;

x = linspace(0,Lx,Nx);
y = linspace(0,Ly,Ny);
[X,Y] = meshgrid(x,y);

D = inf(size(X));

for i = 1:numel(cx)
    dx = abs(X - cx(i));
    dy = abs(Y - cy(i));

    dx = min(dx, Lx - dx);
    dy = min(dy, Ly - dy);

    di = hypot(dx, dy) - rr(i);
    D = min(D, di);
end

admissible = D >= r;

climMin = -0.5;
climMax = 1.65;

%% ========================================================================
% PANEL (a): Explanation box
% ========================================================================

axBox1 = axes('Position',[xBox1 topY boxW rowH]);
draw_explanation_box(axBox1, ...
    'BASELINE: RSA-PBC', ...
    {'Independent baseline'}, ...
    {'Random candidates', ...
     'are generated over', ...
     'the full domain.', ...
     '', ...
     'Many are rejected', ...
     'because of overlap.'}, ...
     fontName, boxFontSize);

%% ========================================================================
% PANEL (a): Conventional RSA
% ========================================================================

ax1 = axes('Position',[xPlot1 topY plotW rowH]);
hold(ax1,'on'); box(ax1,'on'); axis(ax1,'equal');
xlim(ax1,[0 Lx]); ylim(ax1,[0 Ly]);

draw_particles(ax1,cx,cy,rr,theta,[0.75 0.75 0.75]);

trialPts = [
    2.2 2.4
    5.0 5.3
    6.8 2.4
    4.5 2.1
];

for i = 1:size(trialPts,1)

    plot(ax1,trialPts(i,1), trialPts(i,2), 'x', ...
        'Color',[0.90 0.00 0.00], ...
        'MarkerSize',13, ...
        'LineWidth',2.2);

    viscircles(ax1,trialPts(i,:), r, ...
        'Color',[0.85 0.10 0.10], ...
        'LineStyle','--', ...
        'LineWidth',0.9);
end

acceptedPt = [8.4 8.2];

% Confirm illustrated RSA trial outcomes by periodic minimum-image checks.
for q = 1:size(trialPts,1)
    dx = abs(cx-trialPts(q,1)); dy = abs(cy-trialPts(q,2));
    dd = hypot(min(dx,Lx-dx),min(dy,Ly-dy));
    assert(any(dd < rr+r), 'Illustrated rejected RSA trial does not overlap');
end
dx = abs(cx-acceptedPt(1)); dy = abs(cy-acceptedPt(2));
dd = hypot(min(dx,Lx-dx),min(dy,Ly-dy));
assert(all(dd >= rr+r), 'Illustrated accepted RSA trial overlaps existing disks');

plot(ax1,acceptedPt(1), acceptedPt(2), 'o', ...
    'MarkerSize',8, ...
    'MarkerFaceColor',[0.10 0.65 0.25], ...
    'MarkerEdgeColor','k', ...
    'LineWidth',1.0);

viscircles(ax1,acceptedPt, r, ...
    'Color',[0.10 0.55 0.10], ...
    'LineStyle','--', ...
    'LineWidth',0.9);

apply_panel_title(ax1, ...
    '(a) RSA-PBC baseline: trial-and-rejection', ...
    fontName, titleFontSize, titleGapY);

xlabel(ax1,'x','FontName',fontName,'FontSize',axisFontSize,'FontWeight','bold');
ylabel(ax1,'y','FontName',fontName,'FontSize',axisFontSize,'FontWeight','bold');

set(ax1,'FontName',fontName, ...
    'FontSize',axisFontSize, ...
    'FontWeight','bold', ...
    'LineWidth',1.0);

% Compact legend inside panel
plot(ax1,0.65,9.35,'o', ...
    'MarkerSize',7, ...
    'MarkerFaceColor',[0.75 0.75 0.75], ...
    'MarkerEdgeColor','k');

text(ax1,0.95,9.35,'existing', ...
    'FontName',fontName, ...
    'FontSize',10, ...
    'FontWeight','normal', ...
    'VerticalAlignment','middle');

plot(ax1,3.10,9.35,'x', ...
    'Color',[0.90 0.00 0.00], ...
    'MarkerSize',8, ...
    'LineWidth',1.6);

text(ax1,3.40,9.35,'rejected', ...
    'FontName',fontName, ...
    'FontSize',10, ...
    'FontWeight','normal', ...
    'VerticalAlignment','middle');

plot(ax1,5.65,9.35,'o', ...
    'MarkerSize',7, ...
    'MarkerFaceColor',[0.10 0.65 0.25], ...
    'MarkerEdgeColor','k');

text(ax1,5.95,9.35,'accepted', ...
    'FontName',fontName, ...
    'FontSize',10, ...
    'FontWeight','normal', ...
    'VerticalAlignment','middle');

%% ========================================================================
% PANEL (b): Explanation box
% ========================================================================

axBox2 = axes('Position',[xBox2 topY boxW rowH]);
draw_explanation_box(axBox2, ...
    'VGPP: Void mapping', ...
    {'Admissible void', 'detection'}, ...
    {'Distance transform', ...
     'maps the remaining', ...
     'free space.', ...
     '', ...
     'Potential centers are', ...
     'identified from D.'}, ...
     fontName, boxFontSize);

%% ========================================================================
% PANEL (b): Admissible void-space map
% ========================================================================

ax2 = axes('Position',[xPlot2 topY plotW rowH]);
hold(ax2,'on'); box(ax2,'on'); axis(ax2,'equal');
xlim(ax2,[0 Lx]); ylim(ax2,[0 Ly]);

imagesc(ax2,x,y,D);
set(ax2,'YDir','normal');
colormap(ax2,parula);
clim(ax2,[climMin climMax]);

contour(ax2,X,Y,double(admissible),[0.5 0.5], ...
    'LineColor','k', ...
    'LineWidth',1.7);

draw_particles(ax2,cx,cy,rr,theta,[0.92 0.92 0.92]);

apply_panel_title(ax2, ...
    '(b) VGPP: admissible void-space map', ...
    fontName, titleFontSize, titleGapY);

xlabel(ax2,'x','FontName',fontName,'FontSize',axisFontSize,'FontWeight','bold');
ylabel(ax2,'y','FontName',fontName,'FontSize',axisFontSize,'FontWeight','bold');

set(ax2,'FontName',fontName, ...
    'FontSize',axisFontSize, ...
    'FontWeight','bold', ...
    'LineWidth',1.0);

cb2 = colorbar(ax2);
cb2.Label.String = 'Free-space distance';
cb2.Label.FontName = fontName;
cb2.Label.FontSize = 12;
cb2.Label.FontWeight = 'bold';
cb2.FontName = fontName;
cb2.FontSize = 11;
cb2.FontWeight = 'bold';

%% ========================================================================
% PANEL (c): Explanation box
% ========================================================================

axBox3 = axes('Position',[xBox1 botY boxW rowH]);
draw_explanation_box(axBox3, ...
    'VGPP: Insertion', ...
    {'Candidate', 'selection'}, ...
    {'High-clearance', ...
     'voids are identified.', ...
     '', ...
     'Weighted random', ...
     'candidates preserve', ...
     'stochastic dispersion.'}, ...
     fontName, boxFontSize);

%% ========================================================================
% PANEL (c): Largest-void-guided candidate selection
% ========================================================================

ax3 = axes('Position',[xPlot1 botY plotW rowH]);
hold(ax3,'on'); box(ax3,'on'); axis(ax3,'equal');
xlim(ax3,[0 Lx]); ylim(ax3,[0 Ly]);

imagesc(ax3,x,y,D);
set(ax3,'YDir','normal');
colormap(ax3,parula);
clim(ax3,[climMin climMax]);

draw_particles(ax3,cx,cy,rr,theta,[0.92 0.92 0.92]);

Dmasked = D;
Dmasked(~admissible) = -inf;

nVoid = 6;
voidPts = zeros(nVoid,2);
Dtemp = Dmasked;

for k = 1:nVoid
    [maxD, idx] = max(Dtemp(:));
    assert(isfinite(maxD),'Insufficient admissible voids for illustration');
    [row, col] = ind2sub(size(Dtemp), idx);

    voidPts(k,:) = [X(row,col), Y(row,col)];

    % Minimum-image separation avoids selecting duplicate voids across
    % opposite periodic boundaries.
    ddx = abs(X-voidPts(k,1));
    ddy = abs(Y-voidPts(k,2));
    ddx = min(ddx,Lx-ddx);
    ddy = min(ddy,Ly-ddy);
    distToSelected = hypot(ddx,ddy);
    Dtemp(distToSelected < 1.5*r) = -inf;
end

for k = 1:nVoid
    plot(ax3,voidPts(k,1), voidPts(k,2), 'ko', ...
        'MarkerFaceColor',[1.00 0.85 0.10], ...
        'MarkerSize',8, ...
        'LineWidth',1.2);

    viscircles(ax3,voidPts(k,:), r, ...
        'Color',[0.05 0.05 0.05], ...
        'LineStyle','--', ...
        'LineWidth',0.8);
end

text(ax3,1,10,'yellow markers: selected admissible voids', ...
    'FontName',fontName, ...
    'FontSize',10, ...
    'FontWeight','normal', ...
    'BackgroundColor','w', ...
    'EdgeColor',[0.70 0.70 0.70]);

apply_panel_title(ax3, ...
    '(c) VGPP: void-guided candidate selection', ...
    fontName, titleFontSize, titleGapY);

xlabel(ax3,'x','FontName',fontName,'FontSize',axisFontSize,'FontWeight','bold');
ylabel(ax3,'y','FontName',fontName,'FontSize',axisFontSize,'FontWeight','bold');

set(ax3,'FontName',fontName, ...
    'FontSize',axisFontSize, ...
    'FontWeight','bold', ...
    'LineWidth',1.0);

%% ========================================================================
% PANEL (d): Explanation box
% ========================================================================

axBox4 = axes('Position',[xBox2 botY boxW rowH]);
draw_explanation_box(axBox4, ...
    'VGPP: Final output', ...
    {'Completed', 'VGPP packing'}, ...
    {'Illustrative disk', ...
     'configuration.', ...
     '', ...
     'Displayed geometry', ...
     'passes an independent', ...
     'periodic overlap audit.'}, ...
     fontName, boxFontSize);

%% ========================================================================
% PANEL (d): Final VGPP packing
% ========================================================================

ax4 = axes('Position',[xPlot2 botY plotW rowH]);
hold(ax4,'on'); box(ax4,'on'); axis(ax4,'equal');
xlim(ax4,[0 Lx]); ylim(ax4,[0 Ly]);

basePts = [
    1.1 1.2
    2.7 1.1
    4.3 1.3
    5.9 1.1
    7.5 1.4
    8.9 1.0

    1.8 2.7
    3.5 2.5
    5.2 2.8
    6.9 2.5
    8.4 2.8

    1.1 4.2
    2.8 4.3
    4.6 4.1
    6.2 4.4
    7.9 4.2

    1.9 5.8
    3.6 5.9
    5.3 5.7
    7.0 5.9
    8.7 5.8

    1.2 7.4
    2.9 7.6
    4.7 7.4
    6.4 7.7
    8.2 7.4

    1.9 8.9
    3.8 8.7
    5.7 8.9
    7.5 8.8
    9.1 8.9
];

rng(21);
jitterAmp = 0.04;  % small perturbation; full minimum-image audit below
basePts = basePts + jitterAmp*(2*rand(size(basePts))-1);

basePts(:,1) = min(max(basePts(:,1),0.75),9.25);
basePts(:,2) = min(max(basePts(:,2),0.75),9.25);

finalPts = [basePts, r*ones(size(basePts,1),1)];

cxF = finalPts(:,1);
cyF = finalPts(:,2);
rF  = finalPts(:,3);

draw_particles(ax4,cxF,cyF,rF,theta,[0.72 0.72 0.72]);

for i = 1:4
    plot(ax4,[0 Lx],[i*2 i*2], ...
        ':', ...
        'Color',[0.85 0.85 0.85], ...
        'LineWidth',0.8);
end

% Actual independent O(N^2) periodic pair audit for the schematic.
% This demonstrates the displayed circles satisfy the non-overlap criterion;
% it does not imply that these schematic coordinates were generated by VGPP.
minimumClearance = inf;
violatingPairs = 0;
for ii = 1:numel(rF)-1
    dx = abs(cxF(ii)-cxF(ii+1:end));
    dy = abs(cyF(ii)-cyF(ii+1:end));
    dx = min(dx,Lx-dx); dy = min(dy,Ly-dy);
    clearances = hypot(dx,dy) - (rF(ii)+rF(ii+1:end));
    minimumClearance = min(minimumClearance,min(clearances));
    violatingPairs = violatingPairs + nnz(clearances < -1e-8);
end
assert(violatingPairs == 0, ...
    'Panel (d) contains periodic overlaps beyond 1e-8; correct coordinates');
fprintf('R5-C1 panel (d): %d disks; minimum periodic clearance %.6g; violations %d\n', ...
    numel(rF), minimumClearance, violatingPairs);

apply_panel_title(ax4, ...
    '(d) Illustrative packing: PBC-audited geometry', ...
    fontName, titleFontSize, titleGapY);

xlabel(ax4,'x','FontName',fontName,'FontSize',axisFontSize,'FontWeight','bold');
ylabel(ax4,'y','FontName',fontName,'FontSize',axisFontSize,'FontWeight','bold');

set(ax4,'FontName',fontName, ...
    'FontSize',axisFontSize, ...
    'FontWeight','bold', ...
    'LineWidth',1.0);

%% ------------------------------------------------------------------------
% Main title
% ------------------------------------------------------------------------

annotation(fig,'textbox', ...
    [0.20 0.955 0.60 0.035], ...
    'String','Schematic: independent RSA-PBC baseline and VGPP procedure', ...
    'FontName',fontName, ...
    'FontSize',mainTitleSize, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'LineStyle','none');

% The four individual panel titles identify the independent RSA-PBC
% comparison (a) and the three conceptual VGPP stages (b)-(d).
% No left/right grouping banner is used because panel (c) is bottom-left.

%% ------------------------------------------------------------------------
% Export
% ------------------------------------------------------------------------

set(fig,'PaperPositionMode','auto');
set(fig,'InvertHardcopy','off');

% Fixed project output folders (MATLAB R2024b / Windows)
resultDir = 'revision_results\R5_C1';
figDir    = 'pics\R5_C1';
if ~exist(resultDir,'dir'), mkdir(resultDir); end
if ~exist(figDir,'dir'), mkdir(figDir); end

% Figure image
outputFile = fullfile(figDir,'Fig8_RSA_PBC_vs_VGPP_revised.png');
exportgraphics(fig,outputFile,'Resolution',600);
fprintf('Revised conceptual comparison saved to: %s\n',outputFile);

% Numerical audit of the displayed illustrative configuration
pbcAudit = table(numel(rF), minimumClearance, violatingPairs, 1e-8, ...
    'VariableNames',{'ParticleCount','MinimumPeriodicClearance', ...
                     'ToleranceViolatingPairs','OverlapTolerance'});
auditFile = fullfile(resultDir,'R5_C1_Fig8_PBC_audit.csv');
writetable(pbcAudit,auditFile);
fprintf('Periodic audit saved to: %s\n',auditFile);

%% ========================================================================
% Local functions
% ========================================================================

function draw_particles(ax,cx,cy,r,theta,faceColor)

    for i = 1:length(cx)
        xp = cx(i) + r(i)*cos(theta);
        yp = cy(i) + r(i)*sin(theta);

        fill(ax,xp,yp,faceColor, ...
            'EdgeColor','k', ...
            'LineWidth',0.9);
    end
end

function draw_explanation_box(ax,stepLabel,headingLines,bodyLines,fontName,boxFontSize)

    cla(ax);
    axis(ax,'off');
    xlim(ax,[0 1]);
    ylim(ax,[0 1]);

    rectangle(ax, ...
        'Position',[0.04 0.05 0.92 0.90], ...
        'FaceColor',[0.97 0.97 0.97], ...
        'EdgeColor',[0.35 0.35 0.35], ...
        'LineWidth',1.1);

    text(ax,0.50,0.83,stepLabel, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontName',fontName, ...
        'FontSize',boxFontSize, ...
        'FontWeight','bold');

    headingText = strjoin(headingLines,newline);

    text(ax,0.50,0.66,headingText, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontName',fontName, ...
        'FontSize',boxFontSize, ...
        'FontWeight','bold', ...
        'Interpreter','none');

    bodyText = strjoin(bodyLines,newline);

    text(ax,0.50,0.37,bodyText, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontName',fontName, ...
        'FontSize',boxFontSize, ...
        'FontWeight','bold', ...
        'Interpreter','none');
end

function apply_panel_title(ax, titleStr, fontName, titleFontSize, yPos)

    ht = title(ax, titleStr, ...
        'FontName', fontName, ...
        'FontSize', titleFontSize, ...
        'FontWeight', 'bold');

    ht.Units = 'normalized';
    pos = ht.Position;
    pos(2) = yPos;   % increase vertical distance from plot
    ht.Position = pos;
end