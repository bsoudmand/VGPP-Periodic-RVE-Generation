%% R6_C4_PeriodicBoundaryCloseup.m
% Reviewer 6, Comment 4: Geometric illustration of periodic continuity.
% MATLAB R2024b; standalone script with local helper functions.
% This script DOES NOT generate new configurations or alter source data.
clc; clear; close all;

%% USER SETTINGS
% Prefer the original 688-particle coordinates used for Fig. 5(c).
% Set to its full path to label the output as Fig. 5(d).
originalFig5cCSV = '';  % example: 'D:\Projects\...\VGPP_P070_S001_coords.csv'

% If the Fig. 5(c) realization is not available, use archived VALID VGPP
% coordinates at phi_A = 0.760353 as a PROVISIONAL supplementary close-up.
% This must NOT be labeled as the Fig. 5(c) configuration.
allowProvisionalFallback = true;

% Windows project folder. Falls back to the script folder when unavailable.
projectFolder = '.';
Lx = 10; Ly = 10; overlapTolerance = 1e-8;
showFigure = 'on';        % 'on' or 'off'
zoomHalfWidth = 1.00;    % physical length visible on either side of seam
zoomHalfHeight = 1.25;   % physical height/half panel
lineWidth = 0.55;

%% OUTPUT FOLDERS
scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder), scriptFolder = pwd; end
if ~isfolder(projectFolder), projectFolder = scriptFolder; end
picsFolder = fullfile(projectFolder,'pics');
resultsFolder = fullfile(projectFolder,'revision_results','R6_C4');
if ~isfolder(picsFolder), mkdir(picsFolder); end
if ~isfolder(resultsFolder), mkdir(resultsFolder); end

%% INPUT DISCOVERY (never substitute LS benchmark coordinates for VGPP)
if ~isempty(originalFig5cCSV)
    assert(isfile(originalFig5cCSV), 'The specified Fig. 5(c) coordinate CSV does not exist.');
    inputCSV = originalFig5cCSV;
    sourceLabel = 'Original Fig. 5(c) VGPP configuration';
    isProvisional = false;
else
    assert(allowProvisionalFallback, ['Original Fig. 5(c) file not provided. ', ...
        'Supply originalFig5cCSV, or explicitly allow provisional fallback.']);
    choices = {fullfile(projectFolder,'revision_results','R1_C4_best_valid_coordinates.csv'), ...
        fullfile(scriptFolder,'R1_C4_best_valid_coordinates.csv')};
    inputCSV = '';
    for j = 1:numel(choices)
        if isfile(choices{j}), inputCSV = choices{j}; break; end
    end
    assert(~isempty(inputCSV), ['VGPP coordinate file not found. ', ...
        'Place R1_C4_best_valid_coordinates.csv next to this script, ', ...
        'or set originalFig5cCSV to the true Fig. 5(c) data.']);
    sourceLabel = 'Provisional VGPP configuration (NOT Fig. 5c)';
    isProvisional = true;
end

T = readtable(inputCSV, 'VariableNamingRule','preserve');
names = lower(string(T.Properties.VariableNames));
ix = find(ismember(names, ["x","centerx","center_x"]), 1);
iy = find(ismember(names, ["y","centery","center_y"]), 1);
ir = find(ismember(names, ["radius","r","rad"]), 1);
assert(~isempty(ix) && ~isempty(iy) && ~isempty(ir), ...
    'CSV must contain x, y and radius (or r) columns.');
x = double(T{:,ix}); y = double(T{:,iy}); r = double(T{:,ir});
x=x(:); y=y(:); r=r(:); N=numel(x);
assert(N >= 2 && numel(y)==N && numel(r)==N, 'Invalid coordinate vector lengths.');
assert(all(isfinite([x;y;r])) && all(r>0), 'Non-finite or invalid values.');
assert(all(x>=0 & x<Lx & y>=0 & y<Ly), 'Particle centers outside periodic base cell.');
assert(max(r)-min(r) < 1e-10, 'This visualization expects monodisperse disks.');
phi = sum(pi*r.^2)/(Lx*Ly);
if ~isProvisional
    assert(N==688, 'Fig. 5(c) was N=688; check the supplied original input.');
end
fprintf('Input: %s\nN=%d; area fraction = %.9f; radius=%.8f\n',inputCSV,N,phi,r(1));

%% EXACT INDEPENDENT PAIRWISE MINIMUM-IMAGE AUDIT
% Explicitly check boundary-crossing pairs, not just overlaps in base cell.
minGap=inf; closest=[NaN NaN]; overlaps=0;
leftRightPairs=0; topBottomPairs=0; cornerPairs=0;
minCrossGap=inf; nearestCross=[NaN NaN];
rows=zeros(0,8); % i,j,dxRaw,dyRaw,dxPBC,dyPBC,gap,crossCode
for i=1:N-1
    dxRaw=x(i+1:N)-x(i); dyRaw=y(i+1:N)-y(i);
    dx=dxRaw-Lx*round(dxRaw/Lx);
    dy=dyRaw-Ly*round(dyRaw/Ly);
    gap=sqrt(dx.^2+dy.^2)-(r(i)+r(i+1:N));
    [g,jloc]=min(gap);
    if g<minGap, minGap=g; closest=[i, i+jloc]; end
    overlaps=overlaps+nnz(gap < -overlapTolerance);
    crossX=abs(dxRaw-dx)>1e-10;
    crossY=abs(dyRaw-dy)>1e-10;
    leftRightPairs=leftRightPairs+nnz(crossX);
    topBottomPairs=topBottomPairs+nnz(crossY);
    cornerPairs=cornerPairs+nnz(crossX & crossY);
    % Only record NEAR-CONTACT PBC pairs for readable audit CSV.
    isNear=(crossX | crossY) & (gap <= 0.15);
    inds=find(isNear);
    if ~isempty(inds)
        code=double(crossX(inds))+2*double(crossY(inds));
        rows=[rows; [repmat(i,numel(inds),1), i+inds, dxRaw(inds), dyRaw(inds), ...
            dx(inds),dy(inds),gap(inds),code]]; %#ok<AGROW>
    end
    atCross=find(crossX|crossY);
    if ~isempty(atCross)
        [cg,jc]=min(gap(atCross));
        if cg<minCrossGap, minCrossGap=cg; nearestCross=[i,i+atCross(jc)]; end
    end
end
assert(overlaps==0, 'FAIL: saved coordinates have overlaps beyond tolerance.');
fprintf('PASS: minimum periodic gap = %.10g; overlaps beyond tolerance = 0\n',minGap);
fprintf('Cross-periodic pairs: LR=%d, TB=%d, both=%d; closest cross gap=%.9g\n', ...
    leftRightPairs,topBottomPairs,cornerPairs,minCrossGap);

%% SELECT MEANINGFUL SEAM CENTERS
% Select a particle near each periodic face, so the zoom shows actual
% particles intersecting/touching the boundary and their exact replicas.
edgeWidth=max(0.6,2.1*r(1));
nearLR=find(x<edgeWidth | x>Lx-edgeWidth);
nearTB=find(y<edgeWidth | y>Ly-edgeWidth);
assert(~isempty(nearLR) && ~isempty(nearTB), 'No particles near periodic boundaries.');
[~,a]=min(min(x(nearLR),Lx-x(nearLR)));
[~,b]=min(min(y(nearTB),Ly-y(nearTB)));
yLR = y(nearLR(a));
xTB = x(nearTB(b));

% In each subplot, original disks are blue and copies outside the base
% cell are gray/orange. Dashed seam is the actual periodic boundary.
fig=figure('Visible',showFigure,'Color','w','Units','pixels','Position',[90 90 1260 560]);
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax1=nexttile(tl,1);
plotSeam(ax1,x,y,r,Lx,Ly,'LR',yLR,zoomHalfWidth,zoomHalfHeight,lineWidth);
title(ax1,'Left/right periodic seam (x = 0 \equiv L_x)','FontWeight','normal');
ax2=nexttile(tl,2);
plotSeam(ax2,x,y,r,Lx,Ly,'TB',xTB,zoomHalfWidth,zoomHalfHeight,lineWidth);
title(ax2,'Bottom/top periodic seam (y = 0 \equiv L_y)','FontWeight','normal');

if isProvisional
    label=sprintf('Provisional VGPP PBC close-up | N=%d, \\phi_A=%.6f | NOT the Fig. 5(c) realization',N,phi);
    baseName='R6_C4_PROVISIONAL_PBC_Closeup_phi076';
else
    label=sprintf('VGPP Fig. 5(c) periodic boundary close-up | N=%d, \\phi_A=%.6f',N,phi);
    baseName='R6_C4_Fig5d_PBC_Closeup_phi070';
end
sgtitle(tl,label,'Interpreter','tex','FontWeight','normal','FontSize',12);
set(findall(fig,'-property','FontName'),'FontName','Times New Roman');
set(findall(fig,'Type','axes'),'FontSize',11);
exportgraphics(fig, fullfile(picsFolder,[baseName,'.png']), 'Resolution',350);

%% SAVE AUDIT CSV AND HUMAN-READABLE RECORD
nearTable=array2table(rows,'VariableNames', ...
    {'i','j','dxRaw','dyRaw','dxPeriodic','dyPeriodic','surfaceGap','crossCode'});
writetable(nearTable,fullfile(resultsFolder,'R6_C4_boundary_near_pairs.csv'));
summary=table(N,phi,r(1),minGap,overlaps,leftRightPairs,topBottomPairs,cornerPairs, ...
    minCrossGap,nearestCross(1),nearestCross(2),isProvisional, ...
    'VariableNames',{'N','AreaFraction','Radius','MinPeriodicGap', ...
    'OverlapsBeyondTolerance','LRCrossPairs','TBCrossPairs', ...
    'CornerCrossPairs','MinCrossPeriodicGap','ClosestCrossID1', ...
    'ClosestCrossID2','IsProvisional'});
writetable(summary,fullfile(resultsFolder,'R6_C4_periodic_geometry_audit.csv'));
fid=fopen(fullfile(resultsFolder,'R6_C4_notes.txt'),'w');
fprintf(fid,'R6-C4 PBC boundary audit\nSource: %s\nStatus: %s\n',inputCSV,sourceLabel);
fprintf(fid,'N=%d; area fraction=%.9f; radius=%.8f; Lx=%.6f; Ly=%.6f\n',N,phi,r(1),Lx,Ly);
fprintf(fid,'Minimum-image global gap=%.12g; overlaps beyond %.1e: %d\n',minGap,overlapTolerance,overlaps);
fprintf(fid,'Boundary-crossing pair count (all distances, not just contacts): LR=%d; TB=%d; both=%d\n', ...
    leftRightPairs,topBottomPairs,cornerPairs);
fprintf(fid,'Nearest boundary-crossing pair: (%d,%d); clearance=%.12g\n', ...
    nearestCross(1),nearestCross(2),minCrossGap);
fprintf(fid,'Pair CSV restricted to PBC crossing pairs with gap <= 0.15 length units.\n');
fprintf(fid,'CrossCode: 1 = left/right; 2 = top/bottom; 3 = both.\n');
if isProvisional
    fprintf(fid,['IMPORTANT: provisional source has 747 disks at phi~0.76; ', ...
        'it is NOT the N=688 phi~0.70 RVE in manuscript Fig. 5(c).\n', ...
        'Do not call it Fig. 5(d) until matching Fig. 5(c) coordinates are supplied.\n']);
end
fclose(fid);
fprintf('PNG: %s\n',fullfile(picsFolder,[baseName,'.png']));
fprintf('Audit CSVs and notes: %s\n',resultsFolder);

%% LOCAL FUNCTIONS
function plotSeam(ax,x,y,r,Lx,Ly,which,mid,halfWidth,halfHeight,lineWidth)
    hold(ax,'on'); axis(ax,'equal'); box(ax,'on'); grid(ax,'off');
    blue=[0.48 0.70 0.88]; edge=[0.22 0.40 0.62];
    imgFill=[0.98 0.78 0.45]; imgEdge=[0.67 0.35 0.08];
    th=linspace(0,2*pi,73);
    if strcmp(which,'LR')
        xlim(ax,[-halfWidth halfWidth]); ylim(ax,[mid-halfHeight mid+halfHeight]);
        % Base RVE x in [0,Lx] is drawn to right of seam; points near x=Lx
        % are shifted by -Lx and rendered as periodic copies to the left.
        for i=1:numel(x)
            yDisplay=y(i)-Ly*round((y(i)-mid)/Ly);
            for shift=-1:0
                xx=x(i)+shift*Lx;
                if abs(xx)>halfWidth+r(i), continue; end
                if abs(yDisplay-mid)>halfHeight+r(i), continue; end
                isImage=(shift~=0);
                if isImage, fc=imgFill; ec=imgEdge; else, fc=blue; ec=edge; end
                patch(ax,xx+r(i)*cos(th),yDisplay+r(i)*sin(th),fc, ...
                    'EdgeColor',ec,'LineWidth',lineWidth,'FaceAlpha',0.65);
            end
        end
        xline(ax,0,'k--','LineWidth',1.3);
        xlabel(ax,'Periodic x-position relative to seam');ylabel(ax,'y-position');
    else
        xlim(ax,[mid-halfHeight mid+halfHeight]); ylim(ax,[-halfWidth halfWidth]);
        for i=1:numel(x)
            xDisplay=x(i)-Lx*round((x(i)-mid)/Lx);
            for shift=-1:0
                yy=y(i)+shift*Ly;
                if abs(yy)>halfWidth+r(i),continue;end
                if abs(xDisplay-mid)>halfHeight+r(i),continue;end
                isImage=(shift~=0);
                if isImage, fc=imgFill;ec=imgEdge; else, fc=blue;ec=edge; end
                patch(ax,xDisplay+r(i)*cos(th),yy+r(i)*sin(th),fc, ...
                    'EdgeColor',ec,'LineWidth',lineWidth,'FaceAlpha',0.65);
            end
        end
        yline(ax,0,'k--','LineWidth',1.3);
        xlabel(ax,'x-position');ylabel(ax,'Periodic y-position relative to seam');
    end
    ax.Layer='top';
end
