%% R4-C10: Regenerate revised Fig. 4 from the ORIGINAL 20-run benchmark
% MATLAB R2024b. No packing simulations or toolboxes required.
% Input: revision_results\Table_R1_fixedRadius_raw.csv under project root.
% Panel (a): VGPP mean NN-distance CV and first RDF peak (mean +/- SD).
% Panel (b): mean runtime +/- SD from successful runs only.
% Reviewer R5-C13: blue left y-axis in (a); linear 0-2000 s axis in (b).
% The large 0.60 VGPP outlier is retained. No data are deleted or winsorized.
clc; clear; close all;
root = fileparts(mfilename('fullpath'));
if isempty(root), root = pwd; end
% Search known original-benchmark locations; never substitute the 40-run data.
projectRoot = '.';
inputCandidates = {fullfile(projectRoot,'revision_results','Table_R1_fixedRadius_raw.csv'), ...
    fullfile(projectRoot,'Table_R1_fixedRadius_raw.csv'), ...
    fullfile(projectRoot,'revision_results','R4_C10','Table_R1_fixedRadius_raw.csv'), ...
    fullfile(root,'Table_R1_fixedRadius_raw.csv'), ...
    fullfile(root,'Table_R1_fixedRadius_raw(1).csv')};
inputFile = '';
for ii=1:numel(inputCandidates)
    if isfile(inputCandidates{ii}), inputFile=inputCandidates{ii}; break; end
end
assert(~isempty(inputFile), ...
    'Original benchmark CSV not found at the project revision_results folder or fallback locations.');
T = readtable(inputFile,'TextType','string');
required = {'Method','TargetVF','ReachedTarget','ElapsedTime_s','NNDCV','RDFFirstPeak'};
assert(all(ismember(required,T.Properties.VariableNames)), 'Input CSV missing required columns.');
assert(all(isfinite(T.ElapsedTime_s) & T.ElapsedTime_s>0),'Runtime must be finite and positive.');
phis = [0.40 0.50 0.55 0.60 0.65 0.70];
method = strtrim(string(T.Method));
% Full datasets are kept; do not replace failed target attempts with zero-time values.
meansCV=nan(size(phis)); sdCV=meansCV; meansRDF=meansCV; sdRDF=meansCV;
meanVG=meansCV; sdVG=meansCV; medianVG=meansCV; maxVG=meansCV;
meanRSA=meansCV; sdRSA=meansCV; nVG=zeros(size(phis)); nRSA=zeros(size(phis));
for k=1:numel(phis)
    p=phis(k);
    v=T(method=="VGPP" & abs(T.TargetVF-p)<1e-9 & T.ReachedTarget==1,:);
    a=T(method=="RSA-PBC" & abs(T.TargetVF-p)<1e-9 & T.ReachedTarget==1,:);
    nVG(k)=height(v); nRSA(k)=height(a);
    assert(~isempty(v),'VGPP successful results missing at phi=%.2f',p);
    meansCV(k)=mean(v.NNDCV); sdCV(k)=std(v.NNDCV);
    meansRDF(k)=mean(v.RDFFirstPeak); sdRDF(k)=std(v.RDFFirstPeak);
    meanVG(k)=mean(v.ElapsedTime_s); sdVG(k)=std(v.ElapsedTime_s);
    medianVG(k)=median(v.ElapsedTime_s); maxVG(k)=max(v.ElapsedTime_s);
    if ~isempty(a)
        meanRSA(k)=mean(a.ElapsedTime_s); sdRSA(k)=std(a.ElapsedTime_s);
    end
end
assert(all(nVG==20),'Expected 20 successful VGPP runs at each phi.');
assert(isequal(nRSA,[20 20 0 0 0 0]), ...
    'RSA-PBC success eligibility differs from original benchmark: inspect input.');

picsDir=fullfile(projectRoot,'pics','R4_C10');
resultsDir=fullfile(projectRoot,'revision_results','R4_C10');
if ~exist(picsDir,'dir'),mkdir(picsDir);end
if ~exist(resultsDir,'dir'),mkdir(resultsDir);end
% Preserve exact calculated statistics for transparent manuscript checking.
S=table(phis',nVG',nRSA',meansCV',sdCV',meansRDF',sdRDF', ...
 meanVG',sdVG',medianVG',maxVG',meanRSA',sdRSA', ...
 'VariableNames',{'TargetAF','VGPP_SuccessCount','RSA_SuccessCount', ...
 'NNDCV_Mean','NNDCV_SD','RDFFirstPeak_Mean','RDFFirstPeak_SD', ...
 'VGPP_RuntimeMean_s','VGPP_RuntimeSD_s','VGPP_RuntimeMedian_s', ...
 'VGPP_RuntimeMaximum_s','RSA_RuntimeMean_s','RSA_RuntimeSD_s'});
writetable(S,fullfile(resultsDir,'R4_C10_Fig4_values.csv'));

blue=[0.08 0.37 0.65]; orange=[0.85 0.35 0.13];
fig=figure('Color','w','Position',[100 100 1450 610]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
% PANEL (a): twin y-axes and error bars. Retain published descriptors.
ax1=nexttile(1); hold(ax1,'on');
yyaxis(ax1,'left');
hCV=errorbar(phis,meansCV,sdCV,'-o','Color',blue,'MarkerFaceColor',blue, ...
 'LineWidth',2,'MarkerSize',5,'CapSize',7);
ylabel('Nearest-neighbor distance CV');
yyaxis(ax1,'right');
hRDF=errorbar(phis,meansRDF,sdRDF,'-s','Color',orange,'MarkerFaceColor',orange, ...
 'LineWidth',2,'MarkerSize',5,'CapSize',7);
ylabel('First RDF peak');
% Explicit axis colors after creating both yyaxis series.
ax1.YAxis(1).Color=blue;
ax1.YAxis(2).Color=orange;
xlabel('Prescribed particle area fraction, \phi_A');
xlim([0.38 0.72]); xticks(phis); grid off; box on;
legend([hCV hRDF],{'VGPP: NND CV','VGPP: first RDF peak'}, ...
 'Location','north','FontSize',16);
text(ax1,0.02,0.05,'(a)','Units','normalized','FontWeight','bold','FontSize',18);

% PANEL (b): linear 0-2000 s y-axis requested under R5-C13.
% Arithmetic means and sample SDs are preserved in S/CSV; only whiskers
% extending outside the displayed axis limits are visually truncated.
ax2=nexttile(2); hold(ax2,'on');
yTop=2000;
lowVG=max(0,meanVG-sdVG);
highVG=meanVG+sdVG;
lowVGdisplay=min(lowVG,yTop);
highVGdisplay=min(highVG,yTop);
% MATLAB errorbar requires positive sizes. Clip at bounds for display only.
hVG=errorbar(ax2,phis,meanVG,meanVG-lowVGdisplay, ...
    highVGdisplay-meanVG,'-o','Color',blue,'MarkerFaceColor',blue, ...
    'LineWidth',2,'MarkerSize',5,'CapSize',7);
idx=isfinite(meanRSA);
meanRA=meanRSA(idx); sdRA=sdRSA(idx);
hRSA=errorbar(ax2,phis(idx),meanRA,meanRA-max(0,meanRA-sdRA), ...
    min(yTop,meanRA+sdRA)-meanRA,'-s', ...
    'Color',orange,'MarkerFaceColor',orange,'LineWidth',2, ...
    'MarkerSize',5,'CapSize',7);
set(ax2,'YScale','linear');
xlim(ax2,[0.38 0.72]); ylim(ax2,[0 yTop]);
xticks(ax2,phis); yticks(ax2,0:500:yTop);
xlabel(ax2,'Prescribed particle area fraction, \phi_A');
ylabel(ax2,'Elapsed time (s)'); grid(ax2,'off'); box(ax2,'on');
legend(ax2,[hVG hRSA],{'VGPP (n=20)','RSA-PBC (successful targets only)'}, ...
    'Location','west','FontSize',16);
legend(ax2,'AutoUpdate','off'); % Prevent auxiliary truncation markers from entering legend.
text(ax2,0.02,0.95,'(b)','Units','normalized', ...
    'FontWeight','bold','FontSize',18,'VerticalAlignment','top');
% Top-truncation symbol corresponds to an SD whisker crossing 2000 s.
clipped=find(highVG>yTop);
for j=clipped
    plot(ax2,phis(j),yTop-20,'^','MarkerSize',7, ...
        'MarkerFaceColor',blue,'MarkerEdgeColor',blue,'Clipping','on', ...
        'HandleVisibility','off');
end
if ~isempty(clipped)
    text(ax2,0.02,0.86,'Triangle: upper SD exceeds 2000 s', ...
        'Units','normalized','FontSize',18,'VerticalAlignment','top');
end
for ax=[ax1 ax2]
 set(ax,'FontName','Times New Roman','FontSize',18,'LineWidth',2, ...
  'TickDir','in','XMinorTick','off');
end
% Preserve the requested blue/orange axis labels and tick marks.
set(findall(fig,'Type','text'),'FontName','Times New Roman');
exportgraphics(fig,fullfile(picsDir,'Fig4_ordering_runtime_R5_C13.png'),'Resolution',400);
close(fig);

% A separate panel (b) is not needed; the combined figure is publication-ready.

fprintf('\nOriginal per-realization data verified. All recorded runtimes are > 0.\n');
fprintf('RSA-PBC means: phi=0.40 %.6f s, phi=0.50 %.6f s.\n',meanRSA(1),meanRSA(2));
fprintf('VGPP phi=0.60 mean %.3f s; median %.3f s; max %.3f s.\n', ...
 meanVG(4),medianVG(4),maxVG(4));
fprintf('Saved: %s\n',fullfile(picsDir,'Fig4_ordering_runtime_R5_C13.png'));
disp(S);
