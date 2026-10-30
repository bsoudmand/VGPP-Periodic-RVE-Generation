%% R2_C39_IndependentSpatialAudit.m
% Independent post-processing of EXISTING particle-coordinate files.
% No VGPP generation, no mechanical homogenization, no modifications to originals.
% MATLAB R2024b. No special toolbox required.
clc; clear; close all;
root = '.';
sourceDir = fullfile(root,'revision_results','R2_C39');
coordDir = fullfile(sourceDir,'coordinates');
outDir = fullfile(sourceDir,'independent_audit');
figDir = fullfile(root,'pics','R2_C39','independent_audit');
if ~exist(outDir,'dir'), mkdir(outDir); end
if ~exist(figDir,'dir'), mkdir(figDir); end
assert(isfolder(coordDir),'Coordinates folder not found: %s',coordDir);
files = dir(fullfile(coordDir,'R2_C39_L*_rep*.csv'));
assert(~isempty(files),'No R2_C39 coordinate CSVs found.');
% Near-contact shell thickness as a fraction of mean particle radius.
thresholdFactors = [0.005 0.01 0.02 0.05 0.10 0.20];
rdfEdges = (0:0.05:1.1)'; rdfMid=(rdfEdges(1:end-1)+rdfEdges(2:end))/2;
fixedWindowSide=2.5; tolerance=1e-8;
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',12);
% Per-run output and curve arrays
runRows=table(); thresholdRows=table(); rdfRows=table();
for k=1:numel(files)
    name=files(k).name;
    tok=regexp(name,'R2_C39_L(\d+)_rep(\d+)\.csv','tokens','once');
    if isempty(tok), warning('Skipped unexpected filename %s',name); continue; end
    L=str2double(tok{1})/100; rep=str2double(tok{2});
    P=readtable(fullfile(files(k).folder,name));
    required={'x','y','radius'};
    assert(all(ismember(required,P.Properties.VariableNames)), 'Missing x/y/radius: %s',name);
    x=double(P.x(:)); y=double(P.y(:)); r=double(P.radius(:)); N=numel(x);
    assert(all(isfinite([x;y;r])) && all(r>0),'Invalid coordinates/radii in %s',name);
    % Independent minimum-image pair distances and per-particle NN distances.
    minNN=inf(N,1); pairD=[]; pairGap=[];
    pairD=zeros(N*(N-1)/2,1); pairGap=zeros(size(pairD)); offset=0;
    for i=1:N-1
        dx=x(i+1:end)-x(i); dy=y(i+1:end)-y(i);
        dx=dx-L*round(dx/L); dy=dy-L*round(dy/L);
        d=hypot(dx,dy); m=numel(d);
        pairD(offset+(1:m))=d; pairGap(offset+(1:m))=d-(r(i)+r(i+1:end)); offset=offset+m;
        minNN(i)=min(minNN(i),min(d));
        minNN(i+1:N)=min(minNN(i+1:N),d);
    end
    meanNN=mean(minNN); cvNN=std(minNN,0)/meanNN;
    badPairs=sum(pairGap < -tolerance); maxOverlap=max([0;-pairGap]);
    % Periodic, fixed physical window count descriptor; centers, not area fraction.
    nb=max(2,round(L/fixedWindowSide));
    ix=min(nb-1,floor(mod(x,L)/L*nb)); iy=min(nb-1,floor(mod(y,L)/L*nb));
    c=accumarray([ix+1 iy+1],1,[nb nb],@sum,0);
    windowCV=std(c(:),0)/mean(c(:));
    % RDF exact shell normalization for unique unordered pairs.
    shell=pi*(rdfEdges(2:end).^2-rdfEdges(1:end-1).^2);
    h=histcounts(pairD,rdfEdges)';
    rdf=h./(N*(N-1)/2*shell/L^2);
    hGrid=L/round(L/(10/249));  % same integer-grid rule used in the generator
    relativeGridSpacing=hGrid/mean(r);
    q=table(L,rep,N,mean(r),hGrid,relativeGridSpacing,meanNN,cvNN,windowCV,...
        min(pairGap),maxOverlap,badPairs, 'VariableNames',...
        {'L','Rep','N','Radius','GridSpacing','GridSpacingOverRadius','MeanNN','NearestCV',...
         'WindowCountCV','MinGap','MaxOverlap','ViolatingPairs'});
    runRows=[runRows;q]; %#ok<AGROW>
    for j=1:numel(thresholdFactors)
        factor=thresholdFactors(j); shellThickness=factor*mean(r);
        nNear=sum(pairGap<=shellThickness & pairGap>=-tolerance);
        meanCoord=2*nNear/N;
        qr=table(L,rep,N,factor,shellThickness,meanCoord,nNear,...
            'VariableNames',{'L','Rep','N','ThresholdOverRadius','ShellThickness','MeanCoord','NearPairs'});
        thresholdRows=[thresholdRows;qr]; %#ok<AGROW>
    end
    for j=1:numel(rdfMid)
        qr=table(L,rep,rdfMid(j),rdf(j),'VariableNames',{'L','Rep','Distance','RDF'});
        rdfRows=[rdfRows;qr]; %#ok<AGROW>
    end
    fprintf('%3d/%3d  L=%4.1f rep=%02d N=%4d NN-CV=%.5f maxOverlap=%.3g\n',...
        k,numel(files),L,rep,N,cvNN,maxOverlap);
end
writetable(runRows,fullfile(outDir,'R2_C39_audit_per_run.csv'));
writetable(thresholdRows,fullfile(outDir,'R2_C39_contact_threshold_per_run.csv'));
writetable(rdfRows,fullfile(outDir,'R2_C39_rdf_per_run.csv'));
% Summary across independent realizations (SEM describes between-run uncertainty).
Ls=unique(runRows.L); ns=numel(Ls); summary=table(); contactSummary=table(); rdfSummary=table();
for i=1:ns
    L=Ls(i); z=runRows(runRows.L==L,:);
    nr=height(z);
    a=table(L,z.N(1),nr,mean(z.MeanNN),std(z.MeanNN)/sqrt(nr),...
        mean(z.NearestCV),std(z.NearestCV)/sqrt(nr),...
        mean(z.WindowCountCV),std(z.WindowCountCV)/sqrt(nr),...
        max(z.MaxOverlap),sum(z.ViolatingPairs),...
        'VariableNames',{'L','N','Runs','MeanNN','SE_MeanNN','MeanNearestCV','SE_NearestCV',...
        'MeanWindowCV','SE_WindowCV','MaxOverlap','TotalViolatingPairs'});
    summary=[summary;a]; %#ok<AGROW>
    for j=1:numel(thresholdFactors)
        f=thresholdFactors(j); t=thresholdRows(thresholdRows.L==L & abs(thresholdRows.ThresholdOverRadius-f)<1e-12,:);
        q=table(L,z.N(1),f,mean(t.MeanCoord),std(t.MeanCoord)/sqrt(height(t)),...
            'VariableNames',{'L','N','ThresholdOverRadius','MeanCoord','SE_Coord'});
        contactSummary=[contactSummary;q]; %#ok<AGROW>
    end
    for j=1:numel(rdfMid)
        t=rdfRows(rdfRows.L==L & abs(rdfRows.Distance-rdfMid(j))<1e-10,:);
        q=table(L,z.N(1),rdfMid(j),mean(t.RDF),std(t.RDF)/sqrt(height(t)),...
            'VariableNames',{'L','N','Distance','MeanRDF','SE_RDF'});
        rdfSummary=[rdfSummary;q]; %#ok<AGROW>
    end
end
writetable(summary,fullfile(outDir,'R2_C39_audit_summary.csv'));
writetable(contactSummary,fullfile(outDir,'R2_C39_contact_threshold_summary.csv'));
writetable(rdfSummary,fullfile(outDir,'R2_C39_rdf_audit_summary.csv'));
% Compare original per-run stats without altering original files.
orig=fullfile(sourceDir,'R2_C39_raw.csv');
if isfile(orig)
    O=readtable(orig);
    diffs=table();
    for k=1:height(runRows)
        q=runRows(k,:);
        idx=find(abs(O.DomainL-q.L)<1e-9 & O.Realization==q.Rep,1);
        if isempty(idx), continue; end
        deltaNN=q.NearestCV-O.NearestCV(idx);
        % Original near-contact uses 0.02*mean radius.
        s=thresholdRows(thresholdRows.L==q.L & thresholdRows.Rep==q.Rep & abs(thresholdRows.ThresholdOverRadius-0.02)<1e-12,:);
        deltaCoord=s.MeanCoord-O.NearContactCoord(idx);
        diffs=[diffs;table(q.L,q.Rep,deltaNN,deltaCoord,...
            'VariableNames',{'L','Rep','DeltaNearestCV','DeltaContactCoord'})]; %#ok<AGROW>
    end
    writetable(diffs,fullfile(outDir,'R2_C39_original_recalculation_differences.csv'));
    fprintf('Max abs NN-CV difference from original: %.3g\n',max(abs(diffs.DeltaNearestCV)));
    fprintf('Max abs coordination difference: %.3g\n',max(abs(diffs.DeltaContactCoord)));
end
% Figure 1: threshold sensitivity for each L
f=figure('Color','w','Position',[100 100 1000 640]); hold on;
for i=1:ns
    s=contactSummary(contactSummary.L==Ls(i),:);
    errorbar(s.ThresholdOverRadius,s.MeanCoord,s.SE_Coord,'-o','LineWidth',1.5,...
        'DisplayName',sprintf('L = %.1f',Ls(i)));
end
xlabel('Near-contact shell thickness / particle radius');ylabel('Mean near-contact coordination');
legend('Location','best');grid on;box on;
exportgraphics(f,fullfile(figDir,'R2_C39_contact_threshold_sensitivity.png'),'Resolution',300);close(f);
% Figure 2: NN-CV and mean NN with SE
f=figure('Color','w','Position',[100 100 1050 490]);
tiledlayout(1,2,'Padding','compact');
nexttile;errorbar(summary.N,summary.MeanNearestCV,summary.SE_NearestCV,'-o','LineWidth',1.7);grid on;
xlabel('Particle count, N');ylabel('Nearest-neighbor CV');
nexttile;errorbar(summary.N,summary.MeanNN,summary.SE_MeanNN,'-o','LineWidth',1.7);grid on;
xlabel('Particle count, N');ylabel('Mean nearest-neighbor distance');
exportgraphics(f,fullfile(figDir,'R2_C39_independent_nearest_neighbor.png'),'Resolution',300);close(f);
% Figure 3: independently calculated RDF per domain
f=figure('Color','w','Position',[100 100 900 550]); hold on;
for i=1:ns
    s=rdfSummary(rdfSummary.L==Ls(i),:);
    plot(s.Distance,s.MeanRDF,'LineWidth',1.6,'DisplayName',sprintf('L = %.1f',Ls(i)));
end
xlabel('Periodic pair distance (model units)');ylabel('g(r)');legend('Location','best');
grid on;box on;
exportgraphics(f,fullfile(figDir,'R2_C39_independent_rdf.png'),'Resolution',300);close(f);
% Figure 4: grid-size audit. h/r nearly fixed across domains by design.
f=figure('Color','w','Position',[100 100 900 550]);
scatter(runRows.GridSpacingOverRadius,runRows.NearestCV,35,runRows.L,'filled');
xlabel('Grid spacing / particle radius');ylabel('Per-realization nearest-neighbor CV');
cb=colorbar;cb.Label.String='Domain side L';grid on;box on;
exportgraphics(f,fullfile(figDir,'R2_C39_grid_spacing_diagnostic.png'),'Resolution',300);close(f);
fprintf('\nAudited %d coordinate files. Summary at: %s\n',height(runRows),outDir);
disp(summary);disp(contactSummary);
% Interpretation: This is geometric/statistical verification only; NOT effective-property RVE convergence.
