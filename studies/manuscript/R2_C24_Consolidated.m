%% R2_C24_Consolidated_V2.m
% Scientific Reports | Reviewer 2, Comment 24 | MATLAB R2024b
% Archived-data validation ONLY. No packing simulations or toolboxes required.
% This script evaluates reproducibility and sensitivity of measured spatial
% descriptors. It DOES NOT enforce prescribed RDF, psi6, or local-RSD targets.
% Distinguishes ensemble evidence (saved per-run tables) from the independent
% audit of ONE saved particle configuration. Do not pool these evidence tiers.
clc; clear; close all;

%% 1. USER SETTINGS
masterSeed = 240024;
Lx = 10; Ly = 10; expectedRadius = 0.18; overlapTol = 1e-8;
nBootstrap = 2000;                  % descriptive percentile bootstrap
quadratDivisions = [5 6 8 10 12];  % sensitivity; 8 x 8 is manuscript setting
lowAreaThresholdFactor = 0.75;      % manuscript definition: 75% of mean
runCoordinateAudit = true;
strictGeometryAudit = false; % diagnostics NEVER imply a failed geometry passes
makePlots = true;
% Set this to the folder containing your archived results, if desired.
projectFolder = '.';
% Optional explicit absolute paths; leave empty for automatic discovery.
rawC1Path = '';
rawC2Path = '';
coordinatePath = '';      % Optional FULL-PRECISION CSV only
matFilePath = '';         % REQUIRED FOR MAT VALIDATION: supply path, or auto-discover below
% If multiple coordinate sources exist, set matFilePath explicitly.
% Never use the rounded R1_C4_best_valid_coordinates.csv by default.

%% 2. DIRECTORIES AND FILE DISCOVERY
scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder), scriptFolder = pwd; end
if ~isfolder(projectFolder)
    driveMirror = '.';
    if isfolder(driveMirror), projectFolder = driveMirror;
    else, projectFolder = scriptFolder; end
end
outDir = fullfile(projectFolder,'revision_results','R2_C24');
picDir = fullfile(projectFolder,'pics','R2_C24');
if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(picDir), mkdir(picDir); end
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',13);
rng(masterSeed,'twister');

sourceNames = {'R1_C1_raw.csv','R1_C2_raw.csv', ...
               'R2_C24_full_precision_coordinates.csv'};
manualPaths = {rawC1Path,rawC2Path,coordinatePath};
sourcePaths = strings(1,3);
searchFolders = {outDir,projectFolder,scriptFolder,pwd, ...
    fullfile(projectFolder,'revision_results'), ...
    fullfile(projectFolder,'revision_results','R1_C1'), ...
    fullfile(projectFolder,'revision_results','R1_C1_benchmark'), ...
    fullfile(projectFolder,'revision_results','R1_C2'), ...
    fullfile(projectFolder,'revision_results','R1_C4_focused')};
for k=1:3
    if ~isempty(manualPaths{k})
        assert(isfile(manualPaths{k}), 'Explicit input not found: %s',manualPaths{k});
        sourcePaths(k) = string(manualPaths{k});
    else
        for q=1:numel(searchFolders)
            trial=fullfile(searchFolders{q},sourceNames{k});
            if isfile(trial),sourcePaths(k)=string(trial);break;end
        end
        % Fallback: only look within the selected project's revision_results.
        if sourcePaths(k)=="" && isfolder(fullfile(projectFolder,'revision_results'))
            matches=dir(fullfile(projectFolder,'revision_results','**',sourceNames{k}));
            if numel(matches)==1
                sourcePaths(k)=string(fullfile(matches(1).folder,matches(1).name));
            elseif numel(matches)>1
                warning('Multiple files named %s found. Specify its full path.',sourceNames{k});
            end
        end
    end
end
%% PREFERRED COORDINATE SOURCE: ORIGINAL MAT, NEVER ROUNDED CSV
% All numerical coordinates below come from MATLAB's saved table object.
matCandidates = {matFilePath, ...
    fullfile(projectFolder,'revision_results','R1_C4_focused','R1_C4_best_valid_configuration.mat'), ...
    fullfile(projectFolder,'revision_results','R1_C4_focused','R1_C4_best_valid_configuration(1).mat'), ...
    fullfile(scriptFolder,'R1_C4_best_valid_configuration.mat'), ...
    fullfile(scriptFolder,'R1_C4_best_valid_configuration(1).mat')};
selectedMAT = '';
for k=1:numel(matCandidates)
    if ~isempty(matCandidates{k}) && isfile(matCandidates{k})
        selectedMAT = matCandidates{k}; break;
    end
end
if ~isempty(matFilePath) && ~isfile(matFilePath)
    error('Explicit MAT file not found: %s',matFilePath);
end
if ~isempty(selectedMAT)
    fprintf('Recovering original full-precision MATLAB coordinates: %s\n',selectedMAT);
    M = load(selectedMAT);
    assert(isfield(M,'bestParticleTable') && istable(M.bestParticleTable), ...
        'MAT file must contain table bestParticleTable.');
    V = M.bestParticleTable;
    names = lower(string(V.Properties.VariableNames));
    ix=find(ismember(names,["x","centerx","center_x"]),1);
    iy=find(ismember(names,["y","centery","center_y"]),1);
    ir=find(ismember(names,["radius","r","rad"]),1);
    assert(~isempty(ix)&&~isempty(iy)&&~isempty(ir), ...
        'MAT particle table must contain x, y and radius columns.');
    cx=double(V{:,ix});cy=double(V{:,iy});cr=double(V{:,ir});
    cx=cx(:);cy=cy(:);cr=cr(:);
    assert(all(isfinite([cx;cy;cr])) && all(cr>0) && ...
        all(cx>=0 & cx<Lx & cy>=0 & cy<Ly), ...
        'Original MAT coordinates failed basic validity checks.');
    exported=fullfile(outDir,'R2_C24_full_precision_coordinates.csv');
    fid=fopen(exported,'w');assert(fid>0,'Could not open full precision CSV.');
    fprintf(fid,'x,y,radius\n');
    for j=1:numel(cx)
        fprintf(fid,'%.17g,%.17g,%.17g\n',cx(j),cy(j),cr(j));
    end
    fclose(fid);
    sourcePaths(3)=string(exported);
    % Audit native doubles DIRECTLY, never reparse a rounded CSV.
    originalX = cx; originalY = cy; originalR = cr;
    % Source MAT is recorded to establish unambiguous provenance.
elseif sourcePaths(3)==""
    warning(['Original MAT coordinates and full-precision CSV missing. ', ...
        'The ensemble analyses will still run without geometry auditing.']);
end
assert(any(sourcePaths~=""), ...
    'No sources found; set rawC1Path, rawC2Path, or matFilePath.');

manifest=fopen(fullfile(outDir,'R2_C24_manifest.txt'),'w');
fprintf(manifest,'R2-C24 independent archived-data audit | %s\n',datestr(now));
fprintf(manifest,'MATLAB version: %s\n',version);
fprintf(manifest,'Lx=%g Ly=%g radius=%g tolerance=%g\n',Lx,Ly,expectedRadius,overlapTol);
fprintf(manifest,'Bootstrap=%d Seed=%d grids=%s\n',nBootstrap,masterSeed,mat2str(quadratDivisions));
fprintf(manifest,'Low-cell threshold = %g times mean local area fraction\n',lowAreaThresholdFactor);
for k=1:3,fprintf(manifest,'%s --> %s\n',sourceNames{k},char(sourcePaths(k)));end
fprintf(manifest,'The per-run metrics below are archived (not independently recomputed).\n');
fprintf(manifest,'Coordinate validation is single-configuration only.\n');
fprintf(manifest,'Original saved MAT source --> %s\n',selectedMAT);
fclose(manifest);

%% 3. ENSEMBLE-LEVEL ARCHIVED DESCRIPTORS: 20-RUN COMPARISON
if sourcePaths(1)~=""
    T=readtable(sourcePaths(1),'TextType','string');
    required={'Method','TargetVF','Run','ReachedTarget','LocalVFRSD'};
    assert(all(ismember(required,T.Properties.VariableNames)), ...
        'R1_C1_raw.csv does not contain expected columns.');
    T.Method=string(T.Method);
    T.TargetVF=str2double(string(T.TargetVF));
    T.Run=str2double(string(T.Run));
    T.LocalVFRSD=str2double(string(T.LocalVFRSD));
    reached=ismember(lower(strtrim(string(T.ReachedTarget))),["true","1","yes"]);
    if ismember('Completed',T.Properties.VariableNames)
        complete=ismember(lower(strtrim(string(T.Completed))),["true","1","yes"]);
    else
        complete=true(height(T),1);
    end
    valid=complete & reached;
    methods=unique(T.Method,'stable');fractions=unique(T.TargetVF(isfinite(T.TargetVF)));
    fprintf('Found %d archived C1 rows, %d methods.\n',height(T),numel(methods));
    % Warn about repeated method/target/run keys, which bias summaries.
    key=string(T.Method)+"|"+compose('%.8f',T.TargetVF)+"|"+string(T.Run);
    if numel(unique(key))~=numel(key)
        error('R1_C1_raw contains duplicate Method/TargetVF/Run IDs: inspect before analysis.');
    end
    S=table();
    for i=1:numel(methods)
        for j=1:numel(fractions)
            f=fractions(j);
            selected=T.Method==methods(i) & abs(T.TargetVF-f)<1e-10;
            if ~any(selected),continue;end
            good=selected & valid & isfinite(T.LocalVFRSD);
            values=T.LocalVFRSD(good);
            n=numel(values);meanRSD=NaN;sdRSD=NaN;ciLo=NaN;ciHi=NaN;
            if n>0,meanRSD=mean(values);sdRSD=std(values);end
            if n>=2
                bmeans=zeros(nBootstrap,1);
                for b=1:nBootstrap
                    bmeans(b)=mean(values(randi(n,n,1)));
                end
                bmeans=sort(bmeans);
                ciLo=bmeans(max(1,ceil(.025*nBootstrap)));
                ciHi=bmeans(min(nBootstrap,ceil(.975*nBootstrap)));
            end
            cvNN=NaN;
            if ismember('NearestNeighborCV',T.Properties.VariableNames)
                nn=str2double(string(T.NearestNeighborCV(good)));
                nn=nn(isfinite(nn));if ~isempty(nn),cvNN=mean(nn);end
            end
            row=table(methods(i),f,sum(selected),sum(selected & valid),n,...
                meanRSD,sdRSD,ciLo,ciHi,cvNN, ...
                'VariableNames',{'Method','NominalPhi','AttemptedRuns', ...
                'CompletedAtTarget','FiniteDescriptorRuns','MeanLocalRSD', ...
                'SDLocalRSD','BootstrapMeanLo','BootstrapMeanHi','MeanNearestCV'});
            S=[S;row]; %#ok<AGROW>
        end
    end
    writetable(S,fullfile(outDir,'R2_C24_ensemble_descriptor_summary.csv'));
    disp(S);
    % Same Run IDs denote matched labels, NOT necessarily statistically
    % paired random streams. Differences are descriptive, not paired tests.
    D=table();refName="Full_VGPP";
    if any(methods==refName)
        for i=1:numel(methods)
            if methods(i)==refName,continue;end
            for j=1:numel(fractions)
                f=fractions(j);d=[];
                ids=T.Run(T.Method==refName & abs(T.TargetVF-f)<1e-10);
                for k=1:numel(ids)
                    ix=find(T.Method==refName & abs(T.TargetVF-f)<1e-10 & T.Run==ids(k) & valid,1);
                    iy=find(T.Method==methods(i) & abs(T.TargetVF-f)<1e-10 & T.Run==ids(k) & valid,1);
                    if isempty(ix)||isempty(iy),continue;end
                    if isfinite(T.LocalVFRSD(ix)) && isfinite(T.LocalVFRSD(iy))
                        d(end+1,1)=T.LocalVFRSD(ix)-T.LocalVFRSD(iy); %#ok<SAGROW>
                    end
                end
                md=NaN;sd=NaN;
                if ~isempty(d),md=mean(d);sd=std(d);end
                row=table(methods(i),f,numel(d),md,sd, ...
                    'VariableNames',{'Alternative','NominalPhi','MatchedLabels', ...
                    'MeanFullMinusAlternativeRSD','SDOfLabelMatchedDifferences'});
                D=[D;row]; %#ok<AGROW>
            end
        end
    end
    writetable(D,fullfile(outDir,'R2_C24_label_matched_descriptive_differences.csv'));
    if makePlots && ~isempty(S)
        fig=figure('Color','w','Visible','off','Position',[50 50 950 550]);hold on;box on;
        colors=lines(numel(methods));
        for i=1:numel(methods)
            w=S.Method==methods(i) & S.FiniteDescriptorRuns>0;
            if any(w)
                errorbar(S.NominalPhi(w),100*S.MeanLocalRSD(w),100*S.SDLocalRSD(w), ...
                    '-o','Color',colors(i,:),'MarkerFaceColor',colors(i,:), ...
                    'LineWidth',1.3,'DisplayName',char(methods(i)));
            end
        end
        xlabel('Nominal particle area fraction');ylabel('Archived local area-fraction RSD (%)');
        legend('Location','best');
        title('Archived per-run descriptors; error bars are across-run SD');
        exportgraphics(fig,fullfile(picDir,'R2_C24_archived_RSD_vs_phi.png'),'Resolution',250);
        close(fig);
    end
else
    warning('C1 raw file missing; ensemble method comparison skipped.');
end

%% 4. ARCHIVED LOCAL SIXFOLD ORDER (SMALL SAMPLE, DESCRIPTIVE ONLY)
if sourcePaths(2)~=""
    U=readtable(sourcePaths(2),'TextType','string');
    required={'Method','NominalVF','Psi6Local','ValidAtTarget'};
    assert(all(ismember(required,U.Properties.VariableNames)), ...
        'R1_C2_raw.csv does not contain expected columns.');
    U.Method=string(U.Method);
    U.NominalVF=str2double(string(U.NominalVF));
    U.Psi6Local=str2double(string(U.Psi6Local));
    U.ValidAtTarget=ismember(lower(strtrim(string(U.ValidAtTarget))),["true","1","yes"]);
    m=unique(U.Method,'stable');f=unique(U.NominalVF(isfinite(U.NominalVF)));
    B=table();
    for i=1:numel(m)
        for j=1:numel(f)
            v=U.Method==m(i) & abs(U.NominalVF-f(j))<1e-10 & ...
                U.ValidAtTarget & isfinite(U.Psi6Local);
            a=U.Psi6Local(v);if isempty(a),continue;end
            row=table(m(i),f(j),numel(a),mean(a),std(a), ...
                'VariableNames',{'Method','NominalPhi','ValidConfigurations', ...
                'MeanLocalPsi6','SDLocalPsi6'});
            B=[B;row]; %#ok<AGROW>
        end
    end
    writetable(B,fullfile(outDir,'R2_C24_archived_psi6_summary.csv'));
    disp(B);
    if makePlots && ~isempty(B)
        fig=figure('Color','w','Visible','off','Position',[60 60 900 520]);hold on;box on;
        colors=lines(numel(m));
        for i=1:numel(m)
            v=B.Method==m(i);
            errorbar(B.NominalPhi(v),B.MeanLocalPsi6(v),B.SDLocalPsi6(v), ...
                '-o','LineWidth',1.3,'Color',colors(i,:), ...
                'MarkerFaceColor',colors(i,:),'DisplayName',char(m(i)));
        end
        xlabel('Nominal particle area fraction');ylabel('Archived mean local |\psi_6|');
        title('Local sixfold order: limited independent configurations');
        legend('Location','best');
        exportgraphics(fig,fullfile(picDir,'R2_C24_archived_psi6.png'),'Resolution',250);
        close(fig);
    end
else
    warning('C2 raw file missing; sixfold-order summary skipped.');
end

%% 5. INDEPENDENT PERIODIC GEOMETRY AND GRID-SCALE SENSITIVITY
if runCoordinateAudit && sourcePaths(3)~=""
    C=readtable(sourcePaths(3),'TextType','string');
    names=lower(string(C.Properties.VariableNames));
    ix=find(ismember(names,["x","centerx","center_x"]),1);
    iy=find(ismember(names,["y","centery","center_y"]),1);
    ir=find(ismember(names,["radius","r","rad"]),1);
    assert(~isempty(ix)&&~isempty(iy)&&~isempty(ir), ...
        'Coordinate CSV requires x, y and radius columns.');
    if exist('originalX','var') && exist('originalY','var') && exist('originalR','var')
        % Original MATLAB double arrays retained in memory. This removes
        % any ambiguity involving CSV parsing or serialization.
        x=originalX(:); y=originalY(:); r=originalR(:);
        auditOrigin='NATIVE MAT DOUBLE ARRAYS';
    else
        % Full precision CSV fallback. Require no assumption of validity.
        x=str2double(string(C{:,ix})); y=str2double(string(C{:,iy}));
        r=str2double(string(C{:,ir}));
        x=x(:); y=y(:); r=r(:);
        auditOrigin='CSV FALLBACK (CHECK SOURCE AND PRECISION)';
    end
    N=numel(x);
    fprintf('Geometry audit origin: %s\n',auditOrigin);
    assert(N>=2 && numel(y)==N && numel(r)==N,'Inconsistent geometry columns.');
    assert(all(isfinite([x;y;r])) && all(r>0) && ...
        all(x>=0 & x<Lx & y>=0 & y<Ly),'Invalid periodic disk geometry.');
    assert(max(abs(r-expectedRadius))<1e-8,'Unexpected particle radius.');
    phi=sum(pi*r.^2)/(Lx*Ly);
    fprintf('\nCOORDINATE INPUT USED: %s\n',char(sourcePaths(3)));
    fprintf('N=%d, phi=%.12f, min radius=%.12g, max radius=%.12g\n', ...
        N,phi,min(r),max(r));
    minGap=inf;badPairs=0;nearest=inf(N,1);boundaryPairs=0;
    badID1=zeros(0,1);badID2=zeros(0,1);badGaps=zeros(0,1);
    badDistances=zeros(0,1);badDx=zeros(0,1);badDy=zeros(0,1);
    for i=1:N-1
        dxRaw=abs(x(i)-x(i+1:N));dyRaw=abs(y(i)-y(i+1:N));
        crossX=dxRaw>Lx/2;crossY=dyRaw>Ly/2;
        dx=min(dxRaw,Lx-dxRaw);dy=min(dyRaw,Ly-dyRaw);
        dist=hypot(dx,dy);gap=dist-(r(i)+r(i+1:N));
        minGap=min(minGap,min(gap));
        offending=find(gap < -overlapTol);
        badPairs=badPairs+numel(offending);
        if ~isempty(offending)
            badID1=[badID1;repmat(i,numel(offending),1)]; %#ok<AGROW>
            badID2=[badID2;i+offending]; %#ok<AGROW>
            badGaps=[badGaps;gap(offending)]; %#ok<AGROW>
            badDistances=[badDistances;dist(offending)]; %#ok<AGROW>
            badDx=[badDx;dx(offending)]; %#ok<AGROW>
            badDy=[badDy;dy(offending)]; %#ok<AGROW>
        end
        boundaryPairs=boundaryPairs+sum((crossX|crossY)&dist <= (r(i)+r(i+1:N)+0.05));
        nearest(i)=min(nearest(i),min(dist));
        nearest(i+1:N)=min(nearest(i+1:N),dist);
    end
    nnCV=std(nearest)/mean(nearest);
    G=table(N,phi,minGap,badPairs,boundaryPairs,nnCV, ...
        'VariableNames',{'Particles','AreaFraction','MinPeriodicClearance', ...
        'PairsExceedingOverlapTol','NearBoundaryImagePairs','NearestNeighborCV'});
    writetable(G,fullfile(outDir,'R2_C24_single_geometry_audit.csv'));
    fprintf('Audit summary: minimum gap %.16g, bad pairs %d, tolerance %.3g\n',minGap,badPairs,overlapTol);
    writetable(G(:,{'Particles','AreaFraction','MinPeriodicClearance','PairsExceedingOverlapTol'}), ...
        fullfile(outDir,'R2_C24_MAT_independent_audit.csv'));
    disp(G);
    BadPairTable=table(badID1,badID2,badGaps,badDistances,badDx,badDy, ...
        'VariableNames',{'ParticleRow1','ParticleRow2','SurfaceGap', ...
        'PeriodicCenterDistance','PeriodicDx','PeriodicDy'});
    writetable(BadPairTable,fullfile(outDir,'R2_C24_overlap_violation_pairs.csv'));
    writetable(BadPairTable(:,{'ParticleRow1','ParticleRow2','SurfaceGap'}), ...
        fullfile(outDir,'R2_C24_MAT_overlap_pairs.csv'));
    fid=fopen(fullfile(outDir,'R2_C24_coordinate_source_and_geometry.txt'),'w');
    fprintf(fid,'Selected coordinate file: %s\n',char(sourcePaths(3)));
    fprintf(fid,'N=%d area_fraction=%.12f minimum_gap=%.12g violating_pairs=%d\n', ...
        N,phi,minGap,badPairs);
    fprintf(fid,'Geometry audit origin: %s\n',auditOrigin);
    fprintf(fid,'Geometry audit status: %s\n',char(string(badPairs==0)));
    fclose(fid);
    if badPairs>0
        warning(['Saved configuration has %d pair(s) with gap < -%.3g. ', ...
            'Input file: %s. Inspect R2_C24_overlap_violation_pairs.csv ', ...
            'and do not use this configuration as validated geometry.'], ...
            badPairs,overlapTol,char(sourcePaths(3)));
        disp(BadPairTable(1:min(10,height(BadPairTable)),:));
        if strictGeometryAudit
            error(['Geometry audit failed. Review the source pathname and ', ...
                'R2_C24_overlap_violation_pairs.csv before interpreting results.']);
        else
            warning(['DIAGNOSTIC MODE: detected geometry remains INVALID ', ...
                'at the requested tolerance. Statistics are computed for ', ...
                'diagnostic purposes only, NOT validated geometry.']);
        end
    else
        fprintf('Independent periodic geometry audit PASSED: 0 overlap violations.\n');
    end

    % Center-assigned disk area matches the published center-bin local-area
    % fraction convention. It is NOT a disk/cell intersection calculation.
    Q=table();
    for z=1:numel(quadratDivisions)
        K=quadratDivisions(z);
        assert(K>=2 && K==round(K),'Each grid size must be an integer >=2.');
        cx=min(K,max(1,floor(K*x/Lx)+1));
        cy=min(K,max(1,floor(K*y/Ly)+1));
        linearIndex=sub2ind([K,K],cx,cy);
        areaPerParticle=pi*r.^2;
        cellParticleArea=accumarray(linearIndex,areaPerParticle,[K*K,1],@sum,0);
        localPhi=cellParticleArea / ((Lx/K)*(Ly/K));
        % LocalPhi may be >1 for small cells: particle area is assigned in full
        % to its center cell; this is a coarse-grained descriptor, not an
        % exact spatial disk/cell area-intersection measurement.
        cellRSD=std(localPhi)/mean(localPhi);
        lowCells=mean(localPhi < lowAreaThresholdFactor*mean(localPhi));
        row=table(K,numel(localPhi),mean(localPhi),cellRSD,lowCells,...
            sum(localPhi>1),badPairs==0, ...
            'VariableNames',{'CellsPerSide','TotalCells','MeanLocalPhi', ...
            'LocalPhiRSD','FractionBelow75PctMean','CellsWithAssignedPhiAbove1', ...
            'GeometryValidated'});
        Q=[Q;row]; %#ok<AGROW>
    end
    writetable(Q,fullfile(outDir,'R2_C24_single_geometry_grid_sensitivity.csv'));
    disp(Q);
    if makePlots
        fig=figure('Color','w','Visible','off','Position',[60 60 850 520]);
        plot(Q.CellsPerSide,100*Q.LocalPhiRSD,'-o','LineWidth',1.5);
        xlabel('Quadrat grid cells per side');ylabel('Center-assigned local area-fraction RSD (%)');
        title(sprintf('Grid sensitivity for one VGPP configuration, N=%d, \\phi_A=%.4f',N,phi));
        exportgraphics(fig,fullfile(picDir,'R2_C24_single_geometry_grid_sensitivity.png'),'Resolution',250);
        close(fig);
    end
else
    warning('No coordinate file: independent geometry/grid-sensitivity audit skipped.');
end

%% 6. INTERPRETATION AND LIMITS
fid=fopen(fullfile(outDir,'R2_C24_interpretation.txt'),'w');
fprintf(fid,'Reviewer 2 Comment 24 | Evidence and limits\n');
fprintf(fid,'1. Descriptive archived statistics, not imposed target-statistic control.\n');
fprintf(fid,'2. No particle-packing simulations were run in this validation.\n');
fprintf(fid,'3. Bootstrap intervals summarize replicate uncertainty, not equivalence.\n');
fprintf(fid,'4. Same Run IDs allow label-matched differences, not causal pairing.\n');
fprintf(fid,'5. n=2 sixfold-order comparisons remain exploratory.\n');
fprintf(fid,'6. One coordinate file supports only one independent geometric audit.\n');
fprintf(fid,'7. Quadrat analysis assigns each full disk area to its center cell.\n');
fprintf(fid,'8. Quadrat RSD is grid-scale dependent; low threshold is 75%% mean.\n');
fprintf(fid,'9. The model does not prescribe RDF, psi6 or RSD target functions.\n');
fprintf(fid,'10. Suggested wording: quantified (not controlled) spatial statistics.\n');
fprintf(fid,'11. Do not change manuscript claims before checking outputs.\n');
fprintf(fid,'12. Any geometry audit failures are documented explicitly; they are NOT silently accepted.\n');
fclose(fid);
fprintf('\nR2-C24 archived-data validation outputs saved in:\n%s\n%s\n',outDir,picDir);
