%% R2_C32_IndependentOverlapAudit.m
% Independent, coordinate-based audit of periodic disk overlaps.
% MATLAB R2024b. Does NOT trust the solver's overlap/pass flags.
% Corrected V2: fixed overlap-pair table indexing; input processing failures are fatal.
% Requires final particle COORDINATES for each realization to be audited.
% Table_R1_fixedRadius_raw.csv contains aggregates, NOT these coordinates.
% Compatible inputs:
%   (a) separate CSV files with columns x, y, radius;
%   (b) a combined CSV with x, y, radius and Method/TargetVF/Run/Seed;
%   (c) MAT files containing particleTable, VGPPParticleTables (cell), or
%       a scalar struct with field particleTable.
% Put input files in revision_results/R2_C32_coordinates (or edit inputDir).
% File naming is not scientifically significant; metadata from names are
% optional and written as 'unknown' where not provided.
% Outputs: run audit, pair audit, summary, PNG histogram, and self-test.
clc; clear; close all;

%% User settings
Lx = 10.0;
Ly = 10.0;
minGap = 0.0;
finalRadius = 0.18;
overlapTolerance = 1e-8;              % same model length units as coordinates
strictThreshold = 0.0;               % mathematical positive-overlap threshold
nearContactThreshold = 0.02*finalRadius;
inputDir = fullfile('revision_results','R2_C32_coordinates');
outputDir = fullfile('revision_results','R2_C32_audit');
picsDir = 'pics';
if ~exist(outputDir,'dir'), mkdir(outputDir); end
if ~exist(picsDir,'dir'), mkdir(picsDir); end

%% Independent sanity checks (exact geometries, including periodic edge)
checks = [0.10 0.10 0.18; 0.47 0.10 0.18];
s = audit_coordinates(checks,Lx,Ly,0,overlapTolerance,nearContactThreshold);
assert(s.nStrictOverlaps == 0 && s.nToleranceExceedances == 0, ...
    'Self-test failed: non-overlapping pair.');
checks = [0.10 0.10 0.18; 0.45 0.10 0.18];
s = audit_coordinates(checks,Lx,Ly,0,overlapTolerance,nearContactThreshold);
assert(s.nStrictOverlaps == 1 && abs(s.maxGeometricOverlap-0.01)<1e-12, ...
    'Self-test failed: overlapping pair.');
checks = [0.10 0.10 0.18; 9.90 0.10 0.18];
s = audit_coordinates(checks,Lx,Ly,0,overlapTolerance,nearContactThreshold);
assert(s.nStrictOverlaps == 1 && abs(s.maxGeometricOverlap-0.16)<1e-12, ...
    'Self-test failed: periodic boundary overlap.');
checks = [0.10 0.10 0.18; 0.46-5e-9 0.10 0.18];
s = audit_coordinates(checks,Lx,Ly,0,overlapTolerance,nearContactThreshold);
assert(s.nStrictOverlaps == 1 && s.nToleranceExceedances == 0, ...
    'Self-test failed: sub-tolerance overlap.');
fprintf('Four independent geometric self-tests passed.\n');

%% Discover coordinate files
if ~exist(inputDir,'dir')
    mkdir(inputDir);
end
csvFiles = dir(fullfile(inputDir,'*.csv'));
matFiles = dir(fullfile(inputDir,'*.mat'));
allResults = table();
allPairs = table();
processed = 0;
processingFailures = strings(0,1);

for k = 1:numel(csvFiles)
    filePath = fullfile(csvFiles(k).folder,csvFiles(k).name);
    try
        T = readtable(filePath,'VariableNamingRule','preserve');
        [xcol,ycol,rcol] = coord_columns(T);
        if isempty(xcol)
            warning('Skipping %s: requires x, y, radius columns.',csvFiles(k).name);
            continue;
        end
        groups = group_indices(T);
        for g = 1:numel(groups)
            rows = groups{g};
            meta = get_metadata(T,rows,csvFiles(k).name);
            coords = [to_numeric(T.(xcol)(rows)), ...
                      to_numeric(T.(ycol)(rows)), ...
                      to_numeric(T.(rcol)(rows))];
            [allResults,allPairs] = record_case(allResults,allPairs,coords,...
                meta,Lx,Ly,minGap,overlapTolerance,nearContactThreshold);
            processed = processed+1;
        end
    catch ME
        warning('Could not process %s: %s',csvFiles(k).name,ME.message);
        processingFailures(end+1,1) = string(csvFiles(k).name) + ': ' + string(ME.message);
    end
end

for k = 1:numel(matFiles)
    filePath = fullfile(matFiles(k).folder,matFiles(k).name);
    try
        vars = load(filePath);
        names = fieldnames(vars);
        for v = 1:numel(names)
            a = vars.(names{v});
            if istable(a)
                [xc,yc,rc] = coord_columns(a);
                if ~isempty(xc)
                    coords = [to_numeric(a.(xc)),to_numeric(a.(yc)),to_numeric(a.(rc))];
                    meta = get_metadata(a,1:height(a),[matFiles(k).name ':' names{v}]);
                    [allResults,allPairs] = record_case(allResults,allPairs,coords,...
                        meta,Lx,Ly,minGap,overlapTolerance,nearContactThreshold);
                    processed=processed+1;
                end
            elseif iscell(a) && strcmpi(names{v},'VGPPParticleTables')
                for it=1:size(a,1)
                    for ir=1:size(a,2)
                        if ~istable(a{it,ir}) || isempty(a{it,ir}), continue; end
                        t=a{it,ir}; [xc,yc,rc]=coord_columns(t);
                        if isempty(xc), continue; end
                        coords=[to_numeric(t.(xc)),to_numeric(t.(yc)),to_numeric(t.(rc))];
                        meta=struct('Source',string(matFiles(k).name),... 
                            'Method',"VGPP",'Target',string(it),'Run',string(ir),'Seed',"unknown");
                        [allResults,allPairs] = record_case(allResults,allPairs,coords,...
                            meta,Lx,Ly,minGap,overlapTolerance,nearContactThreshold);
                        processed=processed+1;
                    end
                end
            elseif isstruct(a) && isscalar(a) && isfield(a,'particleTable') && istable(a.particleTable)
                t=a.particleTable; [xc,yc,rc]=coord_columns(t);
                if isempty(xc),continue; end
                coords=[to_numeric(t.(xc)),to_numeric(t.(yc)),to_numeric(t.(rc))];
                meta=get_metadata(t,1:height(t),[matFiles(k).name ':' names{v}]);
                [allResults,allPairs] = record_case(allResults,allPairs,coords,...
                    meta,Lx,Ly,minGap,overlapTolerance,nearContactThreshold);
                processed=processed+1;
            end
        end
    catch ME
        warning('Could not process %s: %s',matFiles(k).name,ME.message);
        processingFailures(end+1,1) = string(matFiles(k).name) + ': ' + string(ME.message);
    end
end

if processed == 0
    fprintf(2,['\nNO FINAL CONFIGURATION COORDINATES FOUND.\n' ...
        'Place per-run coordinate CSV/MAT files in:\n  %s\n' ...
        'Each CSV must include x,y,radius (and optionally Method, TargetVF, Run, Seed).\n' ...
        'The raw run summary Table_R1_fixedRadius_raw.csv ALONE cannot independently verify overlap.\n' ...
        'No manuscript-level geometric result can be concluded yet.\n'],inputDir);
    fid=fopen(fullfile(outputDir,'README_missing_coordinates.txt'),'w');
    fprintf(fid,'No coordinate files found. Export final particle coordinates for each run.\n');
    fclose(fid);
else
    writetable(allResults,fullfile(outputDir,'R2_C32_run_audit.csv'));
    if isempty(allPairs)
        allPairs=table(strings(0,1),strings(0,1),strings(0,1),strings(0,1),...
            zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),false(0,1),...
            'VariableNames',{'Source','Method','Target','Run','i','j',...
            'distancePBC','overlapLength','exceedsTolerance'});
    end
    writetable(allPairs,fullfile(outputDir,'R2_C32_positive_overlap_pairs.csv'));
    summary = groupsummary(allResults,{'Method','Target'},'mean',...
        {'maxGeometricOverlap','nStrictOverlaps','nToleranceExceedances','minSignedSurfaceGap'});
    writetable(summary,fullfile(outputDir,'R2_C32_target_summary.csv'));
    fig=figure('Color','w','Position',[80 80 900 550]);
    histogram(allResults.maxGeometricOverlap,20);
    xlabel('Maximum geometric overlap (model-length units)','FontName','Times New Roman');
    ylabel('Number of audited configurations','FontName','Times New Roman');
    title('Independent PBC coordinate audit','FontName','Times New Roman');
    set(gca,'FontName','Times New Roman','FontSize',13);
    exportgraphics(fig,fullfile(picsDir,'R2_C32_max_overlap_distribution.png'),'Resolution',300);
    fprintf('\nAudited %d final configurations.\n',height(allResults));
    fprintf('Configurations with strictly positive overlap: %d\n',sum(allResults.nStrictOverlaps>0));
    fprintf('Configurations exceeding adopted tolerance: %d\n',sum(allResults.nToleranceExceedances>0));
    fprintf('Maximum measured overlap: %.12g\n',max(allResults.maxGeometricOverlap));
    fprintf('NOTE: strict overlap and tolerance acceptance are different outcomes.\n');
end

if ~isempty(processingFailures)
    writelines(processingFailures,fullfile(outputDir,'R2_C32_processing_failures.txt'));
    error('R2-C32 audit incomplete: %d input file(s) could not be processed. See R2_C32_processing_failures.txt.',numel(processingFailures));
end

%% Local helpers -----------------------------------------------------------
function [xc,yc,rc] = coord_columns(T)
    cols=string(T.Properties.VariableNames);
    xc=pick(cols,["x","cx","centerx","center_x"]);
    yc=pick(cols,["y","cy","centery","center_y"]);
    rc=pick(cols,["radius","r","particle_radius"]);
end
function name=pick(cols,choices)
    name='';
    key=lower(regexprep(cols,'[^a-zA-Z0-9]',''));
    opts=lower(regexprep(choices,'[^a-zA-Z0-9]',''));
    for j=1:numel(opts)
        ix=find(key==opts(j),1);
        if ~isempty(ix),name=char(cols(ix));return;end
    end
end
function z=to_numeric(x)
    if isnumeric(x),z=double(x(:));else,z=str2double(string(x(:)));end
end
function groups=group_indices(T)
    cols=string(T.Properties.VariableNames);
    keys={};
    for wanted=["Method","TargetVF","TargetAreaFraction","Run","Seed"]
        ix=find(strcmpi(cols,wanted),1);
        if ~isempty(ix),keys{end+1}=char(cols(ix));end %#ok<AGROW>
    end
    if isempty(keys),groups={ (1:height(T))' };return;end
    [~,~,id]=unique(T(:,keys),'rows','stable');
    groups=accumarray(id,(1:height(T))',[],@(v){v});
end
function meta=get_metadata(T,rows,file)
    meta=struct('Source',string(file),'Method',"unknown",'Target',"unknown",...
        'Run',"unknown",'Seed',"unknown");
    keys={'Method',{'Method'};'Target',{'TargetVF','TargetAreaFraction','Target'};...
          'Run',{'Run'};'Seed',{'Seed'}};
    for j=1:size(keys,1)
        for k=1:numel(keys{j,2})
            col=pick(string(T.Properties.VariableNames),string(keys{j,2}{k}));
            if ~isempty(col),meta.(keys{j,1})=string(T.(col)(rows(1)));break;end
        end
    end
end
function [runs,pairs]=record_case(runs,pairs,C,meta,Lx,Ly,minGap,tol,nearContact)
    s=audit_coordinates(C,Lx,Ly,minGap,tol,nearContact);
    run=table(meta.Source,meta.Method,meta.Target,meta.Run,meta.Seed,...
        s.N,s.realizedAreaFraction,s.maxGeometricOverlap,...
        s.nStrictOverlaps,s.nToleranceExceedances,s.minSignedSurfaceGap,...
        s.nNearContacts,s.worstI,s.worstJ,...
        'VariableNames',{'Source','Method','Target','Run','Seed','N',...
        'realizedAreaFraction','maxGeometricOverlap','nStrictOverlaps',...
        'nToleranceExceedances','minSignedSurfaceGap','nNearContacts','worstI','worstJ'});
    runs=[runs;run]; %#ok<AGROW>
    if ~isempty(s.positivePairs)
        p=s.positivePairs;
        n=size(p,1);
        ptable=table(repmat(meta.Source,n,1),repmat(meta.Method,n,1),...
            repmat(meta.Target,n,1),repmat(meta.Run,n,1),...
            p(:,1),p(:,2),p(:,3),p(:,4),p(:,4)>tol,...
            'VariableNames',{'Source','Method','Target','Run','i','j',...
            'distancePBC','overlapLength','exceedsTolerance'});
        pairs=[pairs;ptable]; %#ok<AGROW>
    end
end
function s=audit_coordinates(C,Lx,Ly,minGap,tol,nearContact)
    assert(size(C,2)==3 && all(isfinite(C(:))) && all(C(:,3)>0),...
        'Coordinates must be finite Nx3 [x,y,radius] with positive radii.');
    n=size(C,1); s.N=n;
    s.realizedAreaFraction=sum(pi*C(:,3).^2)/(Lx*Ly);
    s.maxGeometricOverlap=0; s.nStrictOverlaps=0;
    s.nToleranceExceedances=0; s.minSignedSurfaceGap=NaN;
    s.nNearContacts=0; s.worstI=NaN; s.worstJ=NaN;
    s.positivePairs=zeros(0,4);
    if n<2,return;end
    minGapFound=inf;
    for i=1:n-1
        dx=abs(C(i,1)-C(i+1:end,1));
        dy=abs(C(i,2)-C(i+1:end,2));
        % mod handles centers outside the nominal cell robustly.
        dx=mod(dx,Lx); dx=min(dx,Lx-dx);
        dy=mod(dy,Ly); dy=min(dy,Ly-dy);
        d=hypot(dx,dy);
        gap=d-(C(i,3)+C(i+1:end,3));
        minGapFound=min(minGapFound,min(gap));
        overlap=max(0,-gap);
        s.nStrictOverlaps=s.nStrictOverlaps+nnz(overlap>0);
        s.nToleranceExceedances=s.nToleranceExceedances+nnz(overlap>tol);
        s.nNearContacts=s.nNearContacts+nnz(gap>=-tol & gap<=nearContact);
        [mx,idx]=max(overlap);
        if mx>s.maxGeometricOverlap
            s.maxGeometricOverlap=mx;
            s.worstI=i;s.worstJ=i+idx;
        end
        ids=find(overlap>0);
        if ~isempty(ids)
            s.positivePairs=[s.positivePairs; [repmat(i,numel(ids),1),i+ids,d(ids),overlap(ids)]]; %#ok<AGROW>
        end
    end
    s.minSignedSurfaceGap=minGapFound;
    % minGap is not part of strict physical overlap, but may be a clearance
    % requirement: report violations of prescribed gap independently.
    if minGap>0 && minGapFound<minGap
        fprintf('Warning: requested minGap %.5g not satisfied (minimum %.5g).\n',minGap,minGapFound);
    end
end
