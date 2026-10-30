%% R2_C32_FullPrecisionCSV_Audit.m
% Reviewer 2, Comment 32: independent periodic geometry audit.
% Reads one original double-precision coordinate CSV with x,y,radius.
% Does not call the solver and does not use its overlap/pass flags.
% MATLAB R2024b. No external functions or toolboxes required.
clc; clear; close all;

%% Inputs
Lx = 10.0;
Ly = 10.0;
acceptanceTolerance = 1e-8; % model length units
inputFile = 'revision_results\R2_C32_coordinates\R2_C24_full_precision_coordinates.csv';
outputDir = 'revision_results\R2_C32_audit';
picsDir = 'pics';
if ~exist(outputDir,'dir'), mkdir(outputDir); end
if ~exist(picsDir,'dir'), mkdir(picsDir); end
assert(isfile(inputFile), 'Missing coordinate CSV: %s', inputFile);

%% Simple, explicit CSV reader; avoids table/group importer
fid = fopen(inputFile,'rt');
assert(fid ~= -1, 'Cannot open coordinate CSV: %s', inputFile);
firstLine = fgetl(fid);
assert(ischar(firstLine), 'The CSV is empty.');
if contains(lower(firstLine),'x') && contains(lower(firstLine),'y')
    vals = textscan(fid,'%f%f%f','Delimiter',',', ...
        'CollectOutput',true,'MultipleDelimsAsOne',false);
else
    frewind(fid);
    vals = textscan(fid,'%f%f%f','Delimiter',',', ...
        'CollectOutput',true,'MultipleDelimsAsOne',false);
end
fclose(fid);
coords = vals{1};
assert(size(coords,2)==3 && size(coords,1)>=2, ...
    'Expected coordinate data in three columns: x,y,radius.');
assert(all(isfinite(coords(:))) && all(coords(:,3)>0), ...
    'Invalid coordinates or particle radii.');

N = size(coords,1);
r = coords(:,3);
phi = sum(pi*r.^2)/(Lx*Ly);
numPairs = N*(N-1)/2;
minGap = inf;
maxOverlap = 0;
nStrict = 0;
nExceed = 0;
pairRows = zeros(0,5); % i,j,distancePBC,overlap,exceedsTolerance

%% Independent continuous-coordinate minimum-image pairwise audit
for i = 1:N-1
    j = (i+1:N)';
    dx = mod(abs(coords(i,1)-coords(j,1)),Lx);
    dy = mod(abs(coords(i,2)-coords(j,2)),Ly);
    dx = min(dx,Lx-dx);
    dy = min(dy,Ly-dy);
    distance = hypot(dx,dy);
    signedGap = distance - (r(i)+r(j));
    minGap = min(minGap,min(signedGap));
    overlap = max(0,-signedGap);
    maxOverlap = max(maxOverlap,max(overlap));
    nStrict = nStrict + nnz(overlap>0);
    nExceed = nExceed + nnz(overlap>acceptanceTolerance);
    ix = find(overlap>0);
    if ~isempty(ix)
        pairRows = [pairRows; ... %#ok<AGROW>
            [repmat(i,numel(ix),1),j(ix),distance(ix), ...
             overlap(ix),double(overlap(ix)>acceptanceTolerance)]];
    end
end

%% Results
status = "PASS";
if nExceed>0, status="FAIL"; end
Run = "FullPrecisionR1C4";
T = table(Run,N,phi,numPairs,minGap,maxOverlap,nStrict,nExceed,status, ...
    'VariableNames',{'Source','N','realizedAreaFraction','uniquePairsChecked', ...
    'minSignedSurfaceGap','maxGeometricOverlap','nStrictOverlaps', ...
    'nToleranceExceedances','ToleranceAudit'});
writetable(T,fullfile(outputDir,'R2_C32_run_audit.csv'));
S = T;
writetable(S,fullfile(outputDir,'R2_C32_target_summary.csv'));
P = array2table(pairRows,'VariableNames', ...
    {'i','j','distancePBC','overlapLength','exceedsTolerance'});
writetable(P,fullfile(outputDir,'R2_C32_positive_overlap_pairs.csv'));

fid = fopen(fullfile(outputDir,'R2_C32_audit_interpretation.txt'),'wt');
fprintf(fid,'Input: %s\nParticles: %d\nRealized area fraction: %.15g\n',inputFile,N,phi);
fprintf(fid,'Checked unique pairs: %d\nMinimum signed gap: %.16g\n',numPairs,minGap);
fprintf(fid,'Max positive overlap: %.16g\nPositive-overlap pairs: %d\n',maxOverlap,nStrict);
fprintf(fid,'Pairs exceeding %.1e: %d\nTolerance audit: %s\n',acceptanceTolerance,nExceed,status);
fclose(fid);

disp(T);
fprintf('Strictly positive overlap pairs: %d\n',nStrict);
fprintf('Overlap tolerance violations: %d\n',nExceed);
fprintf('Audit status: %s\n',status);
fprintf('Results saved in: %s\n',outputDir);

% For publication, retain the distinction between positive overlaps and
% tolerance-violating overlaps. One configuration does not validate 120 runs.
