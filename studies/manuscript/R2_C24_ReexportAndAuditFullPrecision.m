%% R2_C24_ReexportAndAuditFullPrecision.m
% MATLAB R2024b - recover the actual particle table stored inside MAT.
% Does not change coordinates or rerun any packing simulation.
clc;
clear;
close all;

matFile = fullfile( ...
    '.', ...
    'revision_results', ...
    'R1_C4_focused', ...
    'R1_C4_best_valid_configuration.mat');

assert(isfile(matFile), ...
    'MAT file not found. Check the specified path.');

S = load(matFile);
assert(isfield(S,'bestParticleTable'), 'Variable bestParticleTable missing.');
T = S.bestParticleTable;
assert(istable(T), 'bestParticleTable is not a MATLAB table.');
disp('Original table columns:'); disp(T.Properties.VariableNames(:));
varNames = lower(string(T.Properties.VariableNames));
ix=find(ismember(varNames,["x","centerx","center_x"]),1);
iy=find(ismember(varNames,["y","centery","center_y"]),1);
ir=find(ismember(varNames,["radius","r","rad"]),1);
assert(~isempty(ix)&&~isempty(iy)&&~isempty(ir), ...
    'Identify x, y, radius columns from the printed column list.');
x=double(T{:,ix}); y=double(T{:,iy}); radius=double(T{:,ir});
x=x(:); y=y(:); radius=radius(:); N=numel(x);
Lx=10; Ly=10; tol=1e-8;
assert(numel(y)==N&&numel(radius)==N&&all(isfinite([x;y;radius])) ...
    &&all(radius>0),'Invalid geometry.');
assert(all(x>=0&x<Lx&y>=0&y<Ly),'Centers outside base cell.');

% Independent double-precision periodic audit; full pair enumeration.
minGap=inf; nBad=0;
pairI=zeros(0,1);pairJ=zeros(0,1);pairGap=zeros(0,1);
for i=1:N-1
    j=(i+1:N)';
    dx=x(i)-x(j); dy=y(i)-y(j);
    dx=dx-Lx*round(dx/Lx);dy=dy-Ly*round(dy/Ly);
    d=hypot(dx,dy); gap=d-radius(i)-radius(j);
    minGap=min(minGap,min(gap));
    mask=gap < -tol;
    nBad=nBad+nnz(mask);
    pairI=[pairI; repmat(i,nnz(mask),1)]; %#ok<AGROW>
    pairJ=[pairJ;j(mask)]; %#ok<AGROW>
    pairGap=[pairGap;gap(mask)]; %#ok<AGROW>
end
phi=sum(pi*radius.^2)/(Lx*Ly);

outDir=fullfile('revision_results','R2_C24');
if ~isfolder(outDir),mkdir(outDir);end
outCSV=fullfile(outDir,'R2_C24_full_precision_coordinates.csv');
% fprintf avoids reliance on MATLAB's default table floating-point format.
fid=fopen(outCSV,'w');assert(fid>0,'Cannot open output CSV.');
fprintf(fid,'x,y,radius\n');
for i=1:N
    fprintf(fid,'%.17g,%.17g,%.17g\n',x(i),y(i),radius(i));
end
fclose(fid);
Audit=table(N,phi,minGap,nBad, ...
    'VariableNames',{'Particles','AreaFraction','MinPeriodicClearance','PairsBeyondTolerance'});
writetable(Audit,fullfile(outDir,'R2_C24_MAT_independent_audit.csv'));
writetable(table(pairI,pairJ,pairGap,'VariableNames', ...
    {'i','j','SurfaceGap'}),fullfile(outDir,'R2_C24_MAT_overlap_pairs.csv'));
disp(Audit);
if isfield(S,'bestRecord')
    disp('Original bestRecord summary:');disp(S.bestRecord);
end
if nBad==0
    fprintf('PASS: MAT coordinates meet tolerance %.1e.\n',tol);
else
    warning('MAT coordinates contain %d violating pairs; inspect CSV.',nBad);
end
fprintf('Saved full-precision coordinates: %s\n',outCSV);
