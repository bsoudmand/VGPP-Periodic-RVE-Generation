%% R2_C27_ParticleCountNormalization.m
% Reviewer 2 Comment 27: particle-count correction of 8x8 quadrat RSD.
% MATLAB R2024b, standalone SCRIPT; uses archived per-run measurements.
% Does NOT generate configurations or claim equilibrium reference validity.
clc; clear; close all;

%% Settings
C = 64; Lx = 10; Ly = 10; radius = 0.18;
projectFolder = '.';
inputFile = 'revision_results\R1_C1_BenchmarkAnalysis\R1_C1_raw.csv'; % Enter full path to R1_C1_raw.csv if automatic search fails
if ~isfolder(projectFolder), projectFolder = pwd; end
outDir = fullfile(projectFolder,'revision_results','R2_C27');
picDir = fullfile(projectFolder,'pics','R2_C27');
if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(picDir), mkdir(picDir); end

%% Locate existing 20-run comparison data
if isempty(inputFile)
    choices = {fullfile(projectFolder,'revision_results','R1_C1_raw.csv'), ...
       fullfile(projectFolder,'revision_results','R1_C1','R1_C1_raw.csv'), ...
       fullfile(projectFolder,'R1_C1_raw.csv'), fullfile(pwd,'R1_C1_raw.csv')};
    for z = 1:numel(choices)
        if isfile(choices{z}), inputFile = choices{z}; break; end
    end
end
if isempty(inputFile) || ~isfile(inputFile)
    [fn,fp] = uigetfile('*.csv','Choose archived R1_C1_raw.csv');
    assert(~isequal(fn,0),'No archive selected.');
    inputFile = fullfile(fp,fn);
end
T = readtable(inputFile,'TextType','string','VariableNamingRule','preserve');
required = {'Method','TargetVF','LocalVFRSD'};
assert(all(ismember(required,T.Properties.VariableNames)), ...
    'The CSV must contain Method, TargetVF and LocalVFRSD.');

%% Normalize field types, filter valid completed runs
methods = string(T.Method);
phiNom = str2double(string(T.TargetVF));
rsd = str2double(string(T.LocalVFRSD));
% str2double(string(numericArray)) also works with numeric MATLAB arrays.
valid = isfinite(phiNom) & isfinite(rsd) & rsd >= 0;
if ismember('ReachedTarget',T.Properties.VariableNames)
    v = lower(strtrim(string(T.ReachedTarget)));
    valid = valid & ismember(v,["1","true"]);
end
if ismember('Completed',T.Properties.VariableNames)
    v = lower(strtrim(string(T.Completed)));
    valid = valid & ismember(v,["1","true"]);
end
assert(any(valid),'No complete valid archive entries available.');

%% A consistent count for the prescribed fixed-radius periodic target
% Integer particle count follows the nominal fraction. It is not inferred
% from a rounded achieved area fraction. Use explicit count if available.
N = round(phiNom*(Lx*Ly)/(pi*radius^2));
if ismember('NFinal',T.Properties.VariableNames)
    maybeN = str2double(string(T.NFinal));
    replace = isfinite(maybeN) & maybeN >= 1;
    N(replace) = maybeN(replace);
elseif ismember('ParticleCount',T.Properties.VariableNames)
    maybeN = str2double(string(T.ParticleCount));
    replace = isfinite(maybeN) & maybeN >= 1;
    N(replace) = maybeN(replace);
end
valid = valid & isfinite(N) & N>=1;

% RSD is a FRACTION in the stored code (e.g., 0.12), not 12.
% Audit prevents accidentally treating a percentage-valued input as fraction.
assert(max(rsd(valid))<2, ['RSD values appear to be in percent. ', ...
    'Convert to fractions before using this script.']);
normalizedRSD = rsd .* sqrt(N/C);
FanoIndex = (N/C) .* rsd.^2;
Nexpected = N;
Archive = table(methods(valid),phiNom(valid),Nexpected(valid),rsd(valid), ...
    normalizedRSD(valid),FanoIndex(valid), ...
    'VariableNames',{'Method','NominalPhi','ParticleCount','LocalRSD', ...
       'CountNormalizedRSD','FanoDispersionIndex'});
writetable(Archive,fullfile(outDir,'R2_C27_per_run_normalized.csv'));

%% Summaries preserve method x target structure and count replicates
meth = unique(Archive.Method,'stable'); vals = unique(Archive.NominalPhi,'sorted');
Smethod = strings(0,1); Sphi = zeros(0,1); Sn = zeros(0,1);
SN = zeros(0,1); SRSD = zeros(0,1); SRSDsd = zeros(0,1);
Sscaled = zeros(0,1); SFano = zeros(0,1); SFanosd = zeros(0,1);
for a = 1:numel(meth)
    for b = 1:numel(vals)
        select = Archive.Method==meth(a) & abs(Archive.NominalPhi-vals(b))<1e-10;
        if ~any(select), continue; end
        v = Archive.FanoDispersionIndex(select);
        u = Archive.LocalRSD(select);
        Smethod(end+1,1) = meth(a);
        Sphi(end+1,1) = vals(b);
        Sn(end+1,1) = numel(v);
        SN(end+1,1) = mean(Archive.ParticleCount(select));
        SRSD(end+1,1) = mean(u);
        SRSDsd(end+1,1) = std(u);
        Sscaled(end+1,1) = mean(Archive.CountNormalizedRSD(select));
        SFano(end+1,1) = mean(v);
        SFanosd(end+1,1) = std(v);
    end
end
Summary = table(Smethod,Sphi,Sn,SN,SRSD,SRSDsd,Sscaled,SFano,SFanosd, ...
    'VariableNames',{'Method','NominalPhi','ValidRuns','MeanParticleCount', ...
    'MeanRSD','SDRSD','MeanCountNormalizedRSD','MeanFano','SDFano'});
writetable(Summary,fullfile(outDir,'R2_C27_normalized_summary.csv'));
disp(Summary);

%% Figure: original RSD and count-adjusted dispersion in separate PNGs
set(groot,'defaultAxesFontName','Times New Roman');
for graph = 1:2
    fig = figure('Color','w','Visible','off','Position',[100 100 950 570]);
    hold on; box on; cmap = lines(numel(meth));
    for a = 1:numel(meth)
        use = Summary.Method==meth(a);
        phi = Summary.NominalPhi(use);
        if graph==1
            y = 100*Summary.MeanRSD(use); e=100*Summary.SDRSD(use);
        else
            y = Summary.MeanFano(use); e=Summary.SDFano(use);
        end
        errorbar(phi,y,e,'-o','LineWidth',1.5,'Color',cmap(a,:), ...
            'MarkerFaceColor',cmap(a,:),'DisplayName',char(meth(a)));
    end
    xlabel('Nominal particle area fraction');
    if graph==1
        ylabel('Local area-fraction RSD (%)');
        filename = 'R2_C27_original_RSD.png';
    else
        ylabel('Count-adjusted quadrat dispersion index');
        filename = 'R2_C27_count_adjusted_Fano.png';
    end
    legend('Location','best'); grid on;
    exportgraphics(fig,fullfile(picDir,filename),'Resolution',300);
    close(fig);
end

%% Transparent record of what can and cannot be inferred
fid=fopen(fullfile(outDir,'R2_C27_interpretation.txt'),'w');
fprintf(fid,'Source: %s\n',inputFile);
fprintf(fid,'8x8 cells; C=%d; domain=%g x %g; disk radius=%g\n',C,Lx,Ly,radius);
fprintf(fid,'Counts inferred by rounding nominal phi unless explicit NFinal/ParticleCount provided.\n');
fprintf(fid,'Fano=(N/C)*(local RSD)^2; scaled RSD=sqrt(Fano).\n');
fprintf(fid,'For independent uniformly placed centers with fixed N, expected Fano=1.\n');
fprintf(fid,'This reference has no hard-core disk exclusion, and is NOT an equilibrium hard-disk baseline.\n');
fprintf(fid,'Values from archived descriptors, not newly recomputed from all source coordinates.\n');
fprintf(fid,'This study does not prescribe or optimize spatial-statistical targets.\n');
fclose(fid);
fprintf('R2-C27 analysis complete: %s\n',outDir);
