%% R2_C28_RSA_RSD_StatisticalComparison.m
% Scientific Reports - Reviewer 2 Comment 28
% Archived-data-only comparison of local area-fraction RSD: VGPP vs RSA-PBC.
% Standalone MATLAB R2024b script; no toolboxes or custom functions.
% No packing simulations and no claims of equivalence based on p > 0.05.
clc; clear; close all;

%% Settings
projectRoot = '.';
if ~isfolder(projectRoot), projectRoot = pwd; end
inputFile = fullfile(projectRoot,'revision_results','Table_R1_fixedRadius_raw.csv');
if ~isfile(inputFile), inputFile = fullfile(pwd,'Table_R1_fixedRadius_raw.csv'); end
assert(isfile(inputFile), 'Cannot find Table_R1_fixedRadius_raw.csv. Set inputFile to its actual location.');
outputDir = fullfile(projectRoot,'revision_results','R2_C28');
picsDir = fullfile(projectRoot,'pics','R2_C28');
if ~isfolder(outputDir), mkdir(outputDir); end
if ~isfolder(picsDir), mkdir(picsDir); end
rng(280028,'twister');
nBootstrap = 10000; nPermutation = 20000; targets = [0.40 0.50];

%% Read original VGPP/RSA archive
T = readtable(inputFile,'TextType','string');
required = {'Method','TargetVF','ReachedTarget','LocalVFRSD'};
assert(all(ismember(required,T.Properties.VariableNames)), 'Required columns missing in original raw CSV.');
method = upper(strtrim(string(T.Method)));
method(contains(method,'VGPP')) = "VGPP";
method(contains(method,'RSA')) = "RSA-PBC";
phi = double(T.TargetVF); rsd = double(T.LocalVFRSD);
if islogical(T.ReachedTarget)
    reached = T.ReachedTarget;
elseif isnumeric(T.ReachedTarget)
    reached = isfinite(T.ReachedTarget) & T.ReachedTarget > 0.5;
else
    reached = ismember(lower(strtrim(string(T.ReachedTarget))), ["1","true","yes"]);
end
% An archived raw CSV may encode local RSD as a fraction or percentage.
validRSD = rsd(isfinite(rsd));
assert(~isempty(validRSD), 'No finite RSD values.');
if median(validRSD) > 1, rsd = rsd/100; unitInfo = 'Input was percentage; converted to fraction';
else, unitInfo = 'Input was fractional RSD'; end

%% Unpaired independent-realization comparison: Welch t + bootstrap + permutation
% Do not automatically pair by run number: seed labels alone do not prove paired design.
Result = table(); Raw = table();
for k=1:numel(targets)
    f=targets(k);
    a = rsd(method=="VGPP" & abs(phi-f)<1e-9 & reached & isfinite(rsd));
    b = rsd(method=="RSA-PBC" & abs(phi-f)<1e-9 & reached & isfinite(rsd));
    na=numel(a); nb=numel(b);
    assert(na==20 && nb==20, 'Expected 20 successful independent realizations for both methods at phi=%.2f; got %d and %d',f,na,nb);
    ma=mean(a); mb=mean(b); sa=std(a); sb=std(b); delta=ma-mb;
    se=sqrt(sa^2/na+sb^2/nb);
    df=(sa^2/na+sb^2/nb)^2 / ((sa^2/na)^2/(na-1)+(sb^2/nb)^2/(nb-1));
    % MATLAB base has no guaranteed Student t CDF / inverse. Use optional tcdf,
    % otherwise permutation p and bootstrap CI remain primary inference.
    tWelch=delta/se; pWelch=NaN;
    if exist('tcdf','file')==2, pWelch=2*(1-tcdf(abs(tWelch),df)); end
    bootDiff=zeros(nBootstrap,1);
    for z=1:nBootstrap
        bootDiff(z)=mean(a(randi(na,na,1)))-mean(b(randi(nb,nb,1)));
    end
    ci = prctile_base(bootDiff,[2.5 97.5]);
    pooled=[a(:);b(:)]; observed=abs(delta); extreme=0;
    for z=1:nPermutation
        idx=randperm(na+nb);
        d=mean(pooled(idx(1:na)))-mean(pooled(idx(na+1:end)));
        extreme=extreme+(abs(d)>=observed-1e-14);
    end
    pPermutation=(extreme+1)/(nPermutation+1);
    sp=sqrt(((na-1)*sa^2+(nb-1)*sb^2)/(na+nb-2));
    hedgesG=(delta/sp)*(1-3/(4*(na+nb)-9));
    row=table(f,na,nb,100*ma,100*sa,100*mb,100*sb,100*delta,100*ci(1),100*ci(2),tWelch,df,pWelch,pPermutation,hedgesG, ...
        'VariableNames',{'NominalPhi','N_VGPP','N_RSA','MeanVGPP_RSD_pct','SDVGPP_pct','MeanRSA_RSD_pct','SDRSA_pct', ...
        'DifferenceVGPPminusRSA_pctPoints','BootstrapCI_Lower_pctPoints','BootstrapCI_Upper_pctPoints', ...
        'WelchT','WelchDF','WelchP_ifAvailable','PermutationP_twoSided','HedgesG'});
    Result=[Result;row]; %#ok<AGROW>
    Raw=[Raw;table(repmat(f,na,1), repmat("VGPP",na,1), (1:na)',100*a(:), 'VariableNames',{'NominalPhi','Method','WithinMethodRow','LocalRSD_pct'})]; %#ok<AGROW>
    Raw=[Raw;table(repmat(f,nb,1), repmat("RSA-PBC",nb,1), (1:nb)',100*b(:), 'VariableNames',{'NominalPhi','Method','WithinMethodRow','LocalRSD_pct'})]; %#ok<AGROW>
end
writetable(Result,fullfile(outputDir,'R2_C28_statistical_comparison.csv'));
writetable(Raw,fullfile(outputDir,'R2_C28_included_RSD_values.csv'));

fid=fopen(fullfile(outputDir,'R2_C28_interpretation.txt'),'w');
fprintf(fid,'Input: %s\n%s\n',inputFile,unitInfo);
fprintf(fid,'Independent two-sample comparison; two-sided random-label permutation test, %d permutations.\n',nPermutation);
fprintf(fid,'Percentile bootstrap CI for difference in means, %d resamples.\n',nBootstrap);
fprintf(fid,'Only successfully attained target fractions 0.40 and 0.50 were compared.\n');
fprintf(fid,'No equivalence margin was prespecified; a nonsignificant p value does NOT establish equivalence.\n');
fprintf(fid,'No multiple-comparison adjustment; two descriptive target-specific comparisons.\n');
for k=1:height(Result)
 fprintf(fid,'phi=%.2f delta(VGPP-RSA)=%.5f percentage points; bootstrap 95%% CI [%.5f,%.5f]; permutation p=%.5g\n', ...
  Result.NominalPhi(k),Result.DifferenceVGPPminusRSA_pctPoints(k),Result.BootstrapCI_Lower_pctPoints(k),Result.BootstrapCI_Upper_pctPoints(k),Result.PermutationP_twoSided(k));
end
fclose(fid);

% One publication-friendly mean/SD figure; preserves older manuscript figures.
figure('Color','w','Position',[80 80 820 520]); hold on; box on;
q=Result.NominalPhi;
h1=errorbar(q-.005,Result.MeanVGPP_RSD_pct,Result.SDVGPP_pct,'-o','LineWidth',1.4,'MarkerSize',7);
h2=errorbar(q+.005,Result.MeanRSA_RSD_pct,Result.SDRSA_pct,'-s','LineWidth',1.4,'MarkerSize',7);
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',1.1);
xlim([0.375 0.525]);xticks(targets);
xlabel('Nominal particle area fraction');ylabel('Local area-fraction RSD (%)');
legend([h1 h2],{'VGPP','RSA-PBC'},'Location','best');
exportgraphics(gcf,fullfile(picsDir,'R2_C28_VGPP_RSA_RSD_comparison.png'),'Resolution',300);
close(gcf);
disp(Result); fprintf('Outputs: %s\n',outputDir);

%% Local helper: percentile without Statistics Toolbox
function q = prctile_base(x,p)
 x=sort(x(:));n=numel(x);q=zeros(size(p));
 for i=1:numel(p)
  h=1+(n-1)*p(i)/100;lo=floor(h);hi=ceil(h);
  q(i)=x(lo)+(h-lo)*(x(hi)-x(lo));
 end
end
