%% R1_C4_PackingCeiling.m
% Scientific Reports | Reviewer #1, Comment 4
% "Reinforce that 0.70 is only the tested limit rather than VGPP's
%  theoretical maximum. Discuss factors determining its practical packing
%  ceiling."
%
% PURPOSE
%  - Test VGPP ABOVE the manuscript's previously highest TESTED fraction 0.70.
%  - Focused resumable diagnostic: 8 baseline + 4 sensitivity attempts in full mode.
%  - Identify the largest EXACT-target, fixed-final-radius, periodically
%    audited packing OBSERVED within defined solver budgets and realizations.
%  - Test how grid resolution, collective solver budget, and prescribed
%    interparticle clearance, and starting seed radius influence
%    *practical* reachability.
%  - Save actual best geometries with coordinates and periodic validation.
%
% This is NOT a mathematical determination of VGPP's global maximum and
% failure of a stochastic trial is NOT proof of an unattainable fraction.
% No theoretical ceiling is inferred by extrapolating a finite sweep.
% Comparative claims about physical packing limits are not made.
%
% CORE PACKING IMPLEMENTATION
%  Instrumented copy of the fixed-final-radius seed -> growth -> periodic
%  audit kernel in the author's code_refined(5).m, via R1_C3_RuntimeScatter.m.
%  Exact same seedScale rule, insertion parameters, growth/backtracking,
%  fallback deletion and audit tolerance as base code. All helper functions
%  are included locally below. No hexagonal fallback; no particle inflation
%  past original final r=0.18. Note: original distance-transform selection
%  is a GRID APPROXIMATION, but final overlap acceptance uses continuous
%  pair distances and periodic minimum images, not the grid.
%  Stable weighted sampling normalization avoids avoidable underflow.
%
% MATLAB R2024b; Image Processing Toolbox required for bwdist/imregionalmax.
% ALL stage durations are exploratory instrumentation; measured runtime
% depends on hardware. No manuscript changes until outputs are validated.
%
% OUTPUTS: revision_results/R1_C4_raw.csv; R1_C4_summary.csv;
%   R1_C4_run_settings.txt; R1_C4_best_valid_coordinates.csv;
%   R1_C4_best_valid_configuration.mat;
%   pics/R1_C4_success.png; R1_C4_runtime.png; R1_C4_best_RVE.png;
%   pics/R1_C4_sensitivity.png
% -------------------------------------------------------------------------
clc; clear; close all;

resultsDir=fullfile('revision_results','R1_C4_focused'); picsDir=fullfile('pics','R1_C4_focused');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end
if ~exist(picsDir,'dir'), mkdir(picsDir); end
if exist('bwdist','file')~=2 || exist('imregionalmax','file')~=2
    error('R1_C4:MissingToolbox', ...
        'Image Processing Toolbox is needed for bwdist and imregionalmax.');
end
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',13);

%% Original parameters: values from code_refined(5).m (multi-target mode)
p=struct();
p.Lx=10; p.Ly=10; p.Nx=250; p.Ny=250;
p.particleMode="monodisperse";
p.rMono=0.18; p.rMean=0.18; p.rStd=0.035;
p.rMin=0.10; p.rMax=0.28;
p.minGap=0;
p.batchSize=80; p.nCandidates=4000;
p.distanceUpdateInterval=5;
p.useLargestVoidInsertion=true; p.useWeightedSampling=true;
p.seedPlacementMode='full_vgpp'; % REQUIRED by insertion kernel
p.largestVoidFraction=0.70;
p.largestVoidMinSeparation=0.75*p.rMono;
p.maxCycles=10000; p.maxFailedCycles=500;
p.relaxationMovesPerCycle=5000;
p.growthCompressionRelaxIterations=800;
p.growthCompressionRelaxFactor=0.45;
p.overlapAuditTolerance=1e-8;
p.growthCompressionOverlapTolerance=p.overlapAuditTolerance;
p.fixedRadiusSeedVF=0.50;
p.fixedRadiusGrowthScale=1.0025;
p.fixedRadiusMaxGrowthSteps=300;
p.fixedRadiusMaxBacktracks=14;
p.fixedRadiusMinGrowthScale=1.00002;
p.targetVFTolerance=1e-6;
p.useHexagonalFallback=false; p.useGrowthCompression=false;
p.stopAtSpecifiedTarget=false;

%% EXPERIMENT CONTROL
% Start with the quick test. Set false for final numerical analysis.
quickTestMode = false;                  % <<< CHANGE TO false AFTER QUICK TEST
runSensitivityStudy = true;           % limited one-factor-at-a-time tests

% 0.70 is a positive control, not a presupposed ceiling. Others test above.
baselineFractions = [0.70 0.72 0.74 0.76];
baselineRuns = 2;

% Each sensitivity changes ONLY the indicated factor while maintaining the
% same target N and, for a given (nominal VF, run), the same RNG seed.
% These are illustrative bounded settings, not solver optimization.
sensitivityFractions = 0.74;
sensitivityRuns = 1;

if quickTestMode
    baselineFractions = [0.70 0.74];
    baselineRuns=1;
    sensitivityFractions=[];
    sensitivityRuns=0;
end

% Conditions: baseline; finer candidate grid; larger growth solver budget;
% and positive *physical* clearance (decreases available center space).
caseNames={'Baseline','FinerGrid','MoreRelaxation','PositiveGap','SmallerSeeds'};
caseGrid=[250 375 250 250 250];
caseRelax=[800 800 3000 800 800];    % kernel uses max(value, 1500)
caseGap=[0 0 0 0.02*p.rMono 0];
caseSeedVF=[0.50 0.50 0.50 0.50 0.40];
if ~runSensitivityStudy || quickTestMode
    caseNames=caseNames(1);
end

seedBase=15400;
rawFile=fullfile(resultsDir,'R1_C4_raw.csv');
if exist(rawFile,'file')==2
    raw=readtable(rawFile,'TextType','string');
    fprintf('Resuming from %d existing attempt(s).\n',height(raw));
else
    raw=table();
end
bestAreaFraction=-inf;
bestRecord=struct();
bestParticleTable=table();
bestMat=fullfile(resultsDir,'R1_C4_best_valid_configuration.mat');
if exist(bestMat,'file')==2
    saved=load(bestMat,'bestRecord','bestParticleTable');
    bestRecord=saved.bestRecord;
    bestParticleTable=saved.bestParticleTable;
    bestAreaFraction=bestRecord.AchievedVF;
end
startAll=tic;

fprintf('\n===========================================================\n');
fprintf('R1 / Comment 4: FIXED-FINAL-RADIUS PRACTICAL PACKING TEST\n');
fprintf('Final r = %.5f, square L = %.2f, audit tolerance = %.2e\n', ...
    p.rMono,p.Lx,p.overlapAuditTolerance);
fprintf('Hexagonal fallback DISABLED. Initial r may be smaller than final r.\n');
fprintf('Targets refer to DISCRETE area fraction from integer N.\n');
fprintf('quickTestMode=%d; sensitivity=%d\n',quickTestMode,runSensitivityStudy);
fprintf('===========================================================\n');

%% Run all cases: save after EVERY attempt to preserve previous progress.
for icase=1:numel(caseNames)
    caseName=string(caseNames{icase});
    pNow=p;
    pNow.Nx=caseGrid(icase);pNow.Ny=caseGrid(icase);
    pNow.growthCompressionRelaxIterations=caseRelax(icase);
    pNow.minGap=caseGap(icase);
    pNow.fixedRadiusSeedVF=caseSeedVF(icase);
    if icase==1
        fractions=baselineFractions; nRuns=baselineRuns;
    else
        fractions=sensitivityFractions; nRuns=sensitivityRuns;
    end
    fprintf('\nCASE %s | grid %d x %d | relax limit %d | gap %.5g | seed VF %.3f\n', ...
        caseName,pNow.Nx,pNow.Ny,...
        max(pNow.growthCompressionRelaxIterations,1500),pNow.minGap, ...
        pNow.fixedRadiusSeedVF);

    for iv=1:numel(fractions)
        targetVF=fractions(iv);
        nTarget=round(targetVF*p.Lx*p.Ly/(pi*p.rMono^2));
        exactVF=nTarget*pi*p.rMono^2/(p.Lx*p.Ly);
        for ir=1:nRuns
            % Resume by exact scenario, target, and replicate; preserve completed runs.
            if ~isempty(raw) && any(raw.Scenario==caseName & ...
                    abs(raw.NominalVF-targetVF)<1e-12 & raw.Run==ir)
                fprintf('  Reused existing %s, phi=%.3f, run=%d\n',caseName,targetVF,ir);
                continue;
            end
            % Match seeds for a given area fraction / run ACROSS conditions.
            seedNow=seedBase+round(10000*targetVF)+ir;
            pr=pNow;pr.targetParticleCount=nTarget;
            pr.activeTargetVF=exactVF;
            rng(seedNow,'twister');
            finished=false; succeeded=false; audited=false;
            status=""; nFinal=NaN; achieved=NaN;
            nOverlapping=NaN; minClearance=NaN; maxOverlap=NaN;
            seedRadius=NaN; growSteps=NaN; validRadius=false;
            diag=R1C4_zero_diagnostics();
            wall=tic;
            try
                result=R1C4_run_instrumented(pr);
                finished=true;
                nFinal=height(result.particleTable);
                achieved=result.maxAreaFraction;
                nOverlapping=result.nOverlappingPairs;
                minClearance=result.minClearance;
                maxOverlap=result.maxOverlap;
                seedRadius=result.seedRadius;
                growSteps=result.growthSteps;
                diag=result.profile;
                validRadius=nFinal>0 && all(abs(result.particleTable.radius-p.rMono)<1e-12);
                % Audit is INDEPENDENT of solver convergence. A nonempty
                % terminal state is valid only if ALL radii and pairs pass.
                audited=(nFinal>0) && validRadius && ~result.hasOverlap && ...
                    isfinite(achieved) && all(isfinite(result.particleTable.x)) && ...
                    all(isfinite(result.particleTable.y));
                succeeded=audited && nFinal==nTarget && ...
                    abs(achieved-exactVF)<1e-10 && logical(result.reachedTarget);
                if succeeded
                    status="Exact target reached; strict periodic audit passed";
                elseif audited
                    status="Target not reached; smaller valid geometry obtained";
                else
                    status="No audited final-radius configuration";
                end
                % Store best independently VERIFIED packing, even if below
                % target due to deletions. Never store invalid geometry.
                if audited && achieved>bestAreaFraction+1e-12
                    bestAreaFraction=achieved;
                    bestRecord=struct('Scenario',caseName,'NominalVF',targetVF,...
                        'ExactTargetVF',exactVF,'Seed',seedNow,'TargetN',nTarget,...
                        'FinalN',nFinal,'AchievedVF',achieved,'TargetSuccess',succeeded,...
                        'GridN',pNow.Nx,'MinGap',pNow.minGap,...
                        'FinalRadius',p.rMono,'AuditTolerance',p.overlapAuditTolerance,...
                        'MinClearance',minClearance,'MaxOverlap',maxOverlap);
                    bestParticleTable=result.particleTable;
                    writetable(bestParticleTable,fullfile(resultsDir, ...
                        'R1_C4_best_valid_coordinates.csv'));
                    save(fullfile(resultsDir,'R1_C4_best_valid_configuration.mat'), ...
                        'bestRecord','bestParticleTable');
                end
            catch ME
                status=string(ME.identifier)+": "+string(ME.message);
                fprintf(2,'ERROR %s phi=%.3f seed=%d: %s\n',...
                    caseName,targetVF,seedNow,ME.message);
            end
            wall_s=toc(wall);
            tab=table(caseName,targetVF,exactVF,nTarget,ir,seedNow,...
                pNow.Nx,pNow.minGap,pNow.fixedRadiusSeedVF,...
                max(pNow.growthCompressionRelaxIterations,1500),...
                finished,succeeded,audited,validRadius,nFinal,achieved,...
                seedRadius,growSteps,minClearance,maxOverlap,nOverlapping,...
                diag.SeedFinalN,diag.SeedCycles,diag.SeedFailedCycles,...
                diag.GrowthAcceptedSteps,diag.GrowthFailedSolves,...
                diag.GrowthBacktrackingAttempts,diag.GrowthSolverIterations,...
                diag.FallbackDeletions,diag.GrowthStopReason,wall_s,status,...
                'VariableNames',{'Scenario','NominalVF','ExactVF','TargetN',...
                'Run','Seed','GridSize','MinGap','InitialSeedVF','RelaxMaxIterations',...
                'Completed','TargetSuccess','AuditedGeometry','RadiusExact',...
                'FinalN','FinalVF','SeedRadius','GrowthSteps',...
                'MinClearance','MaxOverlap','NOverlappingPairs',...
                'SeedFinalN','SeedCycles','SeedFailedCycles',...
                'GrowthAcceptedSteps','FailedSolverCalls',...
                'BacktrackingAttempts','GrowthSolverIterations',...
                'FallbackDeletions','GrowthStopReason','WallTime_s','Status'});
            raw=[raw;tab]; %#ok<AGROW>
            writetable(raw,rawFile);
            fprintf('  phi=%.3f (exact %.6f) seed=%d success=%d geom=%d N=%d/%d time=%.1fs\n', ...
                targetVF,exactVF,seedNow,succeeded,audited,nFinal,nTarget,wall_s);
        end
    end
end

%% Empirical summary. Failure is NOT evidence of an absolute upper bound.
summary=table();
uniqueCases=unique(raw.Scenario,'stable');
for ic=1:numel(uniqueCases)
  sub=raw(raw.Scenario==uniqueCases(ic),:);
  vfs=unique(sub.NominalVF);
  for j=1:numel(vfs)
    a=sub(abs(sub.NominalVF-vfs(j))<1e-12,:);
    n=height(a);
    validA=a.FinalVF(a.AuditedGeometry);
    if isempty(validA)
        largestAudited=NaN;
    else
        largestAudited=max(validA);
    end
    t=table(uniqueCases(ic),vfs(j),a.ExactVF(1),a.TargetN(1),...
        n,sum(a.TargetSuccess),100*mean(a.TargetSuccess),...
        sum(a.AuditedGeometry),largestAudited,mean(a.WallTime_s),...
        median(a.WallTime_s),max(a.WallTime_s),...
        mean(a.SeedCycles,'omitnan'),...
        mean(a.FailedSolverCalls,'omitnan'),...
        mean(a.BacktrackingAttempts,'omitnan'),...
        mean(a.FallbackDeletions,'omitnan'),...
        'VariableNames',{'Scenario','NominalVF','ExactVF','TargetN',...
        'NTrials','NExactSuccess','ExactSuccessRate_pct',...
        'NAuditedGeometry','HighestAuditedVF','MeanRuntime_s',...
        'MedianRuntime_s','MaxRuntime_s','MeanSeedCycles',...
        'MeanFailedSolverCalls','MeanBacktracks','MeanFallbackDeletions'});
    summary=[summary;t]; %#ok<AGROW>
  end
end
writetable(summary,fullfile(resultsDir,'R1_C4_summary.csv'));

%% Figures: success rate vs tested targets; cost; parameter sensitivity.
base=summary(summary.Scenario=="Baseline",:);
fig=figure('Color','w','Position',[90 90 980 530]);
plot(base.NominalVF,base.ExactSuccessRate_pct,'-ko','LineWidth',2,...
  'MarkerFaceColor','k','MarkerSize',7);hold on;grid on;box on;
xticks(base.NominalVF);ylim([0 105]);
xlabel('Target area fraction');ylabel('Exact-target success rate (%)');
title('VGPP: tested target reachability, not a mathematical ceiling');
exportgraphics(fig,fullfile(picsDir,'R1_C4_success.png'),'Resolution',300);
close(fig);

fig=figure('Color','w','Position',[90 90 980 530]);hold on;grid on;box on;
for j=1:height(base)
    a=raw(raw.Scenario=="Baseline" & ...
        abs(raw.NominalVF-base.NominalVF(j))<1e-12,:);
    scatter(base.NominalVF(j)*ones(height(a),1),a.WallTime_s/60,...
        35,'filled','MarkerFaceAlpha',0.65,'HandleVisibility','off');
end
plot(base.NominalVF,base.MedianRuntime_s/60,'-ks','LineWidth',1.8,...
  'MarkerFaceColor','k','DisplayName','Median');
xlabel('Target area fraction');ylabel('Wall-clock runtime (minutes)');
title('Computational cost over attempted dense configurations');
legend('Location','northwest');
exportgraphics(fig,fullfile(picsDir,'R1_C4_runtime.png'),'Resolution',300);
close(fig);

if numel(uniqueCases)>1
    fig=figure('Color','w','Position',[90 90 1050 530]);hold on;grid on;box on;
    for ic=1:numel(uniqueCases)
        t=summary(summary.Scenario==uniqueCases(ic),:);
        plot(t.NominalVF,t.ExactSuccessRate_pct,'-o',...
            'LineWidth',1.6,'DisplayName',char(uniqueCases(ic)));
    end
    xlabel('Target area fraction');ylabel('Exact-target success (%)');
    ylim([0 105]);title('Sensitivity to numerical controls and physical clearance');
    legend('Location','best');
    exportgraphics(fig,fullfile(picsDir,'R1_C4_sensitivity.png'),'Resolution',300);
    close(fig);
end

% Plot a genuine result, NEVER a schematic or synthetic lattice.
if ~isempty(bestParticleTable) && ~isempty(bestRecord)
  fig=figure('Color','w','Position',[80 80 760 760]);
  ax=axes(fig); hold(ax,'on');axis(ax,'equal');box(ax,'on');
  rectangle(ax,'Position',[0 0 p.Lx p.Ly],'EdgeColor','k','LineWidth',1.4);
  theta=linspace(0,2*pi,65);
  for k=1:height(bestParticleTable)
      x=bestParticleTable.x(k); y=bestParticleTable.y(k);
      rr=bestParticleTable.radius(k);
      for sx=-1:1
          for sy=-1:1
              xx=x+sx*p.Lx+rr*cos(theta);
              yy=y+sy*p.Ly+rr*sin(theta);
              if any(xx>=0 & xx<=p.Lx & yy>=0 & yy<=p.Ly)
                  patch(ax,xx,yy,[0.35 0.55 0.75],'EdgeColor','none',...
                    'FaceAlpha',0.75);
              end
          end
      end
  end
  xlim(ax,[0 p.Lx]);ylim(ax,[0 p.Ly]);
  xlabel(ax,'x');ylabel(ax,'y');
  title(ax,sprintf('Best audited VGPP: %.5f (N=%d, target success=%d)',...
    bestRecord.AchievedVF,bestRecord.FinalN,bestRecord.TargetSuccess));
  exportgraphics(fig,fullfile(picsDir,'R1_C4_best_RVE.png'),'Resolution',300);
  close(fig);
end

%% Concise audit report: do not label failures as an upper packing bound.
fid=fopen(fullfile(resultsDir,'R1_C4_run_settings.txt'),'w');
if fid>0
 fprintf(fid,'MATLAB: %s\n',version);
 fprintf(fid,'Platform: %s\n',computer);
 fprintf(fid,'Machine details: ADD CPU/RAM manually for reproducibility\n');
 fprintf(fid,'quickTestMode=%d; runSensitivityStudy=%d\n',quickTestMode,runSensitivityStudy);
 fprintf(fid,'Baseline targets: %s; nRuns=%d\n',mat2str(baselineFractions),baselineRuns);
 fprintf(fid,'Sensitivity targets: %s; nRuns=%d\n',mat2str(sensitivityFractions),sensitivityRuns);
 fprintf(fid,'L=[%.4f, %.4f]; radius=%.6f; baseline grid=%dx%d\n',...
    p.Lx,p.Ly,p.rMono,p.Nx,p.Ny);
 fprintf(fid,'hexagonal fallback disabled; no inflation above prescribed r\n');
 fprintf(fid,'Pairwise overlap tolerance %.3e (in length units)\n',p.overlapAuditTolerance);
 fprintf(fid,'baseline seed VF %.5f; growth steps %d; retries %d\n',...
    p.fixedRadiusSeedVF,p.fixedRadiusMaxGrowthSteps,p.fixedRadiusMaxBacktracks);
 fprintf(fid,'Highest successful EXACT TARGET (all cases): ');
 successes=raw(raw.TargetSuccess,:);
 if isempty(successes)
     fprintf(fid,'NONE\n');
 else
     fprintf(fid,'%.8f\n',max(successes.ExactVF));
 end
 if isfinite(bestAreaFraction)
     fprintf(fid,'Highest final VERIFIED GEOMETRY (all cases): %.8f\n',bestAreaFraction);
 else
     fprintf(fid,'No verified final geometry found.\n');
 end
 fprintf(fid,'Failure never establishes theoretical unreachability.\n');
 fprintf(fid,'PositiveGap modifies a physical constraint, not merely numerical cost.\n');
 fclose(fid);
end
fprintf('\nCOMPLETED in %.1f min. Raw/summary saved to %s.\n',toc(startAll)/60,resultsDir);
if isfinite(bestAreaFraction)
 fprintf('Highest audited actual area fraction = %.6f\n',bestAreaFraction);
end
if quickTestMode
 fprintf('QUICK TEST ONLY: set quickTestMode=false for research results.\n');
end

%% Below: standalone local kernels, derived from code_refined(5).m.
function d=R1C4_zero_diagnostics()
d=struct('SeedInsertion_s',0,'SeedMap_s',0,'SeedCandidate_s',0,...
 'SeedChecks_s',0,'SeedPaint_s',0,'SeedCycles',0,'SeedFailedCycles',0,...
 'SeedConsecutiveFailuresMax',0,'SeedMapRebuilds',0,...
 'SeedCandidateProposals',0,'SeedRejectExisting',0,'SeedRejectBatch',0,...
 'SeedRelaxEvents',0,'SeedRelax_s',0,'SeedFinalN',0,...
 'Growth_s',0,'GrowthTrialSolves',0,'GrowthFailedSolves',0,...
 'GrowthBacktrackingAttempts',0,'GrowthAcceptedSteps',0,...
 'GrowthSolver_s',0,'GrowthSolverIterations',0,'GrowthPairChecks',0,...
 'GrowthMaxSolveIterations',0,'GrowthStopReason',"",...
 'ReachedFinalRadiusBeforeFallback',false,...
 'Fallback_s',0,'FallbackDeletions',0,'FallbackSolverIterations',0,...
 'FallbackSolverPairChecks',0,'FallbackReachedFinalRadius',false,...
 'Audit_s',0,'TotalUnaccounted_s',NaN);
end

function v=R1C4_percentile(x,p)
x=sort(x(isfinite(x)));n=numel(x);
if n==0,v=NaN;return;end
q=1+(n-1)*p/100;i=floor(q);j=ceil(q);
v=x(i)+(x(j)-x(i))*(q-i);
end

function rho=R1C4_spearman(a,b)
idx=isfinite(a)&isfinite(b);a=a(idx);b=b(idx);
if numel(a)<5 || numel(unique(a))<2 || numel(unique(b))<2
  rho=NaN;return;
end
x=R1C4_rank(a);y=R1C4_rank(b);cc=corrcoef(x,y);
rho=cc(1,2);
end

function r=R1C4_rank(x)
x=x(:);[s,idx]=sort(x);r=zeros(size(x));n=numel(x);i=1;
while i<=n
 j=i;while j<n && s(j+1)==s(i),j=j+1;end
 r(idx(i:j))=(i+j)/2;i=j+1;
end
end

function result=R1C4_run_instrumented(params)
% SAME BASELINE seed / growth / fallback decisions as uploaded code_refined(5).m.
% Only add measurements and counters around existing operations.
Lx=params.Lx; Ly=params.Ly; A=Lx*Ly; rFinal=params.rMono;
nTarget=params.targetParticleCount;
targetVF=nTarget*pi*rFinal^2/A;
seedScale=min(1,sqrt(params.fixedRadiusSeedVF/targetVF));
seedScale=max(seedScale,0.65); rSeed=rFinal*seedScale;
ps=params;ps.rMono=rSeed;ps.rMean=rSeed;
ps.useGrowthCompression=false;ps.useHexagonalFallback=false;
ps.stopAtSpecifiedTarget=false;
Nx=params.Nx;Ny=params.Ny;
x=linspace(0,Lx,Nx);y=linspace(0,Ly,Ny);
[X,Y]=meshgrid(x,y);
pixelSize=mean([Lx/(Nx-1),Ly/(Ny-1)]);
cx=[];cy=[];r=[];
occ=false(Ny,Nx);dist=inf(Ny,Nx);area=0;
failed=0;cyc=0; d=R1C4_zero_diagnostics();

% Seed packing: record all failures, including those later recovered.
tSeed=tic;
while numel(r)<nTarget && cyc<params.maxCycles && failed<params.maxFailedCycles
  cyc=cyc+1;
  updateMap=cyc==1 || mod(cyc-1,params.distanceUpdateInterval)==0;
  ps.batchSize=min(params.batchSize,nTarget-numel(r));
  [cx,cy,r,occ,dist,nAdded,area,z]=attempt_DTPAA_batch_insertion(...
      cx,cy,r,occ,dist,area,X,Y,pixelSize,ps,updateMap);
  d.SeedMap_s=d.SeedMap_s+z.Map_s;
  d.SeedCandidate_s=d.SeedCandidate_s+z.Candidate_s;
  d.SeedChecks_s=d.SeedChecks_s+z.Checks_s;
  d.SeedPaint_s=d.SeedPaint_s+z.Paint_s;
  d.SeedMapRebuilds=d.SeedMapRebuilds+z.MapRebuilds;
  d.SeedCandidateProposals=d.SeedCandidateProposals+z.Proposed;
  d.SeedRejectExisting=d.SeedRejectExisting+z.RejectExisting;
  d.SeedRejectBatch=d.SeedRejectBatch+z.RejectBatch;
  if nAdded==0
    failed=failed+1;d.SeedFailedCycles=d.SeedFailedCycles+1;
    d.SeedConsecutiveFailuresMax=max(d.SeedConsecutiveFailuresMax,failed);
    if mod(failed,10)==0 && numel(r)>1
      oldx=cx;oldy=cy;
      tRelax=tic;
      [cx,cy]=stochastic_relaxation_PBC(cx,cy,r,ps,0.08*rSeed,...
        min(3000,params.relaxationMovesPerCycle));
      occ=rebuild_periodic_occupancy(X,Y,cx,cy,r,ps);
      dist=inf(size(occ));
      d.SeedRelax_s=d.SeedRelax_s+toc(tRelax);
      d.SeedRelaxEvents=d.SeedRelaxEvents+1;
      if isempty(cx),cx=oldx;cy=oldy;end
    end
  else
    failed=0;
  end
end
d.SeedInsertion_s=toc(tSeed)-d.SeedRelax_s;
d.SeedCycles=cyc;d.SeedFinalN=numel(r);

% Fixed-radius collective solver and backtracking; count explicit pair checks.
pg=params;
pg.growthCompressionOverlapTolerance=params.overlapAuditTolerance;
pg.growthCompressionRelaxIterations=max(params.growthCompressionRelaxIterations,1500);
reachedFinal=false;growthSteps=0;
if isempty(r),r=zeros(0,1);end
tGrow=tic;
while ~isempty(r) && mean(r)<rFinal*(1-1e-12) && ...
      growthSteps<params.fixedRadiusMaxGrowthSteps
    growthSteps=growthSteps+1;
    trialScale=min(params.fixedRadiusGrowthScale,rFinal/mean(r));
    accepted=false;
    for ib=0:params.fixedRadiusMaxBacktracks
      if trialScale<=1,break;end
      rTrial=min(r*trialScale,rFinal);
      cxTrial=cx;cyTrial=cy;
      ts=tic;
      [cxTrial,cyTrial,success,~,itUsed,checks]=...
        resolve_overlaps_repulsive_PBC(cxTrial,cyTrial,rTrial,pg);
      d.GrowthSolver_s=d.GrowthSolver_s+toc(ts);
      d.GrowthTrialSolves=d.GrowthTrialSolves+1;
      d.GrowthSolverIterations=d.GrowthSolverIterations+itUsed;
      d.GrowthPairChecks=d.GrowthPairChecks+checks;
      d.GrowthMaxSolveIterations=max(d.GrowthMaxSolveIterations,itUsed);
      if success
        cx=cxTrial;cy=cyTrial;r=rTrial;accepted=true;
        d.GrowthAcceptedSteps=d.GrowthAcceptedSteps+1;
        break;
      end
      d.GrowthFailedSolves=d.GrowthFailedSolves+1;
      d.GrowthBacktrackingAttempts=d.GrowthBacktrackingAttempts+1;
      trialScale=1+(trialScale-1)*0.5;
      if trialScale<params.fixedRadiusMinGrowthScale,break;end
    end
    if ~accepted
      d.GrowthStopReason="Radius increment rejected despite backtracking";
      break;
    end
end
if growthSteps>=params.fixedRadiusMaxGrowthSteps && ...
  ~isempty(r) && mean(r)<rFinal*(1-1e-12)
 d.GrowthStopReason="Max growth steps reached";
end
if ~isempty(r) && mean(r)>=rFinal*(1-1e-10)
    r(:)=rFinal;
    ts=tic;
    [cx,cy,sFinal,~,itUsed,checks]=resolve_overlaps_repulsive_PBC(cx,cy,r,pg);
    d.GrowthSolver_s=d.GrowthSolver_s+toc(ts);
    d.GrowthTrialSolves=d.GrowthTrialSolves+1;
    d.GrowthSolverIterations=d.GrowthSolverIterations+itUsed;
    d.GrowthPairChecks=d.GrowthPairChecks+checks;
    d.GrowthMaxSolveIterations=max(d.GrowthMaxSolveIterations,itUsed);
    if ~sFinal,d.GrowthFailedSolves=d.GrowthFailedSolves+1;end
    reachedFinal=sFinal;
end
d.Growth_s=toc(tGrow);
d.ReachedFinalRadiusBeforeFallback=reachedFinal;
if strlength(d.GrowthStopReason)==0
    if reachedFinal,d.GrowthStopReason="Reached final radius";
    else,d.GrowthStopReason="Final radius not reached";end
end

% Original deletion fallback (can lose exact requested N): independent time.
tFallback=tic;
while ~reachedFinal && numel(r)>1
    kRemove=most_crowded_particle_PBC(cx,cy,rFinal,Lx,Ly);
    cx(kRemove)=[];cy(kRemove)=[];
    r=rFinal*ones(numel(cx),1);
    d.FallbackDeletions=d.FallbackDeletions+1;
    [cx,cy,sFinal,~,itUsed,checks]=resolve_overlaps_repulsive_PBC(cx,cy,r,pg);
    d.FallbackSolverIterations=d.FallbackSolverIterations+itUsed;
    d.FallbackSolverPairChecks=d.FallbackSolverPairChecks+checks;
    if sFinal,reachedFinal=true;end
    if numel(r)<0.90*nTarget,break;end
end
d.Fallback_s=toc(tFallback);
d.FallbackReachedFinalRadius=reachedFinal;

% Independent final exact PBC audit, all particle pairs.
ta=tic;
[hasOverlap,maxOverlap,minClearance,nPairs]=final_overlap_audit_PBC(cx,cy,r,params);
d.Audit_s=toc(ta);
reachedTarget=reachedFinal && ~hasOverlap && numel(r)>=nTarget;
particleTable=table((1:numel(r))',cx,cy,r,...
    'VariableNames',{'ID','x','y','radius'});
result=struct('particleTable',particleTable,'maxAreaFraction',...
    sum(pi*r.^2)/A,'reachedTarget',reachedTarget,...
    'targetN',nTarget,'seedRadius',rSeed,'seedScale',seedScale,...
    'growthSteps',growthSteps,'hasOverlap',hasOverlap,...
    'maxOverlap',maxOverlap,'minClearance',minClearance,...
    'nOverlappingPairs',nPairs,'profile',d);
end

%% Baseline packing kernels retained below (fully embedded)
function [cxList, cyList, rList, occCore, distOcc, nInserted, currentArea, stage] = ...
    attempt_DTPAA_batch_insertion(cxList, cyList, rList, occCore, distOcc, ...
    currentArea, X, Y, pixelSize, params, updateDistanceMap)

    Lx = params.Lx;
    Ly = params.Ly;

    nInserted = 0;
    stage = struct('Map_s',0,'Candidate_s',0,'Checks_s',0,'Paint_s',0,...
        'MapRebuilds',0,'Proposed',0,'RejectExisting',0,'RejectBatch',0);

    %% Trial particle radius
    rCurrent = generate_particle_radius(params);

    %% Update distance transform if needed
    timer = tic;
    % Uniform-domain proposals do not require a distance transform. Avoid
    % charging that comparator for a map it does not use.
    if updateDistanceMap && ~strcmp(params.seedPlacementMode,'uniform_domain')
        if any(occCore(:))
            stage.MapRebuilds = stage.MapRebuilds+1;
            distOcc = periodic_distance_transform(occCore, pixelSize);
        else
            distOcc = inf(size(occCore));
        end
    end

    stage.Map_s = toc(timer);
    timer = tic;

    %% Identify eligible proposals. Uniform arm intentionally uses the whole
    % periodic Cartesian grid (not the occupancy mask); every proposed
    % center must still pass the identical exact periodic non-overlap test.
    if strcmp(params.seedPlacementMode, 'uniform_domain')
        admissible = true(size(occCore));
    elseif any(occCore(:))
        admissible = distOcc >= (rCurrent + params.minGap);
    else
        admissible = true(size(occCore));
    end

    admissibleIndices = find(admissible);

    if isempty(admissibleIndices)
        return;
    end

    %% Select candidate points
    nAvailable = numel(admissibleIndices);
    nCandidates = min(nAvailable, params.nCandidates);

    if nCandidates < 1
        return;
    end

    if strcmp(params.seedPlacementMode, 'uniform_domain')
        % Full-domain trial/rejection arm, on the common candidate grid.
        selectedIndices = admissibleIndices(randperm(nAvailable, nCandidates));
    else
        selectedIndices = select_candidate_indices_LV_hybrid( ...
            admissibleIndices, distOcc, nCandidates, params, X, Y);
    end

    candX = X(selectedIndices);
    candY = Y(selectedIndices);
    stage.Candidate_s = toc(timer);
    timer = tic;

    %% Build one insertion batch
    batchX = [];
    batchY = [];
    batchR = [];

    for k = 1:numel(candX)
        stage.Proposed = stage.Proposed + 1;

        xNew = candX(k);
        yNew = candY(k);
        rNew = rCurrent;

        %% Exact periodic check against existing particles
        if ~isempty(cxList)

            dExisting = periodic_center_distance( ...
                xNew, yNew, cxList, cyList, Lx, Ly);

            requiredDistanceExisting = rList + rNew + params.minGap;

            if any(dExisting < requiredDistanceExisting)
                stage.RejectExisting = stage.RejectExisting+1;
                continue;
            end
        end

        %% Exact periodic check against current batch particles
        if ~isempty(batchX)

            dBatch = periodic_center_distance( ...
                xNew, yNew, batchX, batchY, Lx, Ly);

            requiredDistanceBatch = batchR + rNew + params.minGap;

            if any(dBatch < requiredDistanceBatch)
                stage.RejectBatch = stage.RejectBatch+1;
                continue;
            end
        end

        %% Optional cap at the specified target VF
        if isfield(params,'stopAtSpecifiedTarget') && params.stopAtSpecifiedTarget && ...
                isfield(params,'activeTargetVF') && isfinite(params.activeTargetVF)
            domainAreaLocal = params.Lx * params.Ly;
            targetAreaLocal = params.activeTargetVF * domainAreaLocal;
            batchAreaNow = sum(pi * batchR.^2);
            projectedArea = currentArea + batchAreaNow + pi*rNew^2;
            if projectedArea > targetAreaLocal + params.targetVFTolerance*domainAreaLocal
                continue;
            end
        end

        %% Accept candidate into the batch
        batchX(end+1,1) = xNew;
        batchY(end+1,1) = yNew;
        batchR(end+1,1) = rNew;

        if numel(batchX) >= params.batchSize
            break;
        end
    end

    stage.Checks_s = toc(timer);
    timer = tic;

    %% Insert accepted batch into occupancy map
    if isempty(batchX)
        return;
    end

    for p = 1:numel(batchX)

        xp = batchX(p);
        yp = batchY(p);
        rp = batchR(p);

        particleMask = periodic_particle_mask(X, Y, xp, yp, rp, Lx, Ly);
        occCore = occCore | particleMask;

        cxList(end+1,1) = xp;
        cyList(end+1,1) = yp;
        rList(end+1,1)  = rp;

        currentArea = currentArea + pi * rp^2;
        nInserted = nInserted + 1;
    end
    stage.Paint_s = toc(timer);
end

%% ------------------------------------------------------------------------

function [cxList, cyList] = stochastic_relaxation_PBC(cxList, cyList, rList, params, moveAmplitude, nMoves)

    Lx = params.Lx;
    Ly = params.Ly;
    minGap = params.minGap;

    nParticles = numel(rList);

    if nParticles < 2
        return;
    end

    for moveID = 1:nMoves

        %% Randomly select one particle
        i = randi(nParticles);

        oldX = cxList(i);
        oldY = cyList(i);

        %% Random displacement
        dxMove = moveAmplitude * (2*rand - 1);
        dyMove = moveAmplitude * (2*rand - 1);

        newX = mod(oldX + dxMove, Lx);
        newY = mod(oldY + dyMove, Ly);

        %% Check overlap with all other particles
        otherIdx = true(nParticles,1);
        otherIdx(i) = false;

        xOthers = cxList(otherIdx);
        yOthers = cyList(otherIdx);
        rOthers = rList(otherIdx);

        d = periodic_center_distance(newX, newY, xOthers, yOthers, Lx, Ly);
        requiredDistance = rList(i) + rOthers + minGap;

        if all(d >= requiredDistance)
            cxList(i) = newX;
            cyList(i) = newY;
        end
    end
end

%% ------------------------------------------------------------------------

function occCore = rebuild_periodic_occupancy(X, Y, cxList, cyList, rList, params)

    occCore = false(size(X));

    for i = 1:numel(rList)
        particleMask = periodic_particle_mask(X, Y, cxList(i), cyList(i), ...
            rList(i), params.Lx, params.Ly);
        occCore = occCore | particleMask;
    end
end


%% ------------------------------------------------------------------------

function [cxList, cyList, success, maxOverlap, iterUsed, pairChecks] = ...
    resolve_overlaps_repulsive_PBC(cxList, cyList, rList, params)

    % Collective overlap resolution using pairwise repulsive displacements.
    % This is the core high-VF compaction step. It repeatedly detects all
    % overlapping pairs under periodic boundary conditions and pushes them
    % apart by a fraction of their overlap distance.

    Lx = params.Lx;
    Ly = params.Ly;
    minGap = params.minGap;
    maxIter = params.growthCompressionRelaxIterations;
    relaxFactor = params.growthCompressionRelaxFactor;
    tol = params.growthCompressionOverlapTolerance;

    nParticles = numel(rList);
    success = false;
    maxOverlap = inf;
    iterUsed = 0;
    pairChecks = 0;

    if nParticles < 2
        success = true;
        maxOverlap = 0.0;
        return;
    end

    % Pair index list is fixed because the particle number does not change.
    [I, J] = find(triu(true(nParticles), 1));

    for iter = 1:maxIter
        iterUsed = iter;
        pairChecks = pairChecks + numel(I);

        dx = cxList(J) - cxList(I);
        dy = cyList(J) - cyList(I);

        % Minimum-image convention for PBC.
        dx = dx - Lx * round(dx / Lx);
        dy = dy - Ly * round(dy / Ly);

        d = sqrt(dx.^2 + dy.^2);
        requiredDistance = rList(I) + rList(J) + minGap;
        overlap = requiredDistance - d;

        active = overlap > tol;

        if ~any(active)
            success = true;
            maxOverlap = 0.0;
            return;
        end

        maxOverlap = max(overlap(active));

        ia = I(active);
        ja = J(active);
        dxa = dx(active);
        dya = dy(active);
        da  = d(active);
        oa  = overlap(active);

        % Avoid division by zero for coincident centers.
        zeroD = da < eps;
        if any(zeroD)
            randomAngle = 2*pi*rand(sum(zeroD),1);
            dxa(zeroD) = cos(randomAngle);
            dya(zeroD) = sin(randomAngle);
            da(zeroD) = 1.0;
        end

        ux = dxa ./ da;
        uy = dya ./ da;

        % Split the required correction equally between the two particles.
        shift = 0.5 * relaxFactor * oa;

        moveX = zeros(nParticles,1);
        moveY = zeros(nParticles,1);

        moveX = moveX + accumarray(ia, -shift .* ux, [nParticles, 1], @sum, 0);
        moveY = moveY + accumarray(ia, -shift .* uy, [nParticles, 1], @sum, 0);
        moveX = moveX + accumarray(ja,  shift .* ux, [nParticles, 1], @sum, 0);
        moveY = moveY + accumarray(ja,  shift .* uy, [nParticles, 1], @sum, 0);

        % Limit very large moves to avoid oscillation.
        moveMagnitude = sqrt(moveX.^2 + moveY.^2);
        maxMove = 0.25 * mean(rList);
        tooLarge = moveMagnitude > maxMove;

        if any(tooLarge)
            moveX(tooLarge) = moveX(tooLarge) ./ moveMagnitude(tooLarge) * maxMove;
            moveY(tooLarge) = moveY(tooLarge) ./ moveMagnitude(tooLarge) * maxMove;
        end

        cxList = mod(cxList + moveX, Lx);
        cyList = mod(cyList + moveY, Ly);
    end

    % Final check after all iterations.
    pairChecks = pairChecks + numel(I);
    dx = cxList(J) - cxList(I);
    dy = cyList(J) - cyList(I);
    dx = dx - Lx * round(dx / Lx);
    dy = dy - Ly * round(dy / Ly);
    d = sqrt(dx.^2 + dy.^2);
    requiredDistance = rList(I) + rList(J) + minGap;
    overlap = requiredDistance - d;
    maxOverlap = max([0; overlap]);
    success = maxOverlap <= tol;
end

%% ------------------------------------------------------------------------

function distCenter = periodic_distance_transform(occCore, pixelSize)

    [Ny, Nx] = size(occCore);

    occTiled = repmat(occCore, 3, 3);
    distTiled = bwdist(occTiled) * pixelSize;

    rowStart = Ny + 1;
    rowEnd   = 2 * Ny;

    colStart = Nx + 1;
    colEnd   = 2 * Nx;

    distCenter = distTiled(rowStart:rowEnd, colStart:colEnd);
end

%% ------------------------------------------------------------------------

function mask = periodic_particle_mask(X, Y, cx, cy, r, Lx, Ly)

    dx = abs(X - cx);
    dy = abs(Y - cy);

    dx = min(dx, Lx - dx);
    dy = min(dy, Ly - dy);

    mask = dx.^2 + dy.^2 <= r^2;
end

%% ------------------------------------------------------------------------

function d = periodic_center_distance(xNew, yNew, xList, yList, Lx, Ly)

    dx = abs(xList - xNew);
    dx = min(dx, Lx - dx);

    dy = abs(yList - yNew);
    dy = min(dy, Ly - dy);

    d = sqrt(dx.^2 + dy.^2);
end

%% ------------------------------------------------------------------------

function selectedIndices = select_candidate_indices_LV_hybrid( ...
    admissibleIndices, distOcc, nCandidates, params, X, Y)

    % Hybrid largest-void + weighted-random candidate selection.
    %
    % A prescribed fraction of candidate centers is selected from local
    % maxima of the periodic distance-transform field. These points
    % represent the largest remaining admissible voids. The remaining
    % candidates are selected using the weighted random selector to preserve
    % stochastic morphology.

    nAvailable = numel(admissibleIndices);
    nCandidates = min(nCandidates, nAvailable);

    if nAvailable == 0 || nCandidates == 0
        selectedIndices = [];
        return;
    end

    %% If largest-void insertion is disabled, use the original selector

    if ~isfield(params, 'useLargestVoidInsertion') || ~params.useLargestVoidInsertion

        selectedIndices = select_candidate_indices( ...
            admissibleIndices, distOcc, nCandidates, params.useWeightedSampling);

        return;
    end

    %% If distance field is not informative, use original selector

    weights = distOcc(admissibleIndices);

    if isempty(weights) || all(isinf(weights)) || all(isnan(weights))

        selectedIndices = select_candidate_indices( ...
            admissibleIndices, distOcc, nCandidates, params.useWeightedSampling);

        return;
    end

    %% Number of largest-void candidates

    if isfield(params, 'largestVoidFraction')
        largestVoidFraction = params.largestVoidFraction;
    else
        largestVoidFraction = 0.70;
    end

    largestVoidFraction = max(0, min(1, largestVoidFraction));

    nVoid = round(largestVoidFraction * nCandidates);
    nVoid = min(max(nVoid, 1), nCandidates);

    %% Largest-void candidate centers

    voidIndices = find_largest_void_candidates( ...
        admissibleIndices, distOcc, nVoid, params, X, Y);

    %% Weighted-random candidate centers for the remaining quota

    nRemaining = nCandidates - numel(voidIndices);

    if nRemaining > 0

        remainingPool = setdiff(admissibleIndices, voidIndices, 'stable');

        if ~isempty(remainingPool)
            randomIndices = select_candidate_indices( ...
                remainingPool, distOcc, nRemaining, params.useWeightedSampling);
        else
            randomIndices = [];
        end
    else
        randomIndices = [];
    end

    %% Combine candidates

    selectedIndices = unique([voidIndices(:); randomIndices(:)], 'stable');

    %% Fill if unique operation reduced candidate count

    if numel(selectedIndices) < nCandidates

        missingCount = nCandidates - numel(selectedIndices);
        remainingPool = setdiff(admissibleIndices, selectedIndices, 'stable');

        if ~isempty(remainingPool)
            missingCount = min(missingCount, numel(remainingPool));
            extraLocal = randperm(numel(remainingPool), missingCount);
            extraIndices = remainingPool(extraLocal);
            selectedIndices = [selectedIndices(:); extraIndices(:)];
        end
    end
end

%% ------------------------------------------------------------------------

function selectedIndices = select_candidate_indices(admissibleIndices, distOcc, ...
    nCandidates, useWeightedSampling)

    nAvailable = numel(admissibleIndices);
    nCandidates = min(nCandidates, nAvailable);

    if nAvailable == 0 || nCandidates == 0
        selectedIndices = [];
        return;
    end

    %% Uniform random selection
    if ~useWeightedSampling
        selectedLocal = randperm(nAvailable, nCandidates);
        selectedIndices = admissibleIndices(selectedLocal);
        return;
    end

    %% Weighted random selection
    weights = distOcc(admissibleIndices);
    weights = weights(:);

    if isempty(weights) || all(isinf(weights)) || all(isnan(weights))
        selectedLocal = randperm(nAvailable, nCandidates);
        selectedIndices = admissibleIndices(selectedLocal);
        return;
    end

    finiteWeights = weights(isfinite(weights));

    if isempty(finiteWeights)
        selectedLocal = randperm(nAvailable, nCandidates);
        selectedIndices = admissibleIndices(selectedLocal);
        return;
    end

    maxFiniteWeight = max(finiteWeights);

    weights(isinf(weights)) = maxFiniteWeight;
    weights(isnan(weights)) = 0;
    weights(weights < 0) = 0;

    if sum(weights) <= 0
        selectedLocal = randperm(nAvailable, nCandidates);
        selectedIndices = admissibleIndices(selectedLocal);
        return;
    end

    selectedLocal = weighted_sample_without_replacement(weights, nCandidates);
    selectedIndices = admissibleIndices(selectedLocal);
end

%% ------------------------------------------------------------------------

function selected = weighted_sample_without_replacement(weights, n)

    weights = weights(:);

    nTotal = numel(weights);
    n = min(n, nTotal);

    weights(isnan(weights)) = 0;

    finiteWeights = weights(isfinite(weights));

    if isempty(finiteWeights)
        selected = randperm(nTotal, n)';
        return;
    end

    maxFiniteWeight = max(finiteWeights);

    weights(isinf(weights)) = maxFiniteWeight;
    weights(weights < 0) = 0;

    if sum(weights) <= 0
        selected = randperm(nTotal, n)';
        return;
    end

    % Stable Efraimidis-Spirakis weighted sampling without replacement.
    % The power form U^(1/w) in the uploaded baseline underflows to zero
    % for the roughly 6e4 grid points used here. The algebraically
    % equivalent log-key form avoids this artifact. Apply identically to
    % all arms that use weighted proposals, including Full_VGPP.
    weights = weights ./ sum(weights);
    positive = weights > 0;
    keys = -inf(nTotal,1);
    u = max(rand(nnz(positive),1), realmin('double'));
    keys(positive) = log(u)./weights(positive);
    [~, order] = sort(keys, 'descend');
    selected = order(1:n);
end


%% ------------------------------------------------------------------------

function voidIndices = find_largest_void_candidates( ...
    admissibleIndices, distOcc, nVoid, params, X, Y)

    % Detects and ranks the largest void centers using local maxima of
    % the periodic Euclidean distance-transform field.

    if nVoid <= 0
        voidIndices = [];
        return;
    end

    %% Build admissible mask

    admissibleMask = false(size(distOcc));
    admissibleMask(admissibleIndices) = true;

    %% Clean distance field

    distClean = distOcc;

    finiteValues = distClean(isfinite(distClean));

    if isempty(finiteValues)
        voidIndices = [];
        return;
    end

    maxFinite = max(finiteValues);

    distClean(isinf(distClean)) = maxFinite;
    distClean(isnan(distClean)) = 0;

    %% Detect local maxima

    if exist('imregionalmax', 'file') == 2
        localMaxima = imregionalmax(distClean);
    else
        % Fallback if imregionalmax is unavailable.
        % The sorted admissible distance list still prioritizes largest voids.
        localMaxima = false(size(distClean));
    end

    voidMask = localMaxima & admissibleMask;
    candidateVoidIndices = find(voidMask);

    %% If local maxima are not found, use largest admissible distance values

    if isempty(candidateVoidIndices)
        [~, order] = sort(distClean(admissibleIndices), 'descend');
        sortedIndices = admissibleIndices(order);
    else
        voidScores = distClean(candidateVoidIndices);
        [~, order] = sort(voidScores, 'descend');
        sortedIndices = candidateVoidIndices(order);
    end

    %% Spatially thin largest void candidates to avoid clustering

    if isfield(params, 'largestVoidMinSeparation') && params.largestVoidMinSeparation > 0

        voidIndices = spatially_thin_void_candidates( ...
            sortedIndices, nVoid, params.largestVoidMinSeparation, ...
            X, Y, params.Lx, params.Ly);
    else
        nSelect = min(nVoid, numel(sortedIndices));
        voidIndices = sortedIndices(1:nSelect);
    end
end

%% ------------------------------------------------------------------------

function selectedVoidIndices = spatially_thin_void_candidates( ...
    sortedIndices, nVoid, minSeparation, X, Y, Lx, Ly)

    % Keeps largest void centers while enforcing a minimum separation
    % between selected candidate centers under periodic boundary conditions.

    selectedVoidIndices = [];
    selectedX = [];
    selectedY = [];

    for k = 1:numel(sortedIndices)

        idx = sortedIndices(k);

        xNew = X(idx);
        yNew = Y(idx);

        if isempty(selectedX)

            selectedVoidIndices(end+1,1) = idx;
            selectedX(end+1,1) = xNew;
            selectedY(end+1,1) = yNew;

        else

            d = periodic_center_distance(xNew, yNew, selectedX, selectedY, Lx, Ly);

            if all(d >= minSeparation)

                selectedVoidIndices(end+1,1) = idx;
                selectedX(end+1,1) = xNew;
                selectedY(end+1,1) = yNew;
            end
        end

        if numel(selectedVoidIndices) >= nVoid
            break;
        end
    end

    %% If thinning was too strict, fill from the ranked list

    if numel(selectedVoidIndices) < nVoid

        remaining = setdiff(sortedIndices, selectedVoidIndices, 'stable');
        nMissing = nVoid - numel(selectedVoidIndices);
        nMissing = min(nMissing, numel(remaining));

        selectedVoidIndices = [selectedVoidIndices(:); remaining(1:nMissing)];
    end
end

%% ------------------------------------------------------------------------

function r = generate_particle_radius(params)

    if params.particleMode == "monodisperse"
        r = params.rMono;

    elseif params.particleMode == "polydisperse"
        accepted = false;

        while ~accepted
            rTry = params.rMean + params.rStd * randn;

            if rTry >= params.rMin && rTry <= params.rMax
                r = rTry;
                accepted = true;
            end
        end

    else
        error('particleMode must be either "monodisperse" or "polydisperse".');
    end
end

%% ------------------------------------------------------------------------

function [hasOverlap, maxOverlap, minClearance, nOverlappingPairs, worstPair] = ...
    final_overlap_audit_PBC(cxList, cyList, rList, params)

    % Explicit all-pair overlap audit under periodic boundary conditions.
    %
    % This function is the final authority for geometric validity. It checks
    % every particle pair using the minimum-image convention, so it detects
    % both ordinary in-domain overlaps and mirrored overlaps across opposite
    % periodic edges.

    Lx = params.Lx;
    Ly = params.Ly;
    minGap = params.minGap;

    nParticles = numel(rList);

    hasOverlap = false;
    maxOverlap = 0;
    minClearance = inf;
    nOverlappingPairs = 0;
    worstPair = [NaN, NaN];

    if nParticles < 2
        minClearance = NaN;
        return;
    end

    for i = 1:nParticles-1

        xNew = cxList(i);
        yNew = cyList(i);

        xOthers = cxList(i+1:end);
        yOthers = cyList(i+1:end);
        rOthers = rList(i+1:end);

        d = periodic_center_distance(xNew, yNew, xOthers, yOthers, Lx, Ly);

        requiredDistance = rList(i) + rOthers + minGap;
        clearance = d - requiredDistance;
        overlap = requiredDistance - d;

        minClearance = min(minClearance, min(clearance));

        if isfield(params,'overlapAuditTolerance')
            auditTol = params.overlapAuditTolerance;
        else
            auditTol = 1.0e-8;
        end
        overlappingLocal = overlap > auditTol;
        nOverlappingPairs = nOverlappingPairs + sum(overlappingLocal);

        [localMaxOverlap, localIdx] = max([0; overlap(:)]);

        if localMaxOverlap > maxOverlap
            maxOverlap = localMaxOverlap;

            if localIdx > 1
                worstPair = [i, i + localIdx - 1];
            end
        end
    end

    hasOverlap = nOverlappingPairs > 0;

    if isfield(params,'overlapAuditTolerance')
        auditTol = params.overlapAuditTolerance;
    else
        auditTol = 1.0e-8;
    end
    if maxOverlap <= auditTol
        maxOverlap = 0;
    end

    if isinf(minClearance)
        minClearance = NaN;
    end
end


%% ------------------------------------------------------------------------
% Literature-inspired statistical validation block
% -------------------------------------------------------------------------
% This block adds post-processing plots commonly used in recent RVE-generation
% studies: nearest-neighbor distance/orientation, Ripley's K/L function,
% radial distribution function, local volume-fraction uniformity, and a
% compact statistical validation table. It does not modify the generated
% geometry; it only analyzes the final VGPP and RSA-PBC particle coordinates.

function idxRemove = most_crowded_particle_PBC(cx,cy,rFinal,Lx,Ly)

    N = numel(cx);
    nearestClearance = inf(N,1);

    for i = 1:N
        d = periodic_center_distance(cx(i),cy(i),cx,cy,Lx,Ly);
        d(i) = inf;
        nearestClearance(i) = min(d - 2*rFinal);
    end

    [~,idxRemove] = min(nearestClearance);
end
