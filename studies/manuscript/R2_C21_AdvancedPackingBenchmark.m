%% R2_C21_AdvancedPackingBenchmark.m
% Scientific Reports - Reviewer 2, Comment 21
% Independent event-driven Lubachevsky-Stillinger (LS) hard-disk growth
% MATLAB R2024b, standalone SCRIPT with local functions; no toolboxes.
%
% IMPORTANT DESIGN/INTERPRETATION
% 1) This is a direct, transparent event-driven LS reference implementation:
%    particles MOVE ballistically, their radii GROW at a fixed rate, and
%    pair collision times are solved analytically. Collisions reflect the
%    *surface-relative* normal velocity. No soft overlap forces or VGPP
%    collective relaxation are employed.
% 2) This is NOT an optimized heap/cell-list LS implementation. For N~700,
%    runtime can be substantial. Do NOT infer LS's general speed from this
%    pedagogical MATLAB implementation. Both codes should run on the same
%    computer, with transparent timing and identical success definitions.
% 3) Event/time-budget exhaustion is CENSORED, NOT proof of unreachability.
% 4) Other published algorithms (Yan, Kou, Tian, Zeng) are NOT reproduced
%    by this file. This file provides one independent LS reference only.
%
% OUTPUT: revision_results/R2_C21_LS/* (CSV, MAT, PNG)
% USAGE: first run quickTestMode=true; then false. Set production settings
%        based on preliminary cost. All saved runs are checkpointed.
%
clc; clear; close all;

%% User settings
quickTestMode = false;         % TRUE first. Switch to FALSE for full targets.
rerunCompleted = false;       % reproducible checkpointing
L = 10;
rFinal = 0.18;
phiTargets = [0.60, 0.70, 0.74];
nSeeds = 3;                   % initial pilot, NOT statistical superiority
masterSeed = 210521;
growthRate = 0.002;          % model lengths per event-driven time unit
initialRadiusFraction = 0.10;
initialSpeed = 1.0;           % model lengths per time unit
velocityRefreshEvents = 0;  % 0 = pure deterministic LS with collision-only velocities
maxEvents = 200000;          % explicit computational censoring threshold
maxWallSeconds = 1800;       % per realization; change only if documented
maxTravelDt = 0.25;          % radius-independent maximum ballistic horizon
contactTol = 1e-10;          % for event root filtering
finalOverlapTol = 1e-8;      % same geometric criterion as VGPP
saveCoordinates = true;
figureDPI = 250;
outputFolder = fullfile('1st revision', 'revision_results', 'R2_C21_LS');
if ~exist(outputFolder,'dir'), mkdir(outputFolder); end

if quickTestMode
    phiTargets = 0.10;        % small regression/validation case
    nSeeds = 1;
    rFinal = 0.18;
    growthRate = 0.01;
    maxEvents = 12000;
    maxWallSeconds = 90;
end

% Save configuration separately for audit.
settingsPath = fullfile(outputFolder,'R2_C21_settings.txt');
fid = fopen(settingsPath,'w');
fprintf(fid,'MATLAB=%s\n',version);
fprintf(fid,'quickTestMode=%d\nL=%.12g\nrFinal=%.12g\n',quickTestMode,L,rFinal);
fprintf(fid,'phiTargets=%s\nnSeeds=%d\nmasterSeed=%d\n',mat2str(phiTargets),nSeeds,masterSeed);
fprintf(fid,'growthRate=%.12g\ninitialRadiusFraction=%.12g\ninitialSpeed=%.12g\n',growthRate,initialRadiusFraction,initialSpeed);
fprintf(fid,'maxEvents=%d\nmaxWallSeconds=%.3f\nmaxTravelDt=%.12g\n',maxEvents,maxWallSeconds,maxTravelDt);
fprintf(fid,'finalOverlapTol=%.12g\n',finalOverlapTol);
fprintf(fid,'LS is collision-driven growing-disk reference, not optimized.\n');
fclose(fid);

% Validate the actual event solver on an analytic two-body contact case.
ls_selftest();

records = table();
checkpointPath = fullfile(outputFolder, 'R2_C21_LS_raw.csv');
if isfile(checkpointPath)
    records = readtable(checkpointPath,'TextType','string');
end

fprintf('\nR2-C21 LS reference benchmark: quickTest=%d\n',quickTestMode);
for ip = 1:numel(phiTargets)
    phiNominal = phiTargets(ip);
    N = round(phiNominal*L^2/(pi*rFinal^2));
    phiRealizable = N*pi*rFinal^2/L^2;
    for k = 1:nSeeds
        seed = masterSeed + 10000*ip + k + 1000000*double(quickTestMode);
        caseID = sprintf('LS_Q%d_P%06d_S%03d',quickTestMode,round(phiNominal*1e5),k);
        if ~rerunCompleted && ~isempty(records) && ismember('CaseID',records.Properties.VariableNames) && any(string(records.CaseID)==string(caseID))
            fprintf('Skipping checkpointed %s\n',caseID);
            continue;
        end
        fprintf('\n%s: N=%d, target phi=%.8f\n',caseID,N,phiRealizable);
        rng(seed,'twister');
        ticRun = tic;
        x = L*rand(N,2);                   % zero-radius point centers initially
        v = randn(N,2);
        v = v-mean(v,1);
        v = initialSpeed*v / sqrt(mean(sum(v.^2,2)));
        r = 0;                           % exact zero radius avoids initial overlaps
        timeNow = 0;
        nColl = 0; nHorizon = 0; nRebuilds = 0;
        status = "UNKNOWN";
        % Conservative kinetic + growth horizon: many close events near high phi
        % A time horizon is *not* a discretized collision step: event times
        % within the horizon are solved analytically.
        while true
            if r >= rFinal - 1e-12
                r = rFinal;
                status = "REACHED_RADIUS";
                break;
            end
            if nColl >= maxEvents
                status = "CENSORED_EVENT_BUDGET";
                break;
            end
            if toc(ticRun) > maxWallSeconds
                status = "CENSORED_WALL_BUDGET";
                break;
            end
            horizon = min(maxTravelDt, (rFinal-r)/growthRate);
            % Compute earliest pair collision across all periodic images.
            % This full global scan is exact for linearly moving centers and
            % uniformly growing disks over each ballistic horizon; no soft
            % overlap corrections occur.
            [dt, ia, ib, nvec] = ls_earliest_collision(x,v,L,r,growthRate,horizon,contactTol);
            if ~isfinite(dt)
                dt = horizon;
                x = mod(x+dt*v,L);
                r = r + growthRate*dt;
                timeNow = timeNow+dt;
                nHorizon = nHorizon+1;
            else
                % Keep contact position and radius at mathematically same t.
                x = mod(x+dt*v,L);
                r = r + growthRate*dt;
                timeNow = timeNow+dt;
                relativeNormal = dot(v(ia,:)-v(ib,:),nvec);
                % Growing-disk elastic collision: relative normal velocity
                % after = 2*(2*growthRate) - relativeNormal before.
                % For equal masses, add/subtract (2*growthRate-w_n)*n.
                impulse = (2*growthRate-relativeNormal)*nvec;
                v(ia,:) = v(ia,:) + impulse;
                v(ib,:) = v(ib,:) - impulse;
                nColl = nColl+1;
                % Velocity refresh is OFF by default; not part of baseline.
                if velocityRefreshEvents>0 && mod(nColl,velocityRefreshEvents)==0
                    v = v-mean(v,1);
                    v = initialSpeed*v/sqrt(mean(sum(v.^2,2)));
                    nRebuilds=nRebuilds+1;
                end
            end
            if any(~isfinite(x),'all') || any(~isfinite(v),'all') || ~isfinite(r)
                status="NUMERICAL_FAILURE";
                break;
            end
        end
        elapsed = toc(ticRun);
        phiAchieved = N*pi*r^2/L^2;
        [minGap, nOverlap] = ls_periodic_audit(x,L,r,finalOverlapTol);
        exactTarget = status=="REACHED_RADIUS" && nOverlap==0 && ...
            abs(r-rFinal)<=1e-12 && abs(phiAchieved-phiRealizable)<1e-9;
        if status=="REACHED_RADIUS" && ~exactTarget
            status="FAILED_GEOMETRIC_AUDIT";
        end
        fprintf('  %s | phi=%.7f | collisions=%d | wall=%.1fs | minGap=%.3g | overlaps=%d\n', ...
            status,phiAchieved,nColl,elapsed,minGap,nOverlap);
        if saveCoordinates
            coordPath=fullfile(outputFolder,[caseID,'_coords.csv']);
            writetable(table(x(:,1),x(:,2),r*ones(N,1), ...
                'VariableNames',{'x','y','radius'}),coordPath);
        end
        matPath=fullfile(outputFolder,[caseID,'.mat']);
        save(matPath,'x','v','r','L','N','phiNominal','phiRealizable', ...
            'phiAchieved','elapsed','nColl','nHorizon','minGap','nOverlap', ...
            'status','seed','growthRate','timeNow','-v7.3');
        newRow=table(string(caseID),logical(quickTestMode),phiNominal,phiRealizable,N,seed, ...
            string(status),logical(exactTarget),phiAchieved,r,elapsed,nColl,nHorizon, ...
            minGap,nOverlap,growthRate,timeNow, ...
            'VariableNames',{'CaseID','QuickTest','TargetPhi','RealizablePhi','N','Seed', ...
            'Status','ExactTargetSuccess','AchievedPhi','FinalRadius','Runtime_s', ...
            'Collisions','FreeFlightSteps','MinPeriodicGap','OverlapsBeyondTol', ...
            'GrowthRate','PhysicalTime'});
        if isempty(records), records=newRow; else, records=[records;newRow]; end %#ok<AGROW>
        writetable(records,checkpointPath);
    end
end

% Summaries only of applicable test mode; censored attempts remain counted.
T=records(records.QuickTest==logical(quickTestMode),:);
if ~isempty(T)
    phiUniq=unique(T.TargetPhi);
    summary=table();
    for j=1:numel(phiUniq)
        S=T(abs(T.TargetPhi-phiUniq(j))<1e-10,:);
        ss=S(S.ExactTargetSuccess,:);
        rt=NaN;
        if ~isempty(ss),rt=mean(ss.Runtime_s);end
        add=table(phiUniq(j),S.RealizablePhi(1),height(S),sum(S.ExactTargetSuccess), ...
            sum(startsWith(S.Status,"CENSORED")),rt, ...
            'VariableNames',{'NominalPhi','RealizablePhi','NRuns','NExactSuccess', ...
            'NCensored','MeanSuccessfulRuntime_s'});
        summary=[summary;add]; %#ok<AGROW>
    end
    writetable(summary,fullfile(outputFolder,'R2_C21_LS_summary.csv'));
    disp(summary);
end
fprintf('\nReference code complete. DO NOT interpret censored as packing failure.\n');
fprintf('Compare VGPP and LS ONLY after matching acceptance, environment, budgets, and samples.\n');

%% ----------- Local helper: analytic first collision of growing disks ---------
function [bestT,besti,bestj,bestN]=ls_earliest_collision(x,v,L,r,gamma,horizon,tol)
% Vectorized per-particle all-image event scan, collision-driven LS.
N=size(x,1);
if max(sqrt(sum(v.^2,2)))*horizon > L/4
    error('Ballistic horizon too large for periodic image enumeration. Reduce maxTravelDt.');
end
bestT=Inf; besti=0; bestj=0; bestN=[NaN NaN];
R=2*r; G=2*gamma;
for i=1:N-1
    jj=(i+1:N)';
    w=v(i,:)-v(jj,:);
    p0=x(i,:)-x(jj,:);
    ww=sum(w.^2,2);
    for sx=-1:1
        for sy=-1:1
            p=p0+L*[sx sy];
            cc=sum(p.^2,2)-R^2;
            if any(cc < -max(tol,1e-12))
                error('Numerical interpenetration detected before event search.');
            end
            aa=ww-G^2;
            bb=2*(sum(p.*w,2)-R*G);
            disc=bb.^2-4*aa.*cc;
            t=Inf(size(aa));
            linear=abs(aa)<1e-14 & bb<0;
            t(linear)=-cc(linear)./bb(linear);
            quad=~(abs(aa)<1e-14) & disc>=0;
            if any(quad)
                sd=sqrt(max(0,disc(quad)));
                a=aa(quad); b=bb(quad);
                t1=(-b-sd)./(2*a);
                t2=(-b+sd)./(2*a);
                t1(t1<=tol)=Inf;
                t2(t2<=tol)=Inf;
                t(quad)=min(t1,t2);
            end
            eligible=t>tol & t<=horizon & t<bestT & isfinite(t);
            ids=find(eligible);
            if isempty(ids),continue;end
            cp=p(ids,:)+t(ids).*w(ids,:);
            dd=sqrt(sum(cp.^2,2));
            approach=sum(cp.*w(ids,:),2)./max(dd,realmin);
            keep=approach<G-1e-12;
            ids=ids(keep);
            if isempty(ids),continue;end
            [newT,k]=min(t(ids));
            pos=p(ids(k),:)+newT*w(ids(k),:);
            bestT=newT; besti=i; bestj=jj(ids(k)); bestN=pos/norm(pos);
        end
    end
end
end

%% ----------- Independent continuous-coordinate pairwise acceptance audit ----
function [minimumGap,nBad]=ls_periodic_audit(x,L,r,tol)
N=size(x,1);minimumGap=Inf;nBad=0;
for i=1:N-1
    dx=abs(x(i+1:end,1)-x(i,1));dx=min(dx,L-dx);
    dy=abs(x(i+1:end,2)-x(i,2));dy=min(dy,L-dy);
    gaps=hypot(dx,dy)-2*r;
    minimumGap=min(minimumGap,min(gaps));
    nBad=nBad+sum(gaps < -tol);
end
end

function ls_selftest()
% 2 equal growing hard disks in one dimension, p=(1,0), w=(-1,0), r=.1,
% gamma=.05. Contact time solves 1-t=0.2+0.1t => 8/11.
x=[1 1;0 1];v=[-0.5 0;0.5 0];
[t,i,j,n]=ls_earliest_collision(x,v,10,.1,.05,1,1e-12);
assert(abs(t-8/11)<1e-8 && i==1 && j==2 && abs(norm(n)-1)<1e-12, ...
    'Analytic LS growing-disk collision regression failed.');
[gap,bad]=ls_periodic_audit([1 1;1.4 1],10,.2,1e-8);
assert(abs(gap)<1e-12 && bad==0,'Periodic overlap audit regression failed.');
end
