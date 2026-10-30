%% R2_C21_Consolidated.m
% Scientific Reports, Reviewer 2 Comment 21.
% Collision-driven, growing-hard-disk LS reference with binary event heap,
% lazy velocity-change invalidation, and finite periodic flight horizons.
% No overlap-force relaxation, spring forces, or optimization-based packing.
% MATLAB R2024b; standalone script with local functions; no toolboxes.
%
% IMPORTANT: This is a research reference implementation, not a validated
% production-grade or fully optimized Lubachevsky-Stillinger package.
% Its collisions and event scheduling MUST pass checks and its outcomes
% MUST be inspected before being reported. Failed/censored trials do NOT
% establish an LS packing ceiling or VGPP performance superiority.
clc; clear; close all;

%% USER CONFIGURATION
% Modes: audit (default) / smoke / pilot / production / resume
% Consolidated: default READ-ONLY audit; optional smoke/pilot/production/resume.
% Existing LS dynamics are retained from validated V5 continuation engine.
% Resume reads an existing MAT checkpoint. It NEVER reconstructs velocity
% information from a coordinate CSV, and it never silently restarts a case.
mode = "audit";       % DEFAULT: read-only, no dynamics, no checkpoint changes
resumeCaseIDs = "LSOPT_production_P074000_S002"; % used ONLY if mode="resume"
rerunCompleted = false;
projectFolder = '.';
if ~isfolder(projectFolder)
    projectFolder = '.';
end
if ~isfolder(projectFolder), projectFolder=fileparts(mfilename('fullpath')); end
outputFolder=fullfile(projectFolder,'revision_results','R2_C21_LS_optimized');
% Override outputFolder here if the archive is stored somewhere else.
L=10; rFinal=.18; masterSeed=210521;
gamma=.002; initialSpeed=1;
maxWallSeconds=1800;                 % additional wall time for the resumed run
maxCollisionEvents=1000000;         % CUMULATIVE total events from t=0
maxHorizon=.50;
geometricTolerance=1e-8; contactEpsilon=2e-11;
writeCoordinates=true; saveMAT=true;
perEventAudit = mode=="smoke";

% An explicit read-only pathway preserves all existing result/checkpoint files.
if mode=="audit"
    ls_c21_audit_archive(outputFolder,projectFolder,L,rFinal,geometricTolerance);
    return;
end

switch mode
    case "smoke"
        phiTargets=.10; nSeeds=1; gamma=.01; maxWallSeconds=90;
        maxCollisionEvents=20000;
    case "pilot"
        phiTargets=.60; nSeeds=1;
    case {"production","resume"}
        phiTargets=[.60 .70 .74]; nSeeds=3;
    otherwise
        error('Invalid mode. Use smoke, pilot, production or resume.');
end
if mode=="resume"
    assert(~isempty(resumeCaseIDs) && all(strlength(resumeCaseIDs)>0), ...
        'Specify one or more resumeCaseIDs.');
    assert(all(startsWith(resumeCaseIDs,"LSOPT_production_")), ...
        'Every resumeCaseID must refer to an existing production case.');
    assert(numel(unique(resumeCaseIDs))==numel(resumeCaseIDs), ...
        'Duplicate resumeCaseIDs are not allowed.');
end

if ~isfolder(outputFolder), mkdir(outputFolder); end
settings=fullfile(outputFolder, "R2_C21_optimized_settings_V5_"+mode+".txt");
fid=fopen(settings,'w'); assert(fid>0,'Could not write settings file.');
fprintf(fid,'MATLAB=%s\nmode=%s\nL=%.12g\nrFinal=%.12g\nphiTargets=%s\nnSeeds=%d\n', ...
    version,char(mode),L,rFinal,mat2str(phiTargets),nSeeds);
fprintf(fid,'masterSeed=%d\ngrowthRate=%.12g\ninitialSpeed=%.12g\n',masterSeed,gamma,initialSpeed);
fprintf(fid,'maxWallSeconds=%.12g\nmaxCollisionEvents=%d\nmaxHorizon=%.12g\n', ...
    maxWallSeconds,maxCollisionEvents,maxHorizon);
fprintf(fid,'geometricTolerance=%.12g\nperEventAudit=%d\n',geometricTolerance,perEventAudit);
fprintf(fid,'resumeCaseIDs=%s\n',strjoin(cellstr(resumeCaseIDs),', '));
fprintf(fid,'Collision-driven growing hard disks; no force-based overlap correction.\n');
fclose(fid);

% Fail fast if the analytic collision, velocity update, periodic audit or heap
% ordering tests fail. These tests are necessary, but not sufficient, for
% full research-grade LS validation.
ls_opt_selftest();

rawPath=fullfile(outputFolder,'R2_C21_optimized_raw.csv');
if mode=="resume"
    assert(isfile(rawPath), 'Production raw results CSV is missing; do not overwrite past records.');
    backupTime = char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
    rawBackup = fullfile(outputFolder,['R2_C21_optimized_raw_pre_resume_',backupTime,'.csv']);
    [ok,msg]=copyfile(rawPath,rawBackup);
    assert(ok,'Unable to back up raw CSV: %s',msg);
end
if isfile(rawPath)
    records=readtable(rawPath,'TextType','string');
else
    records=table();
end

fprintf('\nR2 C21 optimized LS: %s mode; output: %s\n',char(mode),outputFolder);
for ip=1:numel(phiTargets)
    phiTarget=phiTargets(ip);
    N=round(phiTarget*L^2/(pi*rFinal^2));
    phiReal=N*pi*rFinal^2/L^2;
    for k=1:nSeeds
        caseMode=mode; if mode=="resume",caseMode="production";end
        caseID=sprintf('LSOPT_%s_P%06d_S%03d',char(caseMode),round(phiTarget*1e5),k);
        if mode=="resume" && ~any(string(caseID)==resumeCaseIDs)
            continue;
        end
        checkpointPath=fullfile(outputFolder,[caseID,'.mat']);
        previousElapsed=0;
        if mode=="resume"
            assert(isfile(checkpointPath), ...
                ['MAT checkpoint missing: %s\n', ...
                 'The coordinate CSV alone CANNOT resume LS dynamics.'],checkpointPath);
            cp=load(checkpointPath);
            needed={'x','v','r','t','L','N','gamma','numCollisions', ...
                'numRebuilds','numDiscarded','elapsed','seed','status'};
            assert(all(isfield(cp,needed)), ...
                'Checkpoint lacks restart variables (cannot safely resume).');
            assert(cp.N==N && abs(cp.L-L)<1e-12 && abs(cp.gamma-gamma)<1e-14, ...
                'Checkpoint geometry or growth rate differs from current settings.');
            assert(all(size(cp.x)==[N 2]) && all(size(cp.v)==[N 2]), ...
                'Checkpoint position/velocity dimensions are invalid.');
            assert(all(isfinite(cp.x),'all') && all(isfinite(cp.v),'all') && ...
                isfinite(cp.r) && isfinite(cp.t), 'Nonfinite checkpoint state.');
            assert(cp.r>=0 && cp.r<rFinal-1e-12, ...
                'Checkpoint already reached the prescribed final radius.');
            assert(startsWith(string(cp.status),'CENSORED'), ...
                'Only computationally censored checkpoints can be resumed.');
            [checkpointGap,checkpointOverlaps]=ls_audit(mod(cp.x,L),L,cp.r,geometricTolerance);
            assert(checkpointOverlaps==0, ...
                'Checkpoint has %d overlaps (minimum gap %.8g).', ...
                checkpointOverlaps,checkpointGap);
            assert(maxCollisionEvents>cp.numCollisions, ...
                'New cumulative event budget must exceed checkpoint collision count.');
            % Preserve original state before overwriting the same case.
            stamp=char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
            backup=fullfile(outputFolder,[char(caseID),'_pre_resume_',stamp,'.mat']);
            [ok,msg]=copyfile(checkpointPath,backup);
            assert(ok,'Unable to back up checkpoint: %s',msg);
            x=mod(cp.x,L);v=cp.v;r=cp.r;t=cp.t;
            numCollisions=cp.numCollisions;
            numRebuilds=cp.numRebuilds;numDiscarded=cp.numDiscarded;
            seed=cp.seed;previousElapsed=cp.elapsed;
            firstBadPair=[NaN NaN];status="UNKNOWN";
            fprintf('Restored checkpoint t=%.9g r=%.9g events=%d, previous runtime=%.1fs\n', ...
                t,r,numCollisions,previousElapsed);
        else
            if ~rerunCompleted && ~isempty(records) && ...
                    any(string(records.CaseID)==string(caseID))
                fprintf('Existing case skipped: %s\n',caseID);continue;
            end
            seed=masterSeed+10000*ip+k+1000000*(mode=="smoke");
            rng(seed,'twister');
            x=L*rand(N,2);
            v=randn(N,2);v=v-mean(v,1);
            v=initialSpeed*v/sqrt(mean(sum(v.^2,2)));
            r=0;t=0;numCollisions=0;numRebuilds=0;numDiscarded=0;
            firstBadPair=[NaN NaN];status="UNKNOWN";
        end
        clear lastTime; % No stale synchronization times from an earlier case.
        startWall=tic;
        fprintf('\n%s | N=%d, target=%.8f\n',caseID,N,phiReal);
        
        while r<rFinal-1e-12
            if toc(startWall)>maxWallSeconds
                status="CENSORED_WALL_BUDGET";break;
            end
            if numCollisions>=maxCollisionEvents
                status="CENSORED_EVENT_BUDGET";break;
            end
            % Bound image movement over the current horizon. Pairwise
            % periodic images -1:1 are included in event predictions.
            vmax=max(hypot(v(:,1),v(:,2)));
            dtH=min([maxHorizon,(rFinal-r)/gamma, L/(8*max(vmax,1e-12))]);
            if ~(isfinite(dtH)&&dtH>0), status="NUMERICAL_FAILURE";break;end
            tEnd=t+dtH;
            vers=zeros(N,1,'uint32'); lastTime=t*ones(N,1);
            H=ls_heap_create(max(256,ceil(3*N)));
            % Initial horizon scan: predict all pair events once, then
            % change only pairs involving collision participants.
            for i=1:N-1
                jj=(i+1:N)';
                [tau,nn]=ls_predict_pairs(x(i,:),v(i,:),x(jj,:),v(jj,:),r,gamma,L,dtH,contactEpsilon);
                keep=find(isfinite(tau));
                for z=1:numel(keep)
                    j=jj(keep(z));
                    H=ls_heap_push(H,t+tau(keep(z)),i,j,vers(i),vers(j),nn(keep(z),:));
                end
            end
            numRebuilds=numRebuilds+1;
            failHorizon=false;
            while H.count>0
                if toc(startWall)>maxWallSeconds || numCollisions>=maxCollisionEvents
                    break;
                end
                [H,ev]=ls_heap_pop(H);
                i=ev.i;j=ev.j;
                if ev.vi~=vers(i) || ev.vj~=vers(j)
                    numDiscarded=numDiscarded+1;continue;
                end
                if ev.time > tEnd+1e-12,break;end
                dt=ev.time-t;
                if dt < -1e-10, status="NUMERICAL_FAILURE";failHorizon=true;break;end
                % Lazily advance only the two colliding disks.
                x(i,:)=mod(x(i,:)+(ev.time-lastTime(i))*v(i,:),L);
                x(j,:)=mod(x(j,:)+(ev.time-lastTime(j))*v(j,:),L);
                lastTime([i j])=ev.time;
                rEvent=r+gamma*max(0,dt);
                t=ev.time;r=rEvent;
                % Recalculate normal from synchronized colliding centers.
                % Normal stored at scheduling can be wrong at an image seam.
                dp=x(i,:)-x(j,:);
                dp=dp-L*round(dp/L);
                separation=hypot(dp(1),dp(2));
                if separation < 1e-13
                    status="NUMERICAL_FAILURE_COINCIDENT_CENTERS";
                    failHorizon=true;break;
                end
                n=dp/separation;
                if separation < 2*r-10*geometricTolerance || ...
                        separation > 2*r+1e-6
                    status="NUMERICAL_FAILURE_EVENT_GEOMETRY";
                    fprintf('Invalid scheduled collision: dist=%g, diameter=%g, t=%g\n', ...
                        separation,2*r,t);
                    failHorizon=true;break;
                end
                % Analytic equal-mass collision for growing hard disks.
                % Surface approach velocity before = dot(v_i-v_j,n)-2 gamma.
                wn=dot(v(i,:)-v(j,:),n);
                if wn >= 2*gamma+2e-8
                    status="NUMERICAL_FAILURE";failHorizon=true;break;
                end
                impulse=(2*gamma-wn)*n;
                v(i,:)=v(i,:)+impulse;
                v(j,:)=v(j,:)-impulse;
                vers(i)=vers(i)+1;vers(j)=vers(j)+1;
                numCollisions=numCollisions+1;
                if perEventAudit
                    % Diagnostic only: synchronize a copy, never modify
                    % the actual lazy event-driven state during auditing.
                    xx=mod(x+(t-lastTime).*v,L);
                    [gEvent,nEvent]=ls_audit(xx,L,r,geometricTolerance);
                    if nEvent>0
                        status="FAILED_EVENT_OVERLAP_AUDIT";
                        fprintf('First bad collision event %d at t=%.12g, r=%.12g, minGap=%.9g, badPairs=%d\n', ...
                            numCollisions,t,r,gEvent,nEvent);
                        failHorizon=true;break;
                    end
                end
                % Schedule events affected by changed velocities, using
                % each particle's lazily propagated position at time t.
                remaining=tEnd-t;
                if remaining > 1e-12
                    participants=[i j];
                    for ii=1:2
                        a=participants(ii);
                        others=(1:N)';others(others==a)=[];
                        xb=mod(x(others,:)+(t-lastTime(others)).*v(others,:),L);
                        [tau,nn]=ls_predict_pairs(x(a,:),v(a,:),xb,v(others,:),r,gamma,L,remaining,contactEpsilon);
                        keep=find(isfinite(tau));
                        for z=1:numel(keep)
                            b=others(keep(z));
                            if ii==2 && b==i,continue;end
                            H=ls_heap_push(H,t+tau(keep(z)),a,b,vers(a),vers(b),nn(keep(z),:));
                        end
                    end
                end
            end
            if failHorizon,break;end
            if toc(startWall)>maxWallSeconds
                status="CENSORED_WALL_BUDGET";break;
            end
            if numCollisions>=maxCollisionEvents
                status="CENSORED_EVENT_BUDGET";break;
            end
            % Advance *all* disks to the horizon boundary; the horizon is
            % merely a safe periodic-image rebuild, not a soft integration.
            x=mod(x+(tEnd-lastTime).*v,L);
            r=r+gamma*(tEnd-t);
            t=tEnd;
            lastTime(:)=tEnd; % prevent double propagation in final save
            if any(~isfinite(x),'all') || any(~isfinite(v),'all') || ~isfinite(r)
                status="NUMERICAL_FAILURE";break;
            end
            % Audit after each horizon to catch missed collisions promptly.
            [gap,nbad]=ls_audit(x,L,r,geometricTolerance);
            if nbad>0
                status="FAILED_INTERMEDIATE_OVERLAP_AUDIT";
                fprintf('Intermediate audit FAILED: %d overlaps, minGap=%g\n',nbad,gap);
                break;
            end
            if mod(numRebuilds,25)==0
                fprintf('  phi %.4f/%.4f | coll %d | rebuild %d | %.1fs\n', ...
                    N*pi*r^2/L^2,phiReal,numCollisions,numRebuilds,toc(startWall));
            end
        end
        % Ensure the saved configuration is synchronized at one physical time.
        % When a wall/event budget terminates the event loop midway through
        % a horizon, untouched disks are still stored at their last event
        % times. Saving them without advancing would make the final audit
        % and coordinate file inconsistent.
        if exist('lastTime','var') && numel(lastTime)==N
            x=mod(x + (t-lastTime).*v,L);
            lastTime(:)=t;
        end
        if status=="UNKNOWN",status="REACHED_RADIUS";end
        if r>=rFinal-1e-12 && status=="REACHED_RADIUS",r=rFinal;end
        segmentElapsed=toc(startWall);
        elapsed=previousElapsed+segmentElapsed;
        [minGap,nBad]=ls_audit(x,L,r,geometricTolerance);
        achieved=N*pi*r*r/L^2;
        success=status=="REACHED_RADIUS" && r==rFinal && nBad==0 && ...
            abs(achieved-phiReal)<1e-9;
        if status=="REACHED_RADIUS" && ~success,status="FAILED_FINAL_AUDIT";end
        fprintf('  %s | phi=%.7f | collisions=%d | rebuilds=%d | wall=%.2fs | overlaps=%d\n', ...
            status,achieved,numCollisions,numRebuilds,elapsed,nBad);
        if writeCoordinates
            writetable(table(x(:,1),x(:,2),r*ones(N,1), ...
                'VariableNames',{'x','y','radius'}), ...
                fullfile(outputFolder,[caseID,'_coords.csv']));
        end
        if saveMAT
            save(fullfile(outputFolder,[caseID,'.mat']), ...
                'x','v','r','L','N','gamma','t','status','numCollisions', ...
                'numRebuilds','numDiscarded','elapsed','seed','-v7.3');
        end
        recordMode=mode; if mode=="resume",recordMode="production";end
        row=table(string(caseID),string(recordMode),phiTarget,phiReal,N,seed, ...
            string(status),logical(success),achieved,r,elapsed,numCollisions, ...
            numRebuilds,numDiscarded,minGap,nBad,t, ...
            'VariableNames',{'CaseID','Mode','TargetPhi','RealizablePhi','N','Seed', ...
            'Status','ExactTargetSuccess','AchievedPhi','FinalRadius','Runtime_s', ...
            'Collisions','EventRebuilds','DiscardedStaleEvents','MinPeriodicGap', ...
            'OverlapsBeyondTol','PhysicalTime'});
        if ~isempty(records)
            records(string(records.CaseID)==string(caseID),:)=[];
            records=[records;row]; %#ok<AGROW>
        else
            records=row;
        end
        writetable(records,rawPath);
    end
end

summaryMode=mode; if mode=="resume",summaryMode="production";end
T=records(string(records.Mode)==summaryMode,:);
if ~isempty(T)
    uniquePhi=unique(T.TargetPhi);out=table();
    for u=1:numel(uniquePhi)
        S=T(abs(T.TargetPhi-uniquePhi(u))<1e-10,:);
        successMask=logical(S.ExactTargetSuccess);
        succ=S(successMask,:);
        m=NaN;if ~isempty(succ),m=mean(succ.Runtime_s);end
        row=table(uniquePhi(u),S.RealizablePhi(1),height(S), ...
            sum(successMask),sum(startsWith(string(S.Status),"CENSORED")), ...
            m,'VariableNames',{'TargetPhi','RealizablePhi','Runs', ...
            'ExactSuccesses','Censored','MeanSuccessfulRuntime_s'});
        out=[out;row]; %#ok<AGROW>
    end
    writetable(out,fullfile(outputFolder,'R2_C21_optimized_summary.csv'));
    disp(out);
end
fprintf('\nDone. Only the requested resume cases were executed.\n');
fprintf('Interpret wall/event-budget stops as censored, not LS packing failures.\n');

%% =================== LOCAL FUNCTIONS ================================
function [tau,nOut]=ls_predict_pairs(xa,va,xb,vb,r,gamma,L,horizon,epsT)
% Vectorized pair prediction for one particle a against many b, considering
% periodic images. Gives earliest approaching contact within horizon.
M=size(xb,1);tau=inf(M,1);nOut=zeros(M,2);
if M==0,return;end
w=va-vb; p0=xa-xb; R=2*r;G=2*gamma;
A=sum(w.^2,2)-G^2;
for sx=-1:1
    for sy=-1:1
        p=p0+L*[sx sy];
        B=2*(sum(p.*w,2)-R*G);
        C=sum(p.^2,2)-R^2;
        % Solve |p + w*t| = R + G*t using BOTH quadratic roots.
        % A smaller *negative* root must not hide a positive candidate.
        % Candidate is accepted only for an approaching contact.
        discr=B.^2-4*A.*C;
        cand=inf(M,1);
        lin=abs(A)<1e-13 & abs(B)>1e-13;
        tl=-C(lin)./B(lin);
        goodL=tl>=-epsT;
        idsL=find(lin);
        cand(idsL(goodL))=max(0,tl(goodL));
        quad=abs(A)>=1e-13 & discr>=0;
        idsQ=find(quad);
        if ~isempty(idsQ)
            sq=sqrt(max(0,discr(idsQ)));
            aq=A(idsQ);bq=B(idsQ);cq=C(idsQ);
            sg=ones(size(bq));sg(bq<0)=-1;
            q=-.5*(bq+sg.*sq);
            t1=q./aq;
            t2=inf(size(t1));
            stable=abs(q)>1e-30;
            t2(stable)=cq(stable)./q(stable);
            t2(~stable)=(-bq(~stable)+sq(~stable))./(2*aq(~stable));
            t1(t1 < -epsT)=inf;
            t2(t2 < -epsT)=inf;
            cand(idsQ)=max(0,min(t1,t2));
        end
        % Contact at current time with closing surface-normal speed.
        immediate=abs(C)<1e-10 & B<-1e-13;
        cand(immediate)=0;
        ok=isfinite(cand) & cand>=0 & cand<=horizon+1e-12 & cand<tau;
        ids=find(ok);
        if isempty(ids),continue;end
        pc=p(ids,:)+cand(ids).*w(ids,:);
        d=sqrt(sum(pc.^2,2));
        nn=pc./max(d,realmin);
        approach=sum(nn.*w(ids,:),2) < G-1e-11;
        ids=ids(approach);nn=nn(approach,:);
        tau(ids)=cand(ids);nOut(ids,:)=nn;
    end
end
end

function H=ls_heap_create(cap)
H.t=zeros(cap,1);H.i=zeros(cap,1,'uint32');H.j=H.i;
H.vi=H.i;H.vj=H.i;H.nx=zeros(cap,1);H.ny=zeros(cap,1);
H.count=0;
end

function H=ls_heap_push(H,t,i,j,vi,vj,n)
if ~isfinite(t),return;end
if H.count>=numel(H.t)
    % Geometric array growth amortizes allocation cost.
    nNew=max(256,2*numel(H.t));
    H.t(nNew)=0;H.i(nNew)=uint32(0);H.j(nNew)=uint32(0);
    H.vi(nNew)=uint32(0);H.vj(nNew)=uint32(0);
    H.nx(nNew)=0;H.ny(nNew)=0;
end
k=H.count+1;H.count=k;
while k>1
    parent=floor(k/2);
    if H.t(parent)<=t,break;end
    H=ls_heap_copy(H,k,parent);k=parent;
end
H.t(k)=t;H.i(k)=uint32(i);H.j(k)=uint32(j);
H.vi(k)=vi;H.vj(k)=vj;H.nx(k)=n(1);H.ny(k)=n(2);
end

function [H,e]=ls_heap_pop(H)
assert(H.count>0,'Cannot pop empty heap');
e.time=H.t(1);e.i=double(H.i(1));e.j=double(H.j(1));
e.vi=H.vi(1);e.vj=H.vj(1);e.normal=[H.nx(1) H.ny(1)];
last=H.count;H.count=last-1;
if last==1,return;end
v=H.t(last);ii=H.i(last);jj=H.j(last);a=H.vi(last);b=H.vj(last);
nx=H.nx(last);ny=H.ny(last);
k=1;
while 2*k<=H.count
    ch=2*k;
    if ch+1<=H.count && H.t(ch+1)<H.t(ch),ch=ch+1;end
    if H.t(ch)>=v,break;end
    H=ls_heap_copy(H,k,ch);k=ch;
end
H.t(k)=v;H.i(k)=ii;H.j(k)=jj;H.vi(k)=a;H.vj(k)=b;
H.nx(k)=nx;H.ny(k)=ny;
end

function H=ls_heap_copy(H,d,s)
H.t(d)=H.t(s);H.i(d)=H.i(s);H.j(d)=H.j(s);
H.vi(d)=H.vi(s);H.vj(d)=H.vj(s);
H.nx(d)=H.nx(s);H.ny(d)=H.ny(s);
end

function [gap,nBad]=ls_audit(x,L,r,tol)
N=size(x,1);gap=inf;nBad=0;
for i=1:N-1
    dp=abs(x(i+1:end,:)-x(i,:));dp=min(dp,L-dp);
    gg=hypot(dp(:,1),dp(:,2))-2*r;
    gap=min(gap,min(gg));nBad=nBad+sum(gg < -tol);
end
end

function ls_opt_selftest()
% Analytic collision: p=(1,0), w=(-1,0), 2r=.2, 2gamma=.1 => t=8/11.
[tt,nn]=ls_predict_pairs([1 1],[-.5 0],[0 1],[.5 0],.1,.05,10,1,1e-12);
assert(abs(tt-8/11)<1e-8 && abs(norm(nn)-1)<1e-9,'Collision root regression failed');
% Regression: earlier root is negative, later root is positive.
% The positive root must be available to the contact-direction filter.
[tp,~]=ls_predict_pairs([1 1],[1 0],[0 1],[0 0],.1,0,10,5,1e-12);
assert(isinf(tp),'Separating particles produced spurious collision');
% Periodic cross-boundary contact.
[t,n]=ls_predict_pairs([.1 1],[-.5 0],[9.6 1],[0 0],.1,0,10,2,1e-12);
assert(isfinite(t)&&abs(t-.6)<1e-8&&abs(norm(n)-1)<1e-8,'Periodic collision regression failed');
[g,bad]=ls_audit([1 1;1.4 1],10,.2,1e-8);
assert(abs(g)<1e-12 && bad==0,'Overlap audit regression failed');
% Heap ordering, including tied timestamps.
H=ls_heap_create(4);
H=ls_heap_push(H,3,1,2,uint32(0),uint32(0),[1 0]);
H=ls_heap_push(H,1,1,3,uint32(0),uint32(0),[1 0]);
H=ls_heap_push(H,2,2,3,uint32(0),uint32(0),[1 0]);
for target=1:3
    [H,e]=ls_heap_pop(H);assert(abs(e.time-target)<1e-12,'Heap ordering regression failed');
end
% At zero growth, an equal-mass head-on collision swaps normal speeds.
na=[1 0];va=[-0.5 0];vb=[0.5 0];g=0;
w=dot(va-vb,na);impulse=(2*g-w)*na;
vap=va+impulse;vbp=vb-impulse;
assert(norm(vap-[.5 0])<1e-12 && norm(vbp-[-.5 0])<1e-12, ...
    'Velocity impulse regression failed');
% Correct quadratic root if negative and positive roots coexist;
% monotonic grow against stationary neighbor must hit at t=0.8.
[tr,~]=ls_predict_pairs([1 1],[0 0],[0 1],[0 0],.1,.5,10,1,1e-12);
assert(abs(tr-.8)<1e-10,'Growing-radius root regression failed');
fprintf('Analytical collision, periodic-image, audit, impulse, and heap tests PASSED.\n');
end


function ls_c21_audit_archive(outputFolder,projectFolder,L,rFinal,tol)
% Read-only audit. Never changes the production raw, summary, coordinates or MATs.
assert(isfolder(outputFolder), ['LS archive folder not found: %s\n', ...
    'Set projectFolder or outputFolder in USER CONFIGURATION.'], outputFolder);
rawFile=fullfile(outputFolder,'R2_C21_optimized_raw.csv');
assert(isfile(rawFile),'Production raw CSV not found: %s',rawFile);
R=readtable(rawFile,'TextType','string');
required={'CaseID','Mode','TargetPhi','RealizablePhi','N','Status', ...
 'ExactTargetSuccess','AchievedPhi','FinalRadius','Runtime_s','Collisions'};
assert(all(ismember(required,R.Properties.VariableNames)), ...
       'Production raw CSV does not contain expected V5 schema.');
R.CaseID=string(R.CaseID);R.Mode=string(R.Mode);R.Status=string(R.Status);
if ~islogical(R.ExactTargetSuccess)
  R.ExactTargetSuccess=ismember(lower(string(R.ExactTargetSuccess)),["true","1"]);
end
% Audit production only; smoke and pilot are excluded from final benchmark.
P=R(R.Mode=="production",:);
assert(~isempty(P),'No production records found.');
% Refuse duplicate cases; the raw file should contain one latest state per case.
assert(numel(unique(P.CaseID))==height(P),'Duplicate production CaseID rows detected.');
runFolder=fullfile(projectFolder,'revision_results','R2_C21_consolidated_audit');
if ~isfolder(runFolder),mkdir(runFolder);end
stamp=char(datetime('now','Format','yyyyMMdd_HHmmss'));
reportFolder=fullfile(runFolder,['audit_',stamp]);
if isfolder(reportFolder),reportFolder=[reportFolder,'_repeat'];end
mkdir(reportFolder);
fid=fopen(fullfile(reportFolder,'audit_manifest.txt'),'w');
assert(fid>0,'Could not create audit manifest.');
fprintf(fid,'R2-C21 consolidated archival audit\nMATLAB: %s\nArchive: %s\n',version,outputFolder);
fprintf(fid,'Raw CSV: %s\nGeometric tolerance: %.17g\n',rawFile,tol);
fprintf(fid,'Audits do not rerun dynamics or alter source archives.\n');
fprintf(fid,'Censored trajectories are not geometric failures.\n');
fclose(fid);
Nrec=height(P);
Geom_CaseID=P.CaseID;Geom_Status=P.Status;
Geom_File=strings(Nrec,1);Geom_Particles=nan(Nrec,1);
Geom_AreaFraction=nan(Nrec,1);Geom_MinGap=nan(Nrec,1);
Geom_OverlapPairs=nan(Nrec,1);Geom_Radius=nan(Nrec,1);
Geom_Pass=false(Nrec,1);Geom_AuditNote=strings(Nrec,1);
for k=1:Nrec
    name=char(P.CaseID(k));
    csv=fullfile(outputFolder,[name,'_coords.csv']);
    if ~isfile(csv)
        Geom_AuditNote(k)="Coordinate file missing; record not independently validated";
        continue;
    end
    Geom_File(k)=string(csv);
    C=readtable(csv,'VariableNamingRule','preserve');
    columns=lower(string(C.Properties.VariableNames));
    ix=find(columns=="x",1);iy=find(columns=="y",1);
    ir=find(ismember(columns,["radius","r"]),1);
    if isempty(ix)||isempty(iy)||isempty(ir)
        Geom_AuditNote(k)="Missing x/y/radius columns";continue;
    end
    x=double(C{:,ix});y=double(C{:,iy});rs=double(C{:,ir});
    n=numel(x);Geom_Particles(k)=n;
    if n~=P.N(k)||any(~isfinite([x;y;rs]))||any(rs<=0)|| ...
       any(x<0|x>=L|y<0|y>=L)
        Geom_AuditNote(k)="Invalid counts, coordinates or radii";continue;
    end
    Geom_Radius(k)=mean(rs);
    Geom_AreaFraction(k)=sum(pi*rs.^2)/L^2;
    worst=inf;bad=0;
    for i=1:n-1
        dx=abs(x(i+1:end)-x(i));dx=min(dx,L-dx);
        dy=abs(y(i+1:end)-y(i));dy=min(dy,L-dy);
        gap=hypot(dx,dy)-(rs(i)+rs(i+1:end));
        worst=min(worst,min(gap));bad=bad+sum(gap < -tol);
    end
    Geom_MinGap(k)=worst;Geom_OverlapPairs(k)=bad;
    sameRadii=all(abs(rs-rs(1))<1e-10);
    finalOK=~P.ExactTargetSuccess(k) || ...
        (abs(rs(1)-rFinal)<=1e-10 && ...
         abs(Geom_AreaFraction(k)-P.RealizablePhi(k))<1e-8);
    Geom_Pass(k)=bad==0 && sameRadii && finalOK && ...
        abs(Geom_AreaFraction(k)-P.AchievedPhi(k))<1e-8;
    if Geom_Pass(k),Geom_AuditNote(k)="PASS";
    else,Geom_AuditNote(k)="Mismatch/overlap; inspect saved coordinates";end
end
G=table(Geom_CaseID,Geom_Status,Geom_File,Geom_Particles, ...
 Geom_AreaFraction,Geom_Radius,Geom_MinGap,Geom_OverlapPairs, ...
 Geom_Pass,Geom_AuditNote, ...
 'VariableNames',{'CaseID','Status','CoordinateFile','Particles', ...
 'RecomputedPhi','MeanRadius','MinPeriodicGap','PairsBeyondTolerance', ...
 'Pass','AuditNote'});
writetable(G,fullfile(reportFolder,'R2_C21_coordinate_audit.csv'));
phiList=unique(P.TargetPhi);
Summary_TargetPhi=zeros(0,1);Summary_RealizablePhi=zeros(0,1);
Summary_Attempted=zeros(0,1);Summary_Completed=zeros(0,1);
Summary_Censored=zeros(0,1);Summary_MeanCompletedTime_s=zeros(0,1);
Summary_CompletedAudited=zeros(0,1);
for q=1:numel(phiList)
  idx=abs(P.TargetPhi-phiList(q))<1e-10;
  tmp=P(idx,:);idxRows=find(idx);
  done=tmp.ExactTargetSuccess & tmp.Status=="REACHED_RADIUS";
  censored=startsWith(tmp.Status,"CENSORED");
  times=tmp.Runtime_s(done);
  avg=NaN;if ~isempty(times),avg=mean(times);end
  Summary_TargetPhi(end+1,1)=phiList(q);
  Summary_RealizablePhi(end+1,1)=tmp.RealizablePhi(1);
  Summary_Attempted(end+1,1)=height(tmp);
  Summary_Completed(end+1,1)=sum(done);
  Summary_Censored(end+1,1)=sum(censored);
  Summary_MeanCompletedTime_s(end+1,1)=avg;
  Summary_CompletedAudited(end+1,1)=sum(done & Geom_Pass(idxRows));
end
S=table(Summary_TargetPhi,Summary_RealizablePhi,Summary_Attempted, ...
 Summary_Completed,Summary_Censored,Summary_MeanCompletedTime_s, ...
 Summary_CompletedAudited, ...
 'VariableNames',{'TargetPhi','RealizablePhi','Attempted','Completed', ...
 'ComputationallyCensored','MeanCompletedTime_s','CompletedWithPassingCoordinateAudit'});
writetable(S,fullfile(reportFolder,'R2_C21_recomputed_summary.csv'));
% Existing VGPP comparison is included only if a genuine archived table is found.
comparison=fullfile(projectFolder,'revision_results','R2_C21_VGPP_LS_final_comparison.csv');
if ~isfile(comparison)
  comparison=fullfile(projectFolder,'R2_C21_VGPP_LS_final_comparison.csv');
end
if isfile(comparison)
    copyfile(comparison,fullfile(reportFolder,'archived_VGPP_LS_comparison_unmodified.csv'));
end
fig=figure('Visible','off','Color','w','Position',[100 100 880 520]);
bar(categorical(string(S.TargetPhi)),[S.Completed,S.ComputationallyCensored],'stacked');
ylabel('Number of LS realizations');xlabel('Nominal particle area fraction');
legend({'Reached final radius','Computationally censored'},'Location','best');
title('LS production outcomes (archived runs)');
set(gca,'FontName','Times New Roman','FontSize',13);
exportgraphics(fig,fullfile(reportFolder,'R2_C21_completion_summary.png'),'Resolution',250);
close(fig);
report=fopen(fullfile(reportFolder,'R2_C21_interpretation.txt'),'w');
fprintf(report,'No dynamics were executed. All outputs were derived from archived results.\n');
fprintf(report,'Completed coordinate audits passed: %d/%d available completion records.\n', ...
    sum(Geom_Pass & P.ExactTargetSuccess),sum(P.ExactTargetSuccess));
fprintf(report,'A missing coordinate file is not proof of an invalid LS configuration.\n');
fprintf(report,'Only completed runtimes enter reported means; censored run costs are excluded.\n');
fprintf(report,'Unequal VGPP/LS realization counts and implementation details prevent general speed ranking.\n');
for q=1:height(S)
 fprintf(report,'phi=%.2f completed=%d attempted=%d censored=%d meanCompletedTime=%.6g s\n', ...
 S.TargetPhi(q),S.Completed(q),S.Attempted(q),S.ComputationallyCensored(q),S.MeanCompletedTime_s(q));
end
fclose(report);
disp(S);
fprintf('Read-only R2-C21 audit saved to: %s\n',reportFolder);
end
