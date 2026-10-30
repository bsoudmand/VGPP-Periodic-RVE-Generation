%% R2_C21_AdvancedPackingBenchmark_Optimized_FIXED.m
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
mode = "smoke"; % "smoke" (fast verification), "pilot" (.60, 1 seed), "production" (.60,.70,.74 x 3)
rerunCompleted = false;
outputFolder = 'revision_results\R2_C21_LS_optimized';
L=10; rFinal=.18; masterSeed=210521;
gamma=.002; initialSpeed=1;                 % radial growth speed, initial mean speed
maxWallSeconds=1800; maxCollisionEvents=300000;
maxHorizon=0.50;                          % each rebuild covers finite free-flight interval
geometricTolerance=1e-8; contactEpsilon=2e-11;
writeCoordinates=true; saveMAT=true;

switch mode
    case "smoke"
        phiTargets=.10; nSeeds=1; gamma=.01; maxWallSeconds=90;
        maxCollisionEvents=20000;
    case "pilot"
        phiTargets=.60; nSeeds=1;
    case "production"
        phiTargets=[.60 .70 .74]; nSeeds=3;
    otherwise
        error('Unknown mode. Use smoke, pilot, or production.');
end

if ~isfolder(outputFolder), mkdir(outputFolder); end
settings=fullfile(outputFolder,'R2_C21_optimized_settings.txt');
fid=fopen(settings,'w'); assert(fid>0,'Could not write settings file.');
fprintf(fid,'MATLAB=%s\nmode=%s\nL=%.12g\nrFinal=%.12g\nphiTargets=%s\nnSeeds=%d\n', ...
    version,char(mode),L,rFinal,mat2str(phiTargets),nSeeds);
fprintf(fid,'masterSeed=%d\ngrowthRate=%.12g\ninitialSpeed=%.12g\n',masterSeed,gamma,initialSpeed);
fprintf(fid,'maxWallSeconds=%.12g\nmaxCollisionEvents=%d\nmaxHorizon=%.12g\n', ...
    maxWallSeconds,maxCollisionEvents,maxHorizon);
fprintf(fid,'geometricTolerance=%.12g\n',geometricTolerance);
fprintf(fid,'Collision-driven growing hard disks; no force-based overlap correction.\n');
fclose(fid);

% Fail fast if the analytic collision, velocity update, periodic audit or heap
% ordering tests fail. These tests are necessary, but not sufficient, for
% full research-grade LS validation.
ls_opt_selftest();

rawPath=fullfile(outputFolder,'R2_C21_optimized_raw.csv');
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
        caseID=sprintf('LSOPT_%s_P%06d_S%03d',char(mode),round(phiTarget*1e5),k);
        if ~rerunCompleted && ~isempty(records) && any(string(records.CaseID)==string(caseID))
            fprintf('Checkpoint: skipped %s\n',caseID); continue;
        end
        seed=masterSeed+10000*ip+k+1000000*(mode=="smoke");
        rng(seed,'twister');
        x=L*rand(N,2);
        v=randn(N,2);v=v-mean(v,1);
        v=initialSpeed*v/sqrt(mean(sum(v.^2,2)));
        r=0; t=0; numCollisions=0; numRebuilds=0; numDiscarded=0;
        status="UNKNOWN"; startWall=tic;
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
                n=ev.normal;
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
        elapsed=toc(startWall);
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
        row=table(string(caseID),string(mode),phiTarget,phiReal,N,seed, ...
            string(status),logical(success),achieved,r,elapsed,numCollisions, ...
            numRebuilds,numDiscarded,minGap,nBad,t, ...
            'VariableNames',{'CaseID','Mode','TargetPhi','RealizablePhi','N','Seed', ...
            'Status','ExactTargetSuccess','AchievedPhi','FinalRadius','Runtime_s', ...
            'Collisions','EventRebuilds','DiscardedStaleEvents','MinPeriodicGap', ...
            'OverlapsBeyondTol','PhysicalTime'});
        if ~isempty(records)
            records(records.CaseID==string(caseID),:)=[];
            records=[records;row]; %#ok<AGROW>
        else
            records=row;
        end
        writetable(records,rawPath);
    end
end

T=records(records.Mode==mode,:);
if ~isempty(T)
    uniquePhi=unique(T.TargetPhi);out=table();
    for u=1:numel(uniquePhi)
        S=T(abs(T.TargetPhi-uniquePhi(u))<1e-10,:);
        succ=S(S.ExactTargetSuccess,:);
        m=NaN;if ~isempty(succ),m=mean(succ.Runtime_s);end
        row=table(uniquePhi(u),S.RealizablePhi(1),height(S), ...
            sum(S.ExactTargetSuccess),sum(startsWith(S.Status,"CENSORED")), ...
            m,'VariableNames',{'TargetPhi','RealizablePhi','Runs', ...
            'ExactSuccesses','Censored','MeanSuccessfulRuntime_s'});
        out=[out;row]; %#ok<AGROW>
    end
    writetable(out,fullfile(outputFolder,'R2_C21_optimized_summary.csv'));
    disp(out);
end
fprintf('\nDone. Interpret wall/event-budget stops as censored, not LS packing failures.\n');

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
        % Roundoff at immediate contacts is allowed; deep overlap is not.
        C(abs(C)<1e-9)=max(0,C(abs(C)<1e-9));
        discr=B.^2-4*A.*C;
        cand=inf(M,1);
        lin=abs(A)<1e-13 & B < 0;
        cand(lin)=-C(lin)./B(lin);
        quad=abs(A)>=1e-13 & discr>=0;
        if any(quad)
            inds=find(quad);
            sq=sqrt(max(0,discr(inds)));
            aa=A(inds);bb=B(inds);
            q=(-bb-sign(bb+epsT).*sq)/2;
            t1=q./aa;t2=C(inds)./q;
            z=min(t1,t2);
            z(z<0)=inf;
            tAlt=max(t1,t2);
            useAlt=~isfinite(z);
            z(useAlt)=tAlt(useAlt);
            cand(inds)=z;
        end
        nearZero=C <= 1e-10 & B < -1e-13;
        cand(nearZero)=0;
        ok=isfinite(cand) & cand>=0 & cand<=horizon+1e-12 & cand<tau;
        ids=find(ok);
        if isempty(ids),continue;end
        pc=p(ids,:)+cand(ids).*w(ids,:);
        d=sqrt(sum(pc.^2,2));
        nn=pc./max(d,realmin);
        approach=sum(nn.*w(ids,:),2) < G+1e-10;
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
fprintf('Analytical collision, periodic-image, audit, impulse, and heap tests PASSED.\n');
end
