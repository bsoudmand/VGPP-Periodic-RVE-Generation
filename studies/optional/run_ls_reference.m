function run_ls_reference(repoRoot, mode)
%RUN_LS_REFERENCE Optional LS numerical reference and archive auditor.
if nargin < 2, mode = "audit"; end
mode = string(mode);
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

%% USER CONFIGURATION
% Modes: audit (default) / smoke / pilot / production / resume
% Consolidated: default READ-ONLY audit; optional smoke/pilot/production/resume.
% Existing LS dynamics are retained from validated V5 continuation engine.
% Resume reads an existing MAT checkpoint. It NEVER reconstructs velocity
% information from a coordinate CSV, and it never silently restarts a case.
% Mode is supplied by the caller.
resumeCaseIDs = "LSOPT_production_P074000_S002"; % used ONLY if mode="resume"
rerunCompleted = false;
projectFolder = repoRoot;
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
end
