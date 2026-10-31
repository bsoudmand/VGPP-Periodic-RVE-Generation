function report=vgpp_release_checks(repoRoot,outDir)
%VGPP_RELEASE_CHECKS Targeted regressions; not the full benchmark validation.
if nargin<1,repoRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));end
if nargin<2,outDir=fullfile(repoRoot,'outputs','validation');end
if ~isfolder(outDir),mkdir(outDir);end
previousPath=path;restorePath=onCleanup(@() path(previousPath)); %#ok<NASGU>
addpath(fullfile(repoRoot,'config'),fullfile(repoRoot,'src','core'), ...
    fullfile(repoRoot,'src','benchmark'),fullfile(repoRoot,'src','analysis'), ...
    fullfile(repoRoot,'studies','optional'));
p=vgpp_config(repoRoot);p.sourceFingerprint=vgpp_source_fingerprint(repoRoot);
vgpp_self_test(p);
tol=p.overlapAuditTolerance;
[bad,~,~,n]=final_overlap_audit_PBC([1;1+.36-.5*tol],[1;1],[.18;.18],p);
assert(~bad && n==0);
[bad,~,~,n]=final_overlap_audit_PBC([1;1+.36-1.5*tol],[1;1],[.18;.18],p);
assert(bad && n==1);
for radius={NaN,Inf,0,-.18}
    rejectGeometry(@() final_overlap_audit_PBC([1;5],[1;5],[radius{1};.18],p));
end
rejectGeometry(@() final_overlap_audit_PBC([1;p.Lx],[1;5],[.18;.18],p));
rejectGeometry(@() final_overlap_audit_PBC([1;5],[1;-1],[.18;.18],p));
rejectGeometry(@() final_overlap_audit_PBC([1;5],1,[.18;.18],p));
n=round(.4*p.Lx*p.Ly/(pi*p.rMono^2));
[row,particles]=run_benchmark_case(p,"VGPP",.4,n,1,2);
assert(row.GeometryValid && row.ReachedTarget);
assert(all(particles.x>=0 & particles.x<p.Lx & particles.y>=0 & particles.y<p.Ly));
writetable(row,fullfile(outDir,'endpoint_regression.csv'));
p.checkpointFolder=fullfile(outDir,'checkpoint-regression');
[r1,~]=benchmark_case_checkpoint(p,"VGPP",.4,n,1,2);
file=fullfile(p.checkpointFolder,'case_VGPP_af400_run002.mat');
before=dir(file);[r2,~]=benchmark_case_checkpoint(p,"VGPP",.4,n,1,2);after=dir(file);
assert(isequaln(r1,r2) && before.bytes==after.bytes && before.datenum==after.datenum);
assert(numel(dir(fullfile(p.checkpointFolder,'case_*.mat')))==1);
q=p;q.batchSize=p.batchSize+1;
rejectCheckpoint(@() benchmark_case_checkpoint(q,"VGPP",.4,n,1,2));
q=p;q.sourceFingerprint=repmat('0',1,64);
rejectCheckpoint(@() benchmark_case_checkpoint(q,"VGPP",.4,n,1,2));
original=load(file);changed=original;changed.row.Seed=999;
save(file,'-struct','changed','-v7.3');
guard=onCleanup(@() restoreCheckpoint(file,original));
rejectCheckpoint(@() benchmark_case_checkpoint(p,"VGPP",.4,n,1,2));
clear guard;
% Version-2 records cannot recover boundary-side metadata automatically.
changed=original;changed.formatVersion=2;
save(file,'-struct','changed','-v7.3');
guard=onCleanup(@() restoreCheckpoint(file,original));
rejectCheckpoint(@() benchmark_case_checkpoint(p,"VGPP",.4,n,1,2));
clear guard;
changed=original;changed=rmfield(changed,'historicalParticleTable');
save(file,'-struct','changed','-v7.3');
guard=onCleanup(@() restoreCheckpoint(file,original));
rejectCheckpoint(@() benchmark_case_checkpoint(p,"VGPP",.4,n,1,2));
clear guard;
q=p;q.resumeCheckpoints=false;
rejectCheckpoint(@() benchmark_case_checkpoint(q,"VGPP",.4,n,1,2));
% Upper-boundary storage maps to zero geometry without changing genuine zero.
c=table([1;2],[0;0],[1;2],[.18;.18], ...
    'VariableNames',{'ID','x','y','radius'});h=c;h.x(1)=p.Lx;
validate_historical_coordinates(c,h,p);
h.x(2)=1;
rejectGeometry(@() validate_historical_coordinates(c,h,p));
ls_opt_selftest();
Check=["Core self-test";"Continuous-coordinate tolerance";"Invalid geometry rejection"; ...
    "Upper-endpoint seed 7402 regression";"Checkpoint v3 resume, legacy rejection, metadata and overwrite protection";"Optional LS analytical tests"];
Passed=true(numel(Check),1);report=table(Check,Passed);
writetable(report,fullfile(outDir,'release_checks.csv'));
disp(report);
end

function rejectGeometry(fn)
rejected=false;try,fn();catch ME
rejected=startsWith(string(ME.identifier),"MATLAB:") || startsWith(string(ME.identifier),"VGPP:");
end
assert(rejected,'Invalid geometry was accepted.');
end

function rejectCheckpoint(fn)
rejected=false;try,fn();catch ME,rejected=strcmp(ME.identifier,'VGPP:CheckpointMismatch');end
assert(rejected,'Incompatible checkpoint was accepted.');
end

function restoreCheckpoint(file,original)
save(file,'-struct','original','-v7.3');
end
