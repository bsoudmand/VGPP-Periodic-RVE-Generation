function [row, particleTable] = benchmark_case_checkpoint(params, method, targetVF, nTarget, it, irun)
%BENCHMARK_CASE_CHECKPOINT Resume an identical recorded case or run and save it.
if ~isfield(params,'enableCheckpoints') || ~params.enableCheckpoints
    [row,particleTable] = run_benchmark_case(params,method,targetVF,nTarget,it,irun);
    return;
end
folder = params.checkpointFolder;
if ~isfolder(folder), mkdir(folder); end
filename = fullfile(folder,sprintf('case_%s_af%03d_run%03d.mat', ...
    char(method),round(1000*targetVF),irun));
% A checkpoint depends on both configuration and the exact active source tree.
if ~isfield(params,'sourceFingerprint') || strlength(string(params.sourceFingerprint)) ~= 64
    error('VGPP:SourceFingerprint', 'Source fingerprint missing: run from main.m.');
end
sourceFingerprint = char(params.sourceFingerprint);
% Retain complete configuration to reject inconsistent resumed experiments.
configuration = params;
if isfile(filename) && isfield(params,'resumeCheckpoints') && params.resumeCheckpoints
    cached = load(filename,'row','particleTable','configuration','formatVersion','sourceFingerprint');
    if isfield(cached,'formatVersion') && isequal(cached.formatVersion,2) && ...
            isfield(cached,'sourceFingerprint') && ...
            strcmp(cached.sourceFingerprint,sourceFingerprint) && ...
            isfield(cached,'configuration') && isequaln(cached.configuration,configuration) && ...
            isfield(cached,'row') && height(cached.row)==1 && ...
            cached.row.Method==method && ...
            abs(cached.row.TargetVF-targetVF)<1e-12 && ...
            cached.row.TargetN==nTarget && cached.row.Run==irun && ...
            cached.row.Seed==7300+100*it+irun
        row = cached.row;
        particleTable = cached.particleTable;
        fprintf('  resumed %s at %.2f, run %d\n',char(method),targetVF,irun);
        return;
    end
    error('VGPP:CheckpointMismatch', ...
        'Checkpoint parameters or source checksum differ: %s. Use a separate output directory or clear old checkpoints.',filename);
end
[row,particleTable] = run_benchmark_case(params,method,targetVF,nTarget,it,irun);
formatVersion = 2;
% Save within the checkpoint directory, then rename after a successful write.
tmp = [tempname(folder) '.mat'];
try
    save(tmp,'row','particleTable','configuration','formatVersion','sourceFingerprint','-v7.3');
    [ok,msg] = movefile(tmp,filename,'f');
    if ~ok, error('VGPP:CheckpointWrite','%s',msg); end
catch ME
    if isfile(tmp), delete(tmp); end
    rethrow(ME);
end
end
