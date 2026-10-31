function [row, particleTable, historicalParticleTable] = benchmark_case_checkpoint(params, method, targetVF, nTarget, it, irun)
%BENCHMARK_CASE_CHECKPOINT Resume an identical recorded case or run and save it.
if ~isfield(params,'enableCheckpoints') || ~params.enableCheckpoints
    [row,particleTable,historicalParticleTable] = run_benchmark_case(params,method,targetVF,nTarget,it,irun);
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
if isfile(filename) && (~isfield(params,'resumeCheckpoints') || ~params.resumeCheckpoints)
    error('VGPP:CheckpointMismatch','Existing checkpoint preserved: %s.',filename);
end
if isfile(filename) && isfield(params,'resumeCheckpoints') && params.resumeCheckpoints
    cached = load(filename,'row','particleTable','configuration','formatVersion','sourceFingerprint','historicalParticleTable','statisticalConvention');
    if isfield(cached,'formatVersion') && isequal(cached.formatVersion,3) && ...
            isfield(cached,'historicalParticleTable') && ...
            isfield(cached,'statisticalConvention') && ...
            strcmp(cached.statisticalConvention,'historical-boundary-v1') && ...
            isfield(cached,'sourceFingerprint') && ...
            strcmp(cached.sourceFingerprint,sourceFingerprint) && ...
            isfield(cached,'configuration') && isequaln(cached.configuration,configuration) && ...
            isfield(cached,'row') && height(cached.row)==1 && ...
            cached.row.Method==method && ...
            abs(cached.row.TargetVF-targetVF)<1e-12 && ...
            cached.row.TargetN==nTarget && cached.row.Run==irun && ...
            cached.row.Seed==7300+100*it+irun
        validate_historical_coordinates(cached.particleTable,cached.historicalParticleTable,params);
        historicalParticleTable = cached.historicalParticleTable;
        row = cached.row;
        particleTable = cached.particleTable;
        fprintf('  resumed %s at %.2f, run %d\n',char(method),targetVF,irun);
        return;
    end
    error('VGPP:CheckpointMismatch', ...
        'Checkpoint parameters or source checksum differ: %s. Preserve old checkpoints and use a compatible output directory.',filename);
end
[row,particleTable,historicalParticleTable] = run_benchmark_case(params,method,targetVF,nTarget,it,irun);
formatVersion = 3;
statisticalConvention = 'historical-boundary-v1';
% Save within the checkpoint directory, then rename after a successful write.
tmp = [tempname(folder) '.mat'];
try
    save(tmp,'row','particleTable','configuration','formatVersion','sourceFingerprint','historicalParticleTable','statisticalConvention','-v7.3');
    [ok,msg] = movefile(tmp,filename,'f');
    if ~ok, error('VGPP:CheckpointWrite','%s',msg); end
catch ME
    if isfile(tmp), delete(tmp); end
    rethrow(ME);
end
end
