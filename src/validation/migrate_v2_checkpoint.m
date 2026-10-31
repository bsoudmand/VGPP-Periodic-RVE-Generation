function status=migrate_v2_checkpoint(oldFile,legacyRoot,p)
%MIGRATE_V2_CHECKPOINT Explicit migration of the audited preview source only.
% Never infer a lost upper-boundary side from a canonical zero. No overwrite.
known='80af3a550fdec5744dbfbb3934fe3ed6db2fc532edbfa6f4591f7ef5931cc391';
assert(strcmp(vgpp_source_fingerprint(legacyRoot),known), ...
    'VGPP:MigrationSource','Legacy tree is not the audited source.');
c=load(oldFile);
assert(c.formatVersion==2 && strcmp(c.sourceFingerprint,known), ...
    'VGPP:MigrationSource','Unsupported checkpoint source or version.');
% Only paths and fingerprints may differ; all numerical parameters must match.
old=c.configuration;new=p;
ignore={'resultsFolder','figuresFolder','referenceFile','checkpointFolder','sourceFingerprint'};
for k=1:numel(ignore)
    if isfield(old,ignore{k}),old=rmfield(old,ignore{k});end
    if isfield(new,ignore{k}),new=rmfield(new,ignore{k});end
end
assert(isequaln(old,new),'VGPP:MigrationParameters','Numerical configuration differs.');
row=c.row;particleTable=c.particleTable;
it=find(abs(p.targetVFSweepList-row.TargetVF)<1e-12);
assert(height(row)==1 && isscalar(it) && row.Run>=1 && ...
    row.Run<=p.targetVFSweepRuns && row.Run==fix(row.Run) && ...
    row.Seed==7300+100*it+row.Run && ...
    row.TargetN==round(row.TargetVF*p.Lx*p.Ly/(pi*p.rMono^2)) && ...
    ismember(row.Method,["VGPP","RSA-PBC"]), ...
    'VGPP:MigrationIdentity','Invalid case identity.');
if row.Method=="VGPP" && any(particleTable.x==0 | particleTable.y==0)
    status="rerun_boundary_information_lost";return;
end
[bad,~,~,~]=final_overlap_audit_PBC(particleTable.x,particleTable.y,particleTable.radius,p);
assert(~bad && all(abs(particleTable.radius-p.rMono)<=1e-12), ...
    'VGPP:MigrationGeometry','Invalid saved geometry.');
historicalParticleTable=particleTable;
validate_historical_coordinates(particleTable,historicalParticleTable,p);
[low,rsd]=local_vf_quality_metrics(historicalParticleTable,p);
assert(abs(low-row.LowVFCellFraction)<=p.validationTolerance && ...
    abs(rsd-row.LocalVFRSD)<=p.validationTolerance, ...
    'VGPP:MigrationStatistics','Saved local statistics disagree.');
row.CanonicalLocalVFRSD=rsd;row.CanonicalLowVFCellFraction=low;
row.LocalVFRSD=rsd;row.LowVFCellFraction=low;
configuration=p;sourceFingerprint=p.sourceFingerprint;
formatVersion=3;statisticalConvention='historical-boundary-v1';
migration=struct('originalSourceFingerprint',known, ...
    'originalFormatVersion',2,'elapsedTimeOrigin','preserved original generation timer');
if ~isfolder(p.checkpointFolder),mkdir(p.checkpointFolder);end
[~,name,ext]=fileparts(oldFile);target=fullfile(p.checkpointFolder,[name ext]);
if isfile(target),status="existing_v3_preserved";return;end
% Existing files are never overwritten, including failed migration attempts.
tmp=[tempname(p.checkpointFolder) '.mat'];
cleanup=onCleanup(@() removeTemporary(tmp));
save(tmp,'row','particleTable','historicalParticleTable','configuration', ...
    'sourceFingerprint','formatVersion','statisticalConvention','migration','-v7.3');
[ok,msg]=movefile(tmp,target);assert(ok,'VGPP:MigrationWrite','%s',msg);
status="migrated_losslessly";
end

function removeTemporary(file)
if isfile(file),delete(file);end
end
