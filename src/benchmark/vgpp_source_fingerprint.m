function fingerprint = vgpp_source_fingerprint(repoRoot)
%VGPP_SOURCE_FINGERPRINT SHA-256 of the active MATLAB benchmark source tree.
% MATLAB R2024b runs with JVM support; this checksum is independent of paths.
sourceDirs = {'config','src/core','src/analysis','src/benchmark'};
relativeFiles = {'main.m'};
for k = 1:numel(sourceDirs)
    listing = dir(fullfile(repoRoot,sourceDirs{k},'*.m'));
    for j = 1:numel(listing)
        relativeFiles{end+1} = [strrep(sourceDirs{k},'\','/') '/' listing(j).name]; %#ok<AGROW>
    end
end
relativeFiles = sort(relativeFiles);
md = java.security.MessageDigest.getInstance('SHA-256');
for k = 1:numel(relativeFiles)
    name = relativeFiles{k};
    fileName = fullfile(repoRoot,strrep(name,'/',filesep));
    fid = fopen(fileName,'rb');
    if fid < 0
        error('VGPP:SourceFingerprint','Cannot read source file: %s',fileName);
    end
    guard = onCleanup(@() fclose(fid));
    bytes = fread(fid,Inf,'*uint8');
    clear guard;
    % Include relative names and byte lengths to delimit concatenated files.
    md.update(typecast(uint8(name),'int8'));
    md.update(int8(0));
    md.update(typecast(uint8(sprintf('%d',numel(bytes))),'int8'));
    md.update(int8(0));
    if ~isempty(bytes)
        md.update(typecast(bytes,'int8'));
    end
end
raw = typecast(md.digest(),'uint8');
fingerprint = lower(reshape(dec2hex(raw,2).',1,[]));
end
