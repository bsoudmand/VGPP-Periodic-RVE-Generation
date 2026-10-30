function report = verify_reference_results(actualFile, referenceFile, outDir, tol)
%VERIFY_REFERENCE_RESULTS Compare rerun records against archived benchmark.
assert(isfile(actualFile) && isfile(referenceFile), 'Benchmark or reference CSV missing.');
a = readtable(actualFile,'TextType','string');
b = readtable(referenceFile,'TextType','string');
keys = {'Method','TargetVF','Run','Seed'};
required = [keys,{'TargetN','ReachedTarget','NumberOfParticles','FinalVF','MeanRadius','MaxOverlap', ...
    'MinClearance','OverlappingPairs','DiscreteTargetVF','RadiusScale', ...
    'LocalVFRSD','LowVFCellFraction','MeanNND','NNDCV','RDFFirstPeak'}];
assert(all(ismember(required,a.Properties.VariableNames)) && ...
    all(ismember(required,b.Properties.VariableNames)), 'Missing validation columns.');
keyA = string(a.Method)+"|"+compose('%.10f',a.TargetVF)+"|"+string(a.Run)+"|"+string(a.Seed);
keyB = string(b.Method)+"|"+compose('%.10f',b.TargetVF)+"|"+string(b.Run)+"|"+string(b.Seed);
assert(numel(unique(keyA))==height(a) && numel(unique(keyB))==height(b), ...
    'Duplicate comparison identifiers found.');
[found,loc] = ismember(keyB,keyA);
assert(all(found) && height(a)==height(b), 'Run identifiers differ from reference.');
if ismember('GeometryValid',a.Properties.VariableNames)
    assert(all(a.GeometryValid),'Invalid geometries found in rerun benchmark.');
end
a=a(loc,:);
exactFields={'TargetN','ReachedTarget','NumberOfParticles','OverlappingPairs'};
metricFields={'FinalVF','DiscreteTargetVF','MeanRadius','RadiusScale', ...
    'MaxOverlap','MinClearance','LocalVFRSD','LowVFCellFraction', ...
    'MeanNND','NNDCV','RDFFirstPeak'};
fields=[exactFields,metricFields];
Field=string(fields(:)); MaxAbsoluteDifference=NaN(numel(fields),1);
Pass=false(numel(fields),1);
NonfinitePatternMatch=false(numel(fields),1);
for k=1:numel(fields)
    x=double(a.(fields{k})); y=double(b.(fields{k}));
    masksMatch = isequal(isfinite(x),isfinite(y)) && ...
        isequal(isnan(x),isnan(y)) && ...
        isequal(isinf(x),isinf(y)) && ...
        isequal(x(isinf(x)),y(isinf(y)));
    NonfinitePatternMatch(k) = masksMatch;
    finite = isfinite(x) & isfinite(y);
    if any(finite)
        MaxAbsoluteDifference(k) = max(abs(x(finite)-y(finite)));
    elseif masksMatch
        MaxAbsoluteDifference(k) = 0;
    end
    threshold = tol;
    if ismember(fields{k},exactFields), threshold=0; end
    Pass(k) = masksMatch && MaxAbsoluteDifference(k)<=threshold;
end
report=table(Field,MaxAbsoluteDifference,NonfinitePatternMatch,Pass);
writetable(report,fullfile(outDir,'reference_comparison.csv'));
if ~all(Pass)
    warning('VGPP:ReferenceMismatch', ...
        'Benchmark differs from archived results. Inspect reference_comparison.csv.');
else
    fprintf('Archived benchmark numeric comparisons passed (excluding elapsed time).\n');
end
end
