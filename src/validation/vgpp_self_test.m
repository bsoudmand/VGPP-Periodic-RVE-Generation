function vgpp_self_test(p)
%VGPP_SELF_TEST Check geometry and weighted-sampling conventions.
d = periodic_center_distance(0.01,0.02,9.99,9.98,p.Lx,p.Ly);
assert(abs(d-hypot(0.02,0.04)) < 1e-12,'Periodic distance failed.');
nTarget = round(0.7*p.Lx*p.Ly/(pi*p.rMono^2));
assert(nTarget==688,'Target-count convention differs from baseline.');
assert(p.overlapAuditTolerance==1e-8,'Acceptance tolerance changed.');
assert(~p.useHexagonalFallback,'Ordered fallback is not part of the benchmark.');
if isfield(p,'sourceFingerprint')
    assert(numel(char(p.sourceFingerprint))==64,'Invalid source checksum length.');
end
small = table((1:2)',[0.01;9.99],[0.5;0.5],[0.18;0.18], ...
    'VariableNames',{'ID','x','y','radius'});
[invalid,overlap,~,nBad] = final_overlap_audit_PBC( ...
    small.x,small.y,small.radius,p);
assert(invalid && overlap>p.overlapAuditTolerance && nBad==1, ...
    'Periodic seam overlap was not detected.');
small.x(2)=5.0;
[invalid,~,~,nBad] = final_overlap_audit_PBC( ...
    small.x,small.y,small.radius,p);
assert(~invalid && nBad==0,'Nonoverlapping pair was rejected.');
% Nonfinite, out-of-domain, and nonpositive radii must never pass auditing.
assertInvalidGeometry(@() final_overlap_audit_PBC([NaN;5],[0.5;0.5],[0.18;0.18],p));
assertInvalidGeometry(@() final_overlap_audit_PBC([Inf;5],[0.5;0.5],[0.18;0.18],p));
assertInvalidGeometry(@() final_overlap_audit_PBC([0.1;p.Lx],[0.5;0.5],[0.18;0.18],p));
assertInvalidGeometry(@() final_overlap_audit_PBC([0.1;5],[0.5;0.5],[0.18;0],p));

oldRng = rng; cleanRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(29);
w = ones(50000,1);
nSelected = 4000;
idx = weighted_sample_without_replacement(w,nSelected,'log_keys');
assert(numel(idx)==nSelected && numel(unique(idx))==nSelected, ...
    'Stable weighted sampler returned duplicate indices.');
assert(median(idx)>20000 && median(idx)<30000, ...
    'Stable equal-weight sampling appears biased toward early indices.');
rng(29);
idxLegacy=weighted_sample_without_replacement(w,nSelected,'legacy');
assert(numel(idxLegacy)==nSelected && numel(unique(idxLegacy))==nSelected, ...
    'Legacy sampler returned duplicate indices.');
% Legacy selection is checked for uniqueness only; its distribution can be biased.
assert(all(idx>=1 & idx<=numel(w)) && all(idxLegacy>=1 & idxLegacy<=numel(w)), ...
    'Weighted sampler returned out-of-range indices.');
end

function assertInvalidGeometry(fun)
rejected = false;
try
    fun();
catch ME
    rejected = startsWith(string(ME.identifier),"MATLAB:") || ...
        startsWith(string(ME.identifier),"VGPP:");
end
assert(rejected,'Malformed particle geometry escaped the audit.');
end
