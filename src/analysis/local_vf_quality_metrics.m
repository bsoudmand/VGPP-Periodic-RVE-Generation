function [lowVFCellFraction, localVFRSD] = local_vf_quality_metrics(particleTable, params)
%LOCAL_VF_QUALITY_METRICS Cell-wise variation and depleted-cell fraction.
nBins = params.matrixRichBins;
[localVFMap, localVFRSD] = local_volume_fraction_map_by_centers( ...
    particleTable.x, particleTable.y, particleTable.radius, ...
    params.Lx, params.Ly, nBins);
vals = localVFMap(:);
threshold = params.matrixRichLowVFThresholdFactor * mean(vals);
lowVFCellFraction = mean(vals < threshold);
end
