function stats = compute_spatial_statistics_PBC(particleTable, params, methodName)

    x = particleTable.x(:);
    y = particleTable.y(:);
    r = particleTable.radius(:);

    Lx = params.Lx;
    Ly = params.Ly;
    A  = Lx * Ly;
    N  = numel(r);

    %% All pair distances and pair indices under PBC
    [pairD, pairI, pairJ] = all_pair_distances_PBC(x, y, Lx, Ly);

    %% Nearest-neighbor distance and orientation
    nearestD = zeros(N,1);
    nearestRequired = zeros(N,1);
    nearestAngle = zeros(N,1);

    for i = 1:N
        dAll = periodic_center_distance(x(i), y(i), x, y, Lx, Ly);
        dAll(i) = inf;
        [nearestD(i), jMin] = min(dAll);
        nearestRequired(i) = r(i) + r(jMin) + params.minGap;

        dx = x(jMin) - x(i);
        dy = y(jMin) - y(i);
        dx = dx - Lx * round(dx / Lx);
        dy = dy - Ly * round(dy / Ly);

        angle = atan2d(dy, dx);
        angle = mod(angle, 180);   % undirected nearest-neighbor orientation
        nearestAngle(i) = angle;
    end

    normalizedNND = nearestD ./ nearestRequired;

    angleGrid = linspace(0, 180, 181)';
    orientationCDF = zeros(size(angleGrid));
    for k = 1:numel(angleGrid)
        orientationCDF(k) = mean(nearestAngle <= angleGrid(k));
    end
    csrCDF = angleGrid / 180;
    orientationKS = max(abs(orientationCDF - csrCDF));

    %% Surface-to-surface clearances for all pairs
    requiredPairD = r(pairI) + r(pairJ) + params.minGap;
    pairClearance = pairD - requiredPairD;

    %% Ripley K and L functions under PBC
    rMaxK = 0.50 * min(Lx, Ly);
    rK = linspace(max(mean(r), eps), rMaxK, 80)';
    K = zeros(size(rK));

    for k = 1:numel(rK)
        % Ordered-pair count equals twice the unordered count.
        orderedCount = 2 * sum(pairD <= rK(k));
        K(k) = A * orderedCount / (N^2);
    end

    Lfun = sqrt(K / pi);
    LminusR = Lfun - rK;
    meanAbsLminusR = mean(abs(LminusR), 'omitnan');

    %% Radial distribution function g(r)
    nBins = 70;
    edges = linspace(0, rMaxK, nBins + 1);
    counts = histcounts(pairD, edges);
    rRDF = 0.5 * (edges(1:end-1) + edges(2:end))';
    shellArea = pi * (edges(2:end).^2 - edges(1:end-1).^2)';
    numberDensity = N / A;

    % Expected unordered-pair count in an annulus for CSR under PBC.
    expectedCounts = 0.5 * N * numberDensity .* shellArea;
    rdf = counts(:) ./ max(expectedCounts, eps);

    contactDistance = mean(2*r + params.minGap);
    firstPeakRegion = rRDF >= 0.90*contactDistance & rRDF <= 2.50*contactDistance;
    if any(firstPeakRegion)
        rdfFirstPeak = max(rdf(firstPeakRegion));
    else
        rdfFirstPeak = max(rdf);
    end

    %% Quadrat/local volume-fraction uniformity map
    nBinsXY = 10;
    [localVFMap, localVFRSD] = local_volume_fraction_map_by_centers( ...
        x, y, r, Lx, Ly, nBinsXY);

    %% Final audit
    [~, maxOverlap, minClearance, nOverlappingPairs, worstPair] = ...
        final_overlap_audit_PBC(x, y, r, params);

    %% Collect output
    stats.methodName = methodName;
    stats.N = N;
    stats.volumeFraction = sum(pi*r.^2) / A;
    stats.meanRadius = mean(r);

    stats.nearestD = nearestD;
    stats.normalizedNND = normalizedNND;
    stats.meanNND = mean(nearestD);
    stats.nndCV = std(nearestD) / mean(nearestD);

    stats.orientationAngleDeg = angleGrid;
    stats.orientationCDF = orientationCDF;
    stats.orientationKS = orientationKS;

    stats.pairClearance = pairClearance;
    stats.minClearance = minClearance;
    stats.maxOverlap = maxOverlap;
    stats.nOverlappingPairs = nOverlappingPairs;
    stats.worstPair = worstPair;

    stats.rK = rK;
    stats.K = K;
    stats.LminusR = LminusR;
    stats.meanAbsLminusR = meanAbsLminusR;

    stats.rRDF = rRDF;
    stats.rdf = rdf;
    stats.rdfFirstPeak = rdfFirstPeak;

    stats.localVFMap = localVFMap;
    stats.localVFRSD = localVFRSD;
end
