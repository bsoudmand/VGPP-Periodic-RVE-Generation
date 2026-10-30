function [hasOverlap, maxOverlap, minClearance, nOverlappingPairs, worstPair] = ...
    final_overlap_audit_PBC(cxList, cyList, rList, params)

    % Explicit all-pair overlap audit under periodic boundary conditions.
    %
    % This function is the final authority for geometric validity. It checks
    % every particle pair using the minimum-image convention, so it detects
    % both ordinary in-domain overlaps and mirrored overlaps across opposite
    % periodic edges.

    % Reject noncanonical or nonfinite coordinates before pairwise checks.
    validateattributes(params.Lx, {'numeric'}, {'real','scalar','finite','positive'});
    validateattributes(params.Ly, {'numeric'}, {'real','scalar','finite','positive'});
    validateattributes(params.minGap, {'numeric'}, {'real','scalar','finite','nonnegative'});
    validateattributes(cxList, {'numeric'}, {'real','vector','finite'});
    validateattributes(cyList, {'numeric'}, {'real','vector','finite'});
    validateattributes(rList, {'numeric'}, {'real','vector','finite'});
    if numel(cxList) ~= numel(cyList) || numel(cxList) ~= numel(rList)
        error('VGPP:GeometryShape', 'Particle center and radius arrays must have equal lengths.');
    end
    if any(rList <= 0)
        error('VGPP:InvalidRadius', 'Particle radii must be strictly positive.');
    end
    if any(cxList < 0 | cxList >= params.Lx) || ...
            any(cyList < 0 | cyList >= params.Ly)
        error('VGPP:InvalidCoordinates', ...
            'Particle centers must use the periodic convention 0 <= x < Lx and 0 <= y < Ly.');
    end
    Lx = params.Lx;
    Ly = params.Ly;
    minGap = params.minGap;
    if isfield(params,'overlapAuditTolerance')
        auditTol = params.overlapAuditTolerance;
    else
        auditTol = 1.0e-8;
    end
    validateattributes(auditTol, {'numeric'}, {'real','scalar','finite','nonnegative'});

    nParticles = numel(rList);

    hasOverlap = false;
    maxOverlap = 0;
    minClearance = inf;
    nOverlappingPairs = 0;
    worstPair = [NaN, NaN];

    if nParticles < 2
        minClearance = NaN;
        return;
    end

    for i = 1:nParticles-1

        xNew = cxList(i);
        yNew = cyList(i);

        xOthers = cxList(i+1:end);
        yOthers = cyList(i+1:end);
        rOthers = rList(i+1:end);

        d = periodic_center_distance(xNew, yNew, xOthers, yOthers, Lx, Ly);

        requiredDistance = rList(i) + rOthers + minGap;
        clearance = d - requiredDistance;
        overlap = requiredDistance - d;

        minClearance = min(minClearance, min(clearance));

        overlappingLocal = overlap > auditTol;
        nOverlappingPairs = nOverlappingPairs + sum(overlappingLocal);

        [localMaxOverlap, localIdx] = max([0; overlap(:)]);

        if localMaxOverlap > maxOverlap
            maxOverlap = localMaxOverlap;

            if localIdx > 1
                worstPair = [i, i + localIdx - 1];
            end
        end
    end

    hasOverlap = nOverlappingPairs > 0;

    if maxOverlap <= auditTol
        maxOverlap = 0;
    end

    if isinf(minClearance)
        minClearance = NaN;
    end
end

