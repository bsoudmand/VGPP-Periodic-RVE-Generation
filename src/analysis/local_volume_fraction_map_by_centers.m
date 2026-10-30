function [localVFMap, localVFRSD] = local_volume_fraction_map_by_centers( ...
    x, y, r, Lx, Ly, nBinsXY)

    xEdges = linspace(0, Lx, nBinsXY + 1);
    yEdges = linspace(0, Ly, nBinsXY + 1);

    localVFMap = zeros(nBinsXY, nBinsXY);
    binArea = (Lx / nBinsXY) * (Ly / nBinsXY);

    for i = 1:numel(r)
        ix = find(x(i) >= xEdges(1:end-1) & x(i) < xEdges(2:end), 1, 'first');
        iy = find(y(i) >= yEdges(1:end-1) & y(i) < yEdges(2:end), 1, 'first');

        if isempty(ix)
            ix = nBinsXY;
        end
        if isempty(iy)
            iy = nBinsXY;
        end

        localVFMap(iy, ix) = localVFMap(iy, ix) + pi*r(i)^2 / binArea;
    end

    localVFValues = localVFMap(:);
    localVFRSD = std(localVFValues) / max(mean(localVFValues), eps);
end
