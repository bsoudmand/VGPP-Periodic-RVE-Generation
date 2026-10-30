function [pairD, pairI, pairJ] = all_pair_distances_PBC(x, y, Lx, Ly)

    N = numel(x);
    nPairs = N * (N - 1) / 2;

    pairD = zeros(nPairs,1);
    pairI = zeros(nPairs,1);
    pairJ = zeros(nPairs,1);

    counter = 0;

    for i = 1:N-1
        jList = (i+1:N)';
        dx = x(jList) - x(i);
        dy = y(jList) - y(i);

        dx = dx - Lx * round(dx / Lx);
        dy = dy - Ly * round(dy / Ly);

        d = sqrt(dx.^2 + dy.^2);
        m = numel(jList);

        idx = counter + (1:m);
        pairD(idx) = d;
        pairI(idx) = i;
        pairJ(idx) = jList;

        counter = counter + m;
    end
end
